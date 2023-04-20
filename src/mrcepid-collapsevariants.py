#!/usr/bin/env python
#
# Author: Eugene Gardner (eugene.gardner at mrc.epid.cam.ac.uk)
#
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/

import os
import sys
import dxpy
import tarfile
import pandas as pd

from dxpy import DXFile
from pathlib import Path
from typing import Tuple, Dict, List

from general_utilities.association_resources import generate_linked_dx_file, run_cmd
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.mrc_logger import MRCLogger

# We have to do this to get modules to run properly on DNANexus while still enabling easy editing in PyCharm
sys.path.append('/')
sys.path.append('/collapsevariants/')
sys.path.append('/collapsevariants/tool_parsers/')

# DO NOT move this. It MUST come after the above 'sys.path.append' code to make sure packages run properly
from collapsevariants.ingest_data import IngestData
from collapsevariants.snp_list_generator import SNPListGenerator
from collapsevariants.filtering import filter_bgen
from collapsevariants.snp_list_merger import SNPMerger
from collapsevariants.collapse_logger import CollapseLOGGER

# Set up the system logger – this is NOT the same as LOG_FILE below that records info about the filtering itself
LOGGER = MRCLogger().get_logger()


def stat_writer(sample_tables: List[pd.DataFrame], per_chromosome_totals: Dict[str, int], LOG_FILE: CollapseLOGGER,
                total_sites: int) -> None:
    """ Writes stats about the various collapsing operations performed by this applet

    :param sample_tables: A list of Pandas dataframes (one for each chromosome queried) containing per sample and per
        gene information.
    :param per_chromosome_totals: Total number of variants found per-chromosome
    :param LOG_FILE: The LOG_FILE for this instance to print to
    :param total_sites: Total number of expected sites based on the original query
    :return: None
    """

    # Write a bunch of stats
    LOG_FILE.write_header('Per-chromosome totals')

    found_total_sites = 0
    for chromosome, total in per_chromosome_totals.items():
        found_total_sites += total
        LOG_FILE.write_int(f'Total sites on chr{chromosome}', total)

    LOG_FILE.write_header('Genome-wide totals')
    LOG_FILE.write_int('Total sites expected from filtering expression', total_sites)
    LOG_FILE.write_int('Total sites extracted from all chromosomes', found_total_sites)
    LOG_FILE.write_string('Total expected and total extracted match', str(total_sites == found_total_sites))
    LOG_FILE.write_spacer()

    # Concatenate and sum sample tables:
    master_sample_table = pd.concat(sample_tables).groupby(['sample_id']).sum().reset_index()
    LOG_FILE.write_header('Per-individual totals')
    LOG_FILE.write_int('Median number of alleles per indv', master_sample_table['ac'].median())
    LOG_FILE.write_int('Median number of genes affected per indv', master_sample_table['ac_gene'].median())
    LOG_FILE.write_float('Mean number of alleles per indv', master_sample_table['ac'].mean())
    LOG_FILE.write_float('Mean number of genes affected per indv', master_sample_table['ac_gene'].mean())
    LOG_FILE.write_int('Max number of alleles', master_sample_table['ac'].max())
    LOG_FILE.write_int('Number of individuals with at least 1 allele',
                       pd.value_counts(master_sample_table['ac'] > 0)[True])
    LOG_FILE.write_spacer()

    LOG_FILE.write_header('AC Histogram')
    LOG_FILE.write_generic('AC_bin\tcount\n')
    ac_counts = master_sample_table.value_counts('ac')
    ac_counts = ac_counts.sort_index()
    for ac, count in ac_counts.items():
        LOG_FILE.write_histogram(ac, count)


@dxpy.entry_point('main')
def main(filtering_expression: str, snplist: dict, genelist: dict, output_prefix: str,
         bgen_index: dict, testing_script: dict, testing_directory: str) -> Dict[str, dict]:
    """The main entrypoint in the DNANexus applet that runs CollapseVariants

    This method collapses variants using a variety of ways in four steps (see individual classes / README for more
    information):

    1. Ingest required data onto the AWS instance (class IngestData)
    2. Generate a list of variants to collapse on (class SNPListGenerator)
    3. Filter and format these variants into files appropriate for various burden testing approaches
    4. Write statistics about the collapsing that has been performed

    :param filtering_expression: A string filtering expression to filter variants (must be compatible with
        pandas.query())
    :param snplist: A DXFile containing a list of varIDs to use as a custom mask
    :param genelist: A DXFile containing a list of gene symbols to collapse into a custom mask
    :param output_prefix: A name to append to beginning of output files.
    :param bgen_index: A DXFile containing information of bgen files to collapse on
    :param testing_script: Script compatible with pytest. If not null, invoke the collapsevariants testing suite
        via :func:`test`.
    :param testing_directory: Directory containing test files if in testing mode.
    :return: A Dictionary with keys of output strings and values of DXIndex
    """

    if testing_script is not None:
        LOGGER.info('Testing mode activated...')
        if testing_directory is None:
            raise ValueError(f'Testing mode invoked but -itesting_directory not provided!')
        output_tarball, linked_log_file = test(output_prefix, bgen_index, testing_script, testing_directory)
        output = {'output_tarball': dxpy.dxlink(generate_linked_dx_file(output_tarball)),
                  'log_file': dxpy.dxlink(linked_log_file)}

    else:
        # Set up our logfile for recording information on
        LOG_FILE = CollapseLOGGER(output_prefix)

        # This loads all data
        ingested_data = IngestData(bgen_index, filtering_expression, snplist, genelist)

        # First generate a list of ALL variants genome-wide that we want to retain:
        snp_list_generator = SNPListGenerator(ingested_data,
                                              LOG_FILE)

        # Now build a thread worker that contains as many threads
        # instance takes a thread and 1 thread for monitoring
        thread_utility = ThreadUtility()

        # Now loop through each chromosome and do the actual filtering...
        # ...launch the requested threads
        for chromosome in snp_list_generator.chromosomes:
            chrom_bgen_index = ingested_data.bgen_index[chromosome]
            thread_utility.launch_job(filter_bgen,
                                      file_prefix=output_prefix,
                                      chromosome=chromosome,
                                      chrom_bgen_index=chrom_bgen_index)

        # And gather the resulting futures
        sample_tables = []
        chromosome_totals = dict()

        for result in thread_utility:
            per_chromosome_total, chromosome, sample_table = result
            sample_tables.append(sample_table)
            chromosome_totals[chromosome] = per_chromosome_total
        LOG_FILE.write_spacer()

        # And write stats about the collected futures
        stat_writer(sample_tables, chromosome_totals, LOG_FILE, snp_list_generator.total_sites)

        # Here we check if we made a SNP-list. If so, we need to merge across all chromosomes into single
        # per-snp/gene files:
        if ingested_data.found_snps or ingested_data.found_genes:
            LOGGER.info('Making merged SNP files for burden testing...')
            SNPMerger(snp_list_generator.chromosomes, output_prefix, ingested_data.found_genes)

        LOGGER.info('Closing LOG file...')
        linked_log_file = LOG_FILE.close_writer()

        # Here we are taking all the files generated by the various functions above and adding them to a single tar
        # to enable easy output. The only output of this applet is thus a single .tar.gz file per VCF file
        LOGGER.info('Generating final tarball...')
        output_tarball = Path(f'{output_prefix}.tar.gz')
        tar = tarfile.open(output_tarball, 'w:gz')
        for file in Path('./').glob(f'{output_prefix}.*'):
            if '.tar.gz' not in file.name:  # Don't want to remove the tar itself... yet...
                tar.add(file)
                file.unlink()
        tar.close()

        # Set output
        output = {'output_tarball': dxpy.dxlink(generate_linked_dx_file(output_tarball)),
                  'log_file': dxpy.dxlink(linked_log_file)}

    return output


def test(output_prefix: str, bgen_index: dict, testing_script: dict, testing_directory: str) -> Tuple[Path, DXFile]:
    """Run the collapsevariants testing suite.

    This method is invisible to the applet and can only be accessed by using API calls via dxpy.DXApplet() on
    a local machine. See the resources in the `./test/` folder for more information on running tests.

    :param output_prefix: A prefix to name the output tarball returned by this method.
    :param bgen_index: A DXFile containing information of bgen files to collapse on
    :param testing_script: The dxfile ID of the pytest-compatible script
    :param testing_directory: The name of the folder containing test resources on the DNANexus platform
    :return: Dict of containing the pytest log in a tar.gz to ensure compatibility with the main() method returns
    """

    LOGGER.info('Launching mrcepid-collapsevariants with the testing suite')
    dxpy.download_dxfile(dxid=testing_script['$dnanexus_link'], filename='test.py')

    # I then set an environment variable that tells pytest where the testing directory is
    os.environ['CI'] = '500'  # Make sure pytest logs aren't truncated
    os.environ['TEST_DIR'] = testing_directory
    LOGGER.info(f'TEST_DIR environment variable set: {os.getenv("TEST_DIR")}')
    os.environ['BGEN_INDEX'] = bgen_index['$dnanexus_link']
    LOGGER.info(f'BGEN_INDEX environment variable set: {os.getenv("BGEN_INDEX")}')

    # pytest always throws an error when a test fails, which causes the entire suite to fall apart (and,
    # problematically, not return the logfile...). So we catch a runtime error if thrown by run_cmd() and then return
    # the log that (hopefully) should already exist. This will fall apart if there is an issue with run_cmd that is
    # outside of running pytest.
    out_log = Path(f'pytest.{output_prefix}.log')
    try:
        run_cmd('pytest test.py', is_docker=False, stdout_file=out_log)
    except RuntimeError:
        pass

    output_tarball = Path(f'{output_prefix}.tar.gz')
    tar = tarfile.open(output_tarball, "w:gz")
    tar.add(out_log)
    tar.close()

    # Write a fake logger as it is a required output
    fake_logger = CollapseLOGGER('fake_logger')
    fake_logger.write_header('FAKE')
    fake_logger_file = fake_logger.close_writer()

    return output_tarball, fake_logger_file


dxpy.run()
