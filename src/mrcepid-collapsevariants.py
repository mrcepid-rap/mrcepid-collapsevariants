#!/usr/bin/env python
#
# Author: Eugene Gardner (eugene.gardner at mrc.epid.cam.ac.uk)
#
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/

import sys
import dxpy
import tarfile
import pandas as pd

from pathlib import Path
from typing import Tuple, Dict, List

from general_utilities.association_resources import generate_linked_dx_file
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.mrc_logger import MRCLogger

# We have to do this to get modules to run properly on DNANexus while still enabling easy editing in PyCharm
sys.path.append('/')
sys.path.append('/collapsevariants/')
sys.path.append('/collapsevariants/tool_parsers/')

# DO NOT move this. It MUST come after the above 'sys.path.append' code to make sure packages run properly
from collapsevariants.ingest_data import IngestData
from collapsevariants.snp_list_generator import SNPListGenerator
from collapsevariants.snp_list_merger import SNPMerger
from collapsevariants.filtering import run_filtering
from collapsevariants.tool_parsers.bolt_parser import parse_filters_BOLT, check_vcf_stats
from collapsevariants.tool_parsers.saige_parser import parse_filters_SAIGE
from collapsevariants.tool_parsers.staar_parser import STAARParser, STAARMergingException
from collapsevariants.collapse_logger import CollapseLOGGER

# Set up the system logger – this is NOT the same as LOG_FILE below that records info about the filtering itself
LOGGER = MRCLogger().get_logger()


def filter_bgen(file_prefix: str, chromosome: str, chrom_bgen_index: dict) -> Tuple[int, str, pd.DataFrame]:
    """Helper method for running separate BGENs through the collapsing process (returns function as a future).

    This method just works through all possible formatting methods (for BOLT, SAIGE, and STAAR) to generate a final
    suite of files for input to the runassociationtesting burden testing suite. These files are:

    1. A BOLT-ready .bgen and .sample file
    2. A SAIGE/REGENIE-ready .bcf, bcf.csi (index), and groupFile
    3. A STAAR-ready sparse matrix saved in .rds format and a .tsv containing variant information

    :param file_prefix: A name to append to beginning of output files.
    :param chromosome: The chromosome currently being processed. This must be the short form of the chromosome name
        (e.g., '1' not 'chr1').
    :param chrom_bgen_index: a BGENIndex TypedDict containing information on filepaths containing filtered & annotated
        variants.
    :return: A Tuple containing the total number of variants found for this chromosome after filtering, the chromosome
        ID, and a pandas.DataFrame containing per-sample and per-ENST totals for log reporting purposes.
    """

    # Ingest the filtered BGEN file into this instance
    LOGGER.info(f'Processing bgen: {chromosome}.filtered.bgen')
    bgenprefix = f'filtered_bgen/chr{chromosome}.filtered'  # Get a prefix name for all files

    # Download the requisite files for this chromosome according to the index dict:
    dxpy.download_dxfile(dxpy.DXFile(chrom_bgen_index['index']).get_id(), f'{bgenprefix}.bgen.bgi')
    dxpy.download_dxfile(dxpy.DXFile(chrom_bgen_index['sample']).get_id(), f'{bgenprefix}.sample')
    dxpy.download_dxfile(dxpy.DXFile(chrom_bgen_index['bgen']).get_id(), f'{bgenprefix}.bgen')

    # Run filtering according to the user-provided filtering expression
    # Also gets the total number of variants retained for this chromosome
    num_variants = run_filtering(bgenprefix, chromosome, file_prefix)
    LOGGER.info(f'Identified {num_variants} variants that match the given filtering expression in file '
                f'{file_prefix}.{chromosome}.bgen')

    # If there are no variants we need to create an empty dummy pd.DataFrame to signify that this process completed BUT
    # that no variants were found
    if num_variants == 0:
        LOGGER.warning(f'No files found for chromosome {chromosome}, excluding from final files...')
        return num_variants, chromosome, pd.DataFrame()
    else:
        # Here we are then taking the file generated by run_filtering() and generating various text/plink/vcf files
        # to generate a merged set of variants we want to test across all VCF files and that will be used as part of
        # mrcepid-mergecollapsevariants.
        #
        # JUST TO BE CLEAR – the names of the functions here are not THAT important (e.g., files generated in the
        # function parse_filters_BOLT() will be used for other tools/workflows). It was just for me (Eugene Gardner)
        # to keep things organised when writing this code
        genes, snp_gene_map = parse_filters_SAIGE(file_prefix, chromosome)

        poss_indv, samples = parse_filters_BOLT(file_prefix, chromosome, genes, snp_gene_map)
        sample_table = check_vcf_stats(poss_indv, samples)

        # STAAR fails sometimes for unknown reasons, so try it twice if it fails before throwing the entire process
        try:
            STAARParser(file_prefix, chromosome).parse_filters_STAAR()
        except STAARMergingException:
            LOGGER.warning(f'STAAR chr {chromosome} failed to merge, trying again...')
            STAARParser(file_prefix, chromosome).parse_filters_STAAR()

        # Purge files that we no longer need:
        Path(f'{file_prefix}.{chromosome}.bgen').unlink()
        Path(f'{file_prefix}.{chromosome}.bgen.bgi').unlink()
        Path(f'{file_prefix}.{chromosome}.parsed.txt').unlink()
        Path(f'{file_prefix}.{chromosome}.snps').unlink()

        LOGGER.info(f'Finished bgen: chr{chromosome}.filtered.bgen')
        return num_variants, chromosome, sample_table


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
def main(filtering_expression: str, snplist: dict, genelist: dict, file_prefix: str,
         bgen_index: dict) -> Dict[str, dict]:
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
    :param file_prefix: A name to append to beginning of output files.
    :param bgen_index: A DXFile containing information of bgen files to collapse on
    :return: A Dictionary with keys of output strings and values of DXIndex
    """

    # Set up our logfile for recording information on
    LOG_FILE = CollapseLOGGER(file_prefix)

    # This loads all data
    ingested_data = IngestData(filtering_expression, bgen_index, snplist, genelist)

    # First generate a list of ALL variants genome-wide that we want to retain:
    snp_list_generator = SNPListGenerator(ingested_data.bgen_index,
                                          ingested_data.filtering_expression,
                                          ingested_data.found_snps,
                                          ingested_data.found_genes,
                                          LOG_FILE)

    # Now build a thread worker that contains as many threads
    # instance takes a thread and 1 thread for monitoring
    thread_utility = ThreadUtility()

    # Now loop through each chromosome and do the actual filtering...
    # ...launch the requested threads
    for chromosome in snp_list_generator.chromosomes:
        chrom_bgen_index = ingested_data.bgen_index[chromosome]
        thread_utility.launch_job(filter_bgen,
                                  file_prefix=file_prefix,
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

    # Here we check if we made a SNP-list. If so, we need to merge across all chromosomes into single per-snp files:
    if ingested_data.found_snps or ingested_data.found_genes:  # run for gene list as well, add
        LOGGER.info('Making merged SNP files for burden testing...')
        SNPMerger(snp_list_generator.chromosomes, file_prefix, ingested_data.found_genes)

    LOGGER.info('Closing LOG file...')
    linked_log_file = LOG_FILE.close_writer()

    # Here we are taking all the files generated by the various functions above and adding them to a single tar
    # to enable easy output. The only output of this applet is thus a single .tar.gz file per VCF file
    LOGGER.info('Generating final tarball...')
    output_tarball = Path(f'{file_prefix}.tar.gz')
    tar = tarfile.open(output_tarball, 'w:gz')
    for file in Path('./').glob(f'{file_prefix}.*'):
        if '.tar.gz' not in file.name:  # Don't want to remove the tar itself... yet...
            tar.add(file)
            file.unlink()
    tar.close()

    # Set output
    output = {'output_tarball': dxpy.dxlink(generate_linked_dx_file(output_tarball)),
              'log_file': dxpy.dxlink(linked_log_file)}

    return output


dxpy.run()
