#!/usr/bin/env python
#
# Author: Eugene Gardner (eugene.gardner at mrc.epid.cam.ac.uk)
#
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/

import dxpy
import tarfile
import subprocess

from pathlib import Path
from typing import Dict

from collapsevariants.parallelization_wrappers import generate_generic_masks, generate_snp_or_gene_masks, \
    generate_genotype_matrices, update_log_file
from general_utilities.association_resources import generate_linked_dx_file
from general_utilities.mrc_logger import MRCLogger

from collapsevariants.ingest_data import IngestData
from collapsevariants.snp_list_generator import SNPListGenerator
from collapsevariants.collapse_logger import CollapseLOGGER

# Set up the system logger – this is NOT the same as LOG_FILE below that records info about the filtering itself
LOGGER = MRCLogger().get_logger()


@dxpy.entry_point('main')
def main(filtering_expression: str, snplist: dict, genelist: dict, output_prefix: str,
         bgen_index: dict) -> Dict[str, dict]:
    """The main entrypoint in the DNANexus applet that runs CollapseVariants

    This method collapses variants using a variety of ways in four steps (see individual classes / README for more
    information):

    1. Ingest required data onto the AWS instance (class IngestData)
    2. Generate a list of variants to collapse on (class SNPListGenerator)
    3. Convert variants into genotype matrices for each BGEN file (method generate_genotype_matrices)
    4. Write statistics about the collapsing that has been performed (method update_log_file)
    5. Filter and format these variants into files appropriate for various burden testing approaches (methods
         generate_snp_or_gene_masks and generate_generic_masks)

    :param filtering_expression: A string filtering expression to filter variants (must be compatible with
        pandas.query())
    :param snplist: A DXFile containing a list of varIDs to use as a custom mask
    :param genelist: A DXFile containing a list of gene symbols to collapse into a custom mask
    :param output_prefix: A name to append to beginning of output files.
    :param bgen_index: A DXFile containing information of bgen files to collapse on
    :return: A Dictionary with keys of output strings and values of DXIndex
    """

    # Set up our logfile for recording information on
    log_file = CollapseLOGGER(Path(f'{output_prefix}.log'))

    # Proceeding in the order of the steps above:
    # 1. Loads all DNANexus-specific data
    LOGGER.info('Ingesting data...')
    ingested_data = IngestData(bgen_index, filtering_expression, snplist, genelist)

    # 2. Generate a list of ALL variants genome-wide that we want to retain:
    LOGGER.info('Filtering variants according to provided inputs...')
    snp_list_generator = SNPListGenerator(ingested_data.vep_dict, ingested_data.filtering_expression,
                                          ingested_data.gene_list_path, ingested_data.snp_list_path, log_file)

    # 3. Generate genotype matrices for each bgen file:
    LOGGER.info('Generating genotype matrices for each bgen file...')
    genotype_index = generate_genotype_matrices(snp_list_generator.genes, ingested_data.bgen_index)

    # 4. Write information about collapsing to the log file
    LOGGER.info('Updating log file with per-sample and per-ENST totals...')
    update_log_file(snp_list_generator.genes, genotype_index, len(ingested_data.sample_ids),
                    snp_list_generator.total_sites, log_file)

    # 5. Filter the bgen files to only include the variants we want to keep
    LOGGER.info('Generating final output files for RunAssociationTesting...')
    if ingested_data.snp_list_path:
        LOGGER.info('Making SNP files for burden testing...')
        output_files = generate_snp_or_gene_masks(snp_list_generator.genes, genotype_index,
                                                  ingested_data.sample_ids, output_prefix, 'SNP')
    elif ingested_data.gene_list_path:
        LOGGER.info('Making GENE files for burden testing...')
        output_files = generate_snp_or_gene_masks(snp_list_generator.genes, genotype_index,
                                                  ingested_data.sample_ids, output_prefix, 'GENE')
    else:
        LOGGER.info('Making standard files for burden testing...')
        output_files = generate_generic_masks(snp_list_generator.genes, genotype_index, ingested_data.sample_ids,
                                              output_prefix)

    LOGGER.info('Closing LOG file...')
    linked_log_file = log_file.close_and_upload()

    # Here we are taking all the files generated by the various functions above and adding them to a single tar
    # to enable easy output. The only output of this applet is thus a single .tar.gz file per VCF file
    LOGGER.info("Generating final tarball...")

    # 1. Build the full tar command with compression (`-z`)
    output_tarball = Path(f"{output_prefix}.tar.gz")
    cmd = ["tar", "-czf", str(output_tarball)]
    # Add each file path to the tar command
    for file_path in output_files:
        cmd.append(str(file_path))

    # 2. Log the exact command for clarity/troubleshooting
    LOGGER.info("Running command: %s", " ".join(cmd))

    # 3. Run the command, raising an error if tar returns non-zero
    subprocess.run(cmd, check=True)

    LOGGER.info("Tarball created at: %s", output_tarball)

    # Set output
    output = {'output_tarball': dxpy.dxlink(generate_linked_dx_file(output_tarball)),
              'log_file': dxpy.dxlink(linked_log_file)}

    return output

dxpy.run()
