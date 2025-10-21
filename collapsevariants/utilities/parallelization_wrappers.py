########################################################################################################################
# This file contains a series of helper methods that are used to parallelize the collapsing process across all BGEN files.
# These methods typically do not contain much functionality themselves, but rather wrap other methods in a way that
# allows for easier calling of the wrapped methods in parallel.
#
# This means that these methods are not unit tested.
########################################################################################################################

from typing import Dict, Tuple

import dxpy
import numpy as np
import pandas as pd
from general_utilities.bgen_utilities.genotype_matrix import generate_csr_matrix_from_bgen, make_variant_list
from general_utilities.import_utils.file_handlers.export_file_handler import ExportFileHandler
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
from general_utilities.import_utils.import_lib import BGENInformation
from general_utilities.job_management.joblauncher_factory import joblauncher_factory
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.mrc_logger import MRCLogger
from scipy.sparse import csr_matrix, hstack, save_npz, load_npz

from collapsevariants.utilities.collapse_logger import CollapseLOGGER
from collapsevariants.utilities.collapse_utils import GenotypeInfo
from collapsevariants.utilities.collapse_utils import check_matrix_stats, \
    stat_writer
from collapsevariants.utilities.ingest_data import download_bgen

LOGGER = MRCLogger(__name__).get_logger()


def generate_genotype_matrices(genes: Dict[str, pd.DataFrame], bgen_index: Dict[str, BGENInformation],
                               should_collapse=True) -> Dict[str, Tuple[csr_matrix, Dict[str, GenotypeInfo]]]:
    """Helper method for parellelizing :func:`generate_genotype_matrix` across all BGEN files with at least one variant.

    This method generates csr_matrices for each BGEN file in the input dictionary of genes. It simply wraps the
    :func:`generate_genotype_matrix` method in a ThreadUtility object to parallelize the process.

    Note that this method is un-tested as it wraps a method that requires DNA Nexus to run.

    :param genes: A dictionary containing the genes to collapse with keys of the BGEN file prefixes and values of
        a Pandas DataFrame containing per-variant information.
    :param bgen_index: A dictionary containing BGENInformation objects for each BGEN file prefix.
    :param should_collapse: If True, collapse the matrix to remove redundant columns. Default is True.
    :return: A pandas.DataFrame containing per-sample and per-ENST totals for log reporting purposes.
    """

    # Generate genotype matrices for each BGEN file in parallel

    # set the launcher
    launcher = joblauncher_factory()

    # set the exporter
    exporter = ExportFileHandler(delete_on_upload=False)

    for bgen_prefix in genes.keys():

        LOGGER.info(f'Getting files ready for {bgen_prefix}')

        # variant list is a df that we need to export and upload
        genes[bgen_prefix].to_csv(f"{bgen_prefix}.tsv", sep='\t', index=False)
        variant_list = exporter.export_files(f"{bgen_prefix}.tsv")

        launcher.launch_job(function=generate_genotype_matrix,
                            inputs={
                                'bgen_prefix': bgen_prefix,
                                'bgen': bgen_index[bgen_prefix]['bgen'].get_input_str(),
                                'index': bgen_index[bgen_prefix]['index'].get_input_str(),
                                'sample': bgen_index[bgen_prefix]['sample'].get_input_str(),
                                'variant_list': variant_list,
                                'should_collapse': should_collapse
                            },
                            outputs=
                            ['bgen_prefix', 'genotypes', 'summary_dict']
                            )
    launcher.submit_and_monitor()

    genotype_index = {}

    for result in launcher:
        bgen_prefix = result['bgen_prefix']
        geno_matrix = load_npz(result['genotypes'])
        summary_dict = result['summary_dict']

        genotype_index[bgen_prefix] = (geno_matrix, summary_dict)

    return genotype_index


@dxpy.entry_point('generate_genotype_matrix')
def generate_genotype_matrix(bgen_prefix: str, bgen: str, index: str, sample: str,
                             variant_list: dict, should_collapse=True, delete_on_complete: bool = True) -> dict:
    """
    Helper method that wraps :func:`generate_csr_matrix_from_bgen` to generate a genotype matrix for a single BGEN file.

    This wrapper method is used to parallelize the generation of genotype matrices across all BGEN files with at least one
    variant. We don't parallelize :func:`generate_csr_matrix_from_bgen` directly to allow for :func:`download_bgen` to
    be separated out and allow for unit testing of :func:`generate_csr_matrix_from_bgen` detached from DNANexus.

    :param bgen_prefix: A string representing the prefix of the BGEN file to run in this current thread.
    :param bgen: dxlink to download the bgen file.
    :param index: dxlink to download the bgen index file.
    :param sample: dxlink to download the sample file.
    :param variant_list: A dxlink to a pandas.DataFrame containing the variants to collapse on.
    :param should_collapse: If True, collapse the matrix to remove redundant columns. Default is True.
    :param delete_on_complete: If True, delete the BGEN, index, and sample files after processing. Required for testing purposes.
        Default is True.
    :return: A tuple containing the BGEN file prefix (for thread tracking) and the csr_matrix generated from the
        BGEN file.
    """

    # download the files needed to run the subjob
    bgen_path = InputFileHandler(bgen, download_now=True).get_file_handle()
    index_path = InputFileHandler(index, download_now=True).get_file_handle()
    sample_path = InputFileHandler(sample, download_now=True).get_file_handle()
    variants_file = InputFileHandler(variant_list, download_now=True).get_file_handle()
    variant_list = pd.read_csv(variants_file, sep='\t')
    print(variant_list.head())

    print('here')
    variant_list = make_variant_list(variant_list)

    print('here2')
    # Generate the CSR matrix from the BGEN file
    summary_dict = {}
    genotypes = []
    current_start = 0

    for gene, gene_information in variant_list.items():
        gene_genotypes, gene_summary_dict = generate_csr_matrix_from_bgen(bgen_path, sample_path,
                                                                          variant_filter_list=gene_information['vars'],
                                                                          chromosome=gene_information['chrom'],
                                                                          start=gene_information['min'],
                                                                          end=gene_information['max'],
                                                                          should_collapse_matrix=should_collapse)

        print('here3')

        # Build the genotype matrix
        genotypes.append(gene_genotypes)

        print('here4')

        # Build the summary dict
        current_end = current_start + gene_summary_dict['n_columns']
        summary_dict[gene] = GenotypeInfo(
            allele_count=gene_summary_dict['allele_count'],
            n_variants=gene_summary_dict['n_variants'],
            n_columns=gene_summary_dict['n_columns'],
            gene_index=[var_n for var_n in range(current_start, current_end)]
        )
        current_start = current_end

        print('here5')

    if delete_on_complete:
        bgen_path.unlink()
        index_path.unlink()
        sample_path.unlink()

    # Finalise matrix creation
    genotypes = hstack(genotypes)

    # save to file
    output_path = f"{bgen_prefix}_genotypes.npz"
    save_npz(output_path, genotypes)

    print('here6')

    print(bgen_prefix)
    print(output_path)
    print(summary_dict)

    return {
        'bgen_prefix': bgen_prefix,
        'genotypes': output_path,
        'summary_dict': summary_dict
    }


def update_log_file(genes: Dict[str, pd.DataFrame],
                    genotype_index: Dict[str, Tuple[csr_matrix, Dict[str, GenotypeInfo]]],
                    n_samples: int, expected_total_sites: int, stat_logger: CollapseLOGGER) -> None:
    """Update the CollapseLOGGER with per-sample and per-ENST totals across all BGEN files.

    This method also handles the parallelization of the :func:`check_matrix_stats` method across all BGEN files and the
    subsequent concatenation of the results into a single set of totals for the entire dataset.

    :param genes: A dictionary containing the genes to collapse with keys of the BGEN file prefixes and values of
        a Pandas DataFrame containing per-variant information
    :param genotype_index: A dictionary containing values of a Tuple(csr_matrix of genotypes, Dict[gene_id] = GenotypeInfo TypedDict)
        and keys of each BGEN file prefix.
    :param n_samples: The number of samples in the BGEN files.
    :param expected_total_sites: The expected total number of sites in the BGEN files provided by
        :func:`SNPListGenerator`.
    :param stat_logger: A CollapseLOGGER object for logging statistics about the collapsing process.
    :return: None
    """

    # Check stats for each genotype matrix
    thread_utility = ThreadUtility(incrementor=10)
    for bgen_prefix in genes.keys():
        thread_utility.launch_job(function=check_matrix_stats,
                                  inputs={
                                      'genotypes': genotype_index[bgen_prefix],
                                      'variant_list': genes[bgen_prefix],
                                  },
                                  outputs=[
                                      'ac_table', 'gene_ac_table', 'gene_totals'
                                  ]
                                  )
        thread_utility.submit_and_monitor()

        ac_table = np.zeros(n_samples)
        gene_ac_table = np.zeros(n_samples)
        gene_totals = dict()
        for result in thread_utility:
            ac_table = np.add(ac_table, result["ac_table"])
            gene_ac_table = np.add(gene_ac_table, result["gene_ac_table"])
            gene_totals.update(result["gene_totals"])

        stat_writer(ac_table, gene_ac_table, gene_totals, expected_total_sites, stat_logger)
