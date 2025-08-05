########################################################################################################################
# This file contains a series of helper methods that are used to parallelize the collapsing process across all BGEN files.
# These methods typically do not contain much functionality themselves, but rather wrap other methods in a way that
# allows for easier calling of the wrapped methods in parallel.
#
# This means that these methods are not unit tested.
########################################################################################################################

from typing import Dict, Tuple

import numpy as np
import pandas as pd
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.mrc_logger import MRCLogger
from scipy.sparse import csr_matrix, hstack

from collapsevariants.utilities.collapse_logger import CollapseLOGGER
from collapsevariants.utilities.collapse_utils import check_matrix_stats, \
    stat_writer
from general_utilities.bgen_utilities.genotype_matrix import generate_csr_matrix_from_bgen, make_variant_list
from collapsevariants.utilities.ingest_data import BGENIndex, download_bgen
from utilities.collapse_utils import GenotypeInfo

LOGGER = MRCLogger(__name__).get_logger()


def generate_genotype_matrices(genes: Dict[str, pd.DataFrame], bgen_index: Dict[str, BGENIndex]) -> Dict[str, Tuple[csr_matrix, Dict[str, GenotypeInfo]]]:
    """Helper method for parellelizing :func:`generate_genotype_matrix` across all BGEN files with at least one variant.

    This method generates csr_matrices for each BGEN file in the input dictionary of genes. It simply wraps the
    :func:`generate_genotype_matrix` method in a ThreadUtility object to parallelize the process.

    Note that this method is un-tested as it wraps a method that requires DNA Nexus to run.

    :return: A pandas.DataFrame containing per-sample and per-ENST totals for log reporting purposes.
    """

    # Generate genotype matrices for each BGEN file
    thread_utility = ThreadUtility(error_message='Error in generation genotype arrays', incrementor=10)
    for bgen_prefix in genes.keys():
        thread_utility.launch_job(generate_genotype_matrix,
                                  bgen_prefix=bgen_prefix,
                                  chrom_bgen_index=bgen_index[bgen_prefix],
                                  variant_list=genes[bgen_prefix])
    genotype_index = {bgen_prefix: (geno_matrix, summary_dict) for bgen_prefix, geno_matrix, summary_dict in
                      thread_utility}
    return genotype_index


def generate_genotype_matrix(bgen_prefix: str, chrom_bgen_index: BGENIndex,
                             variant_list: pd.DataFrame, delete_on_complete: bool = True) -> Tuple[str, csr_matrix, Dict[str, GenotypeInfo]]:
    """
    Helper method that wraps :func:`generate_csr_matrix_from_bgen` to generate a genotype matrix for a single BGEN file.

    This wrapper method is used to parallelize the generation of genotype matrices across all BGEN files with at least one
    variant. We don't parallelize :func:`generate_csr_matrix_from_bgen` directly to allow for :func:`download_bgen` to
    be separated out and allow for unit testing of :func:`generate_csr_matrix_from_bgen` detached from DNANexus.

    :param bgen_prefix: A string representing the prefix of the BGEN file to run in this current thread.
    :param chrom_bgen_index: A BGENIndex object containing the paths to the BGEN file, BGEN index file, and BGEN sample file.
    :param variant_list: A pandas.DataFrame containing the variants to collapse on.
    :param delete_on_complete: If True, delete the BGEN, index, and sample files after processing. Required for testing purposes.
        Default is True.
    :return: A tuple containing the BGEN file prefix (for thread tracking) and the csr_matrix generated from the
        BGEN file.
    """
    # Note that index is required but is not explicitly taken as input by BgenReader. It MUST have the same
    # name as the bgen file, but with a .bgi suffix.
    bgen_path, index_path, sample_path = download_bgen(chrom_bgen_index)

    variant_list = make_variant_list(variant_list)

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
                                                                          should_collapse_matrix=True)

        # Build the genotype matrix
        genotypes.append(gene_genotypes)

        # Build the summary dict
        current_end = current_start + gene_summary_dict['n_columns']
        summary_dict[gene] = GenotypeInfo(
            allele_count=gene_summary_dict['allele_count'],
            n_variants=gene_summary_dict['n_variants'],
            n_columns=gene_summary_dict['n_columns'],
            gene_index=[var_n for var_n in range(current_start, current_end)]
        )
        current_start = current_end

    if delete_on_complete:

        bgen_path.unlink()
        index_path.unlink()
        sample_path.unlink()

    # Finalise matrix creation
    genotypes = hstack(genotypes)

    return bgen_prefix, genotypes, summary_dict


def update_log_file(genes: Dict[str, pd.DataFrame], genotype_index: Dict[str, Tuple[csr_matrix, Dict[str, GenotypeInfo]]],
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
    thread_utility = ThreadUtility(error_message='Error in generation of per-bgen totals', incrementor=10)
    for bgen_prefix in genes.keys():
        thread_utility.launch_job(check_matrix_stats,
                                  genotypes=genotype_index[bgen_prefix],
                                  variant_list=genes[bgen_prefix])

    ac_table = np.zeros(n_samples)
    gene_ac_table = np.zeros(n_samples)
    gene_totals = dict()
    for result in thread_utility:
        ac_table = np.add(ac_table, result[0])
        gene_ac_table = np.add(gene_ac_table, result[1])
        gene_totals.update(result[2])

    stat_writer(ac_table, gene_ac_table, gene_totals, expected_total_sites, stat_logger)


