########################################################################################################################
# This file contains a series of helper methods that are used to parallelize the collapsing process across all BGEN files.
# These methods typically do not contain much functionality themselves, but rather wrap other methods in a way that
# allows for easier calling of the wrapped methods in parallel.
#
# This means that some of these methods are not unit tested as they are either 1) wrappers around methods that are or
# 2) are methods that directly require DNANexus to run.
########################################################################################################################

from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.mrc_logger import MRCLogger
from scipy.io import mmwrite
from scipy.sparse import csr_matrix, hstack

from utilities.collapse_logger import CollapseLOGGER
from utilities.collapse_utils import check_matrix_stats, \
    stat_writer
from collapsevariants.tool_parsers.bolt_parser import BOLTParser
from collapsevariants.tool_parsers.regenie_parser import REGENIEParser
from collapsevariants.tool_parsers.saige_parser import SAIGEParser
from collapsevariants.tool_parsers.staar_parser import STAARParser

LOGGER = MRCLogger(__name__).get_logger()


def update_log_file(genes: Dict[str, pd.DataFrame], genotype_index: Dict[str, csr_matrix], n_samples: int,
                    expected_total_sites: int, stat_logger: CollapseLOGGER) -> None:
    """Update the CollapseLOGGER with per-sample and per-ENST totals across all BGEN files.

    This method also handles the parallelization of the :func:`check_matrix_stats` method across all BGEN files and the
    subsequent concatenation of the results into a single set of totals for the entire dataset.

    :param genes: A dictionary containing the genes to collapse with keys of the BGEN file prefixes and values of
        a Pandas DataFrame containing per-variant information
    :param genotype_index: A dictionary containing values of csr_matrix and keys of each BGEN file prefix.
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


def generate_generic_masks(genes: Dict[str, pd.DataFrame], genotype_index: Dict[str, csr_matrix],
                           sample_ids: List[str], output_prefix: str) -> List[Path]:
    """Wrapper to help generate output files for each tool.

    :param genes: A dictionary containing the genes to collapse with keys of the BGEN file prefixes and values of
        a Pandas DataFrame containing per-variant information.
    :param genotype_index: A dictionary containing values of csr_matrix and keys of each BGEN file prefix.
    :param sample_ids: A list of sample IDs for processing this file.
    :param output_prefix: A string representing the prefix of the output files.
    :return: A list of Path objects representing the output files created by the implementing classes.
    """

    # Generate output files for each tool
    output_files = []
    tool_methods = [BOLTParser, SAIGEParser, REGENIEParser, STAARParser]
    for tool in tool_methods:
        tool_instance = tool(genes, genotype_index, sample_ids, output_prefix)
        output_files.extend(tool_instance.get_output_files())

    return output_files


def generate_snp_or_gene_masks(genes: Dict[str, pd.DataFrame], genotype_index: Dict[str, csr_matrix],
                               sample_ids: List[str], output_prefix: str, bgen_type: str) -> List[Path]:
    """
    Wrapper similar to generate_generic_masks, but for SNP and GENE masks, creating output inputs for various tools.

    This method takes the output of the collapsing process and 'stacks' the resulting matrices. Since only one matrix
    is required for each GENE / SNP mask, we do not need to create multiple outputs for each .bgen as when running
    a filtering expression.

    :param genes: A dictionary containing the genes to collapse with keys of the BGEN file prefixes and values of
        a Pandas DataFrame containing per-variant information.
    :param genotype_index: A dictionary containing values of csr_matrix and keys of each BGEN file prefix.
    :param sample_ids: A list of sample IDs for processing this file.
    :param output_prefix: A string representing the prefix of the output files.
    :param bgen_type: A string representing the type of BGEN file (e.g., 'SNP' or 'GENE').
    :return: A list of Path objects representing the output files created by the implementing classes.
    """
    # Need to concatenate all matrices into a single matrix while ensuring concatenation is done in the same order for
    # variant indices AND genotype matrices

    # 1. Collect submatrices and variant indices
    final_variant_index_list = []
    matrix_list = []

    for bgen_prefix, variant_index in genes.items():
        current_matrix = genotype_index[bgen_prefix]
        # Collect the variant index DataFrame
        final_variant_index_list.append(variant_index)
        # Collect the sparse matrix for later concatenation
        matrix_list.append(current_matrix[0])

    # 2. Concatenate all variant indices
    final_variant_index = pd.concat(final_variant_index_list)

    # 3. Perform one hstack on the list of sparse matrices
    final_genotype_matrix = hstack(matrix_list)

    # 4. Write the final data to disk (Matrix Market format)
    matrix_output_path = Path(f'{output_prefix}.{bgen_type}.STAAR.mtx')

    ## make sure the output is in the correct format
    mmwrite(matrix_output_path, final_genotype_matrix)

    sample_output_path = STAARParser.make_samples_dict(output_prefix, bgen_type, sample_ids)
    variant_output_path = STAARParser.make_variants_dict(output_prefix, bgen_type, final_variant_index)

    return [matrix_output_path, sample_output_path, variant_output_path]
