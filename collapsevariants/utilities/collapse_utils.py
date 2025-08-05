from pathlib import Path
from typing import Tuple, Dict, List, TypedDict

import numpy as np
import pandas as pd
from general_utilities.mrc_logger import MRCLogger
from scipy.sparse import csr_matrix

from collapsevariants.utilities.collapse_logger import CollapseLOGGER

LOGGER = MRCLogger().get_logger()


def get_sample_ids(sample_path: Path) -> List[str]:
    """
    Extract sample IDs from a sample file.

    This function reads a sample file and extracts the sample IDs, skipping the first two header lines.

    :param sample_path: Path to the sample file.
    :return: A list of sample IDs.
    """
    with sample_path.open('r') as sample_reader:
        sample_ids = [line.split()[0] for line in sample_reader.readlines()]
        sample_ids = sample_ids[2:]

    return sample_ids


class GenotypeInfo(TypedDict):
    """
    TypedDict to hold information about the genotype matrix and summary statistics.

    :attr allele_count: The total number of alleles in the matrix.
    :attr n_variants: The total number of variants in the matrix.
    :attr n_columns: The total number of columns representing this gene (may not == n_variants if collapsing the matrix).
    :attr gene_index: A list of indices corresponding to the variants in the matrix.
    """

    allele_count: int
    n_variants: int
    n_columns: int
    gene_index: List[int]


def check_matrix_stats(genotypes: Tuple[csr_matrix, Dict[str, GenotypeInfo]], variant_list: pd.DataFrame) -> Tuple[
    np.ndarray, np.ndarray, Dict[str, int]]:
    """
    Get information relating to included variants in csr_matrix format.

    This method calculates per-sample and per-gene totals for this chromosome.

    Remember: We don't need the actual samples here, we just need to ensure sample order when calculations are performed
    and stored is consistent with the order of the samples in the BGEN file.

    :param genotypes: A tuple containing:
        - A csr_matrix with genotypes for this BGEN file (samples as rows, variants as columns).
        - A summary dictionary with information for each gene.
    :param variant_list: A pandas DataFrame containing the list of variants for calculating per-gene totals.
    :return: A tuple containing three items:
        1) A numpy array with per-sample allele counts.
        2) A numpy array with the number of genes affected per sample.
        3) A dictionary with the total number of variants per transcript.
    """

    genotype_matrix, summary_dict = genotypes

    ac_table = np.zeros(genotype_matrix.shape[0])
    gene_ac_table = np.zeros(genotype_matrix.shape[0])

    ENSTs = variant_list['ENST'].unique()
    gene_totals = dict.fromkeys(variant_list['ENST'].unique(), 0)  # Quickly make a blank dictionary

    # We iterate per-ENST here, as we need to calculate both per-sample and per-gene totals, so no reason to iterate
    # twice.
    for ENST in ENSTs:
        # Extract genotypes for the current gene
        current_genotypes = genotype_matrix[:, summary_dict[ENST]['gene_index']]

        # Calculate per-sample genotype totals
        sample_sums = current_genotypes.sum(axis=1).A1

        # Determine the number of samples with at least one allele
        gene_sums = np.where(sample_sums > 0., 1., 0.)

        # Record the total number of variants for the current gene
        gene_totals[ENST] = summary_dict[ENST]['n_variants']

        # Update the allele count table with per-sample totals
        ac_table = np.add(ac_table, sample_sums)

        # Update the gene count table with the number of genes affected per sample
        gene_ac_table = np.add(gene_ac_table, gene_sums)

    return ac_table, gene_ac_table, gene_totals


def stat_writer(ac_table: np.ndarray, gene_ac_table: np.ndarray, gene_totals: Dict[str, int],
                expected_total_sites: int, log_file: CollapseLOGGER) -> None:
    """ Writes stats about the various collapsing operations performed by this applet.

    :param ac_table: A numpy array containing the number of alleles per individual
    :param gene_ac_table: A numpy array containing the number of genes affected per individual
    :param gene_totals: A dictionary containing the total number of variants per gene
    :param expected_total_sites: Total number of expected sites based on the original query
    :param log_file: The LOG_FILE for this instance to print to
    :return: None
    """

    found_total_sites = sum(gene_totals.values())

    log_file.write_header('Genome-wide totals')
    log_file.write_int('Total sites expected from filtering expression', expected_total_sites)
    log_file.write_int('Total sites extracted from all chromosomes', found_total_sites)
    log_file.write_string('Total expected and total extracted match', str(expected_total_sites == found_total_sites))
    log_file.write_spacer()

    # Concatenate and sum sample tables:
    log_file.write_header('Per-individual totals')
    log_file.write_int('Median number of alleles per indv', np.median(ac_table))
    log_file.write_int('Median number of genes affected per indv', np.median(gene_ac_table))
    log_file.write_float('Mean number of alleles per indv', np.mean(ac_table))
    log_file.write_float('Mean number of genes affected per indv', np.mean(gene_ac_table))
    log_file.write_int('Max number of alleles', np.max(ac_table))
    log_file.write_int('Number of individuals with at least 1 allele', len(ac_table.nonzero()[0]))
    log_file.write_spacer()

    log_file.write_header('AC Histogram')
    log_file.write_generic('AC_bin\tcount')
    acs, counts = np.unique(ac_table, return_counts=True)
    for i in range(len(acs)):
        log_file.write_histogram(acs[i], counts[i])
    log_file.write_spacer()

    log_file.write_header('per-ENST Totals')
    log_file.write_generic('ENST\tTotal_variants')
    for gene, total in gene_totals.items():
        log_file.write_histogram(gene, total)
