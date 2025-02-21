from pathlib import Path
from pathlib import Path
from typing import Tuple, Dict, List

import numpy as np
import pandas as pd
from bgen import BgenReader
from general_utilities.mrc_logger import MRCLogger
from scipy.sparse import csr_matrix

from collapsevariants.collapse_logger import CollapseLOGGER

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


def generate_csr_matrix_from_bgen(variant_list: pd.DataFrame, bgen_path: Path, sample_path: Path,
                                  uncollapsed_matrix: bool = False) -> tuple:
    """Convert BGEN genotypes into a sparse matrix format.

    Creates a CSR matrix from BGEN file genotypes for use in STAAR/GLM association testing.
    The matrix represents samples as rows and either genes or variants as columns.
    Non-alternate genotypes are filtered out during conversion.

    Args:
        variant_list: DataFrame containing variants to extract
        bgen_path: Path to BGEN file
        sample_path: Path to sample file
        uncollapsed_matrix: If True, keep variants separate. If False, sum variants per gene.

    Returns:
        tuple: (csr_matrix, summary_dict) where csr_matrix has shape (n_samples, n_genes_or_variants)
    """

    # First aggregate across genes to generate a list of genes and their respective variants
    search_list = variant_list.groupby('ENST').aggregate(
        CHROM=('CHROM', 'first'),
        MIN=('POS', 'min'),
        MAX=('POS', 'max'),
        VARS=('varID', list)
    )

    with BgenReader(bgen_path, sample_path=sample_path, delay_parsing=True) as bgen_reader:
        current_samples = np.array(bgen_reader.samples)

        # create a list to store genotype arrays
        genotype_arrays = []
        summary_dict = {}

        # iterate through each gene in our search list
        for gene_n, current_gene in enumerate(search_list.itertuples()):

            # implement a fix to ensure we are pulling out chromosomes as integers
            # chrom = current_gene.CHROM
            # if isinstance(chrom, str):
            #     chrom = chrom.replace("chr", "").strip()
            #     if chrom not in ['X', 'Y']:
            #         chrom = int(chrom)

            # get the actual data from the bgen file
            variants = bgen_reader.fetch(current_gene.CHROM, current_gene.MIN, current_gene.MAX)

            # create a store for the variant level information
            variant_arrays = []

            # collect genotype arrays for each variant
            for current_variant in variants:
                LOGGER.info('current variant')

                if current_variant.rsid in current_gene.VARS:
                    # pull out the actual genotypes
                    current_probabilities = current_variant.probabilities

                    # store variant codings
                    variant_array = np.where(current_probabilities[:, 1] == 1, 1.,
                                             np.where(current_probabilities[:, 2] == 1, 2., 0.))

                    # store variant level information in the array we created
                    variant_arrays.append(variant_array)

            # stack the variant information for all variants in the gene
            stacked_variants = np.column_stack(variant_arrays)
            LOGGER.info('finished appending variants')

            # if we are collapsing here (to save on memory), leave as False. If set to True, we won't collapse
            # and instead the uncollapsed stacked variants will be appended
            # note the vector naming convention in this small section is a bit hacky but we want the vectors naming
            # to be consistent so it works for the rest of the function
            if uncollapsed_matrix is not True:
                stacked_variants = stacked_variants.sum(axis=1)

            # append the variant arrays to the genotype array
            genotype_arrays.append(stacked_variants)
            LOGGER.info('finished stacking variants')

            # record the per-gene stats in a dict
            summary_dict[current_gene.Index] = {
                'sum': np.sum(stacked_variants),  # use np.sum to get total of all values
                'variants_per_gene': len(variant_arrays),  # get number of variants
                'gene_index': gene_n,
            }

        # stack all genotype arrays into a matrix (samples × variants)
        final_genotypes = np.column_stack(genotype_arrays)
        LOGGER.info('finished stacking genotypes')

        # convert this to a csr matrix
        LOGGER.info('finished making csr matrix')

        final_genotypes = csr_matrix(final_genotypes, shape=(len(current_samples), len(search_list)))

    return final_genotypes, summary_dict


def check_matrix_stats(genotypes: tuple, variant_list: pd.DataFrame) -> Tuple[
    np.ndarray, np.ndarray, Dict[str, int]]:
    """Get information relating to included variants in csr_matrix format.

    This method calculates per-sample and per-gene totals for this chromosome

    Remember: We don't need the actual samples here, we just need to ensure sample order when calculations are performed
    and stored is consistent with the order of the samples in the BGEN file.

    :param genotypes: A csr_matrix containing the genotypes for this bgen file. i = samples, j = variants.
    :param variant_list: A pandas.DataFrame containing the list of variants for the purposes of calculating per-gene
        totals.
    :return: A Tuple containing three items: 1) per-sample allele counts, 2) Number of genes affected per-sample,
        3) A dictionary containing the total number of variants per transcript. The first two are numpy arrays to enable
        fast additions across multiple bgen files
    """
    genotype_matrix, summary_dict = genotypes

    ac_table = np.zeros(genotype_matrix.shape[0])
    gene_ac_table = np.zeros(genotype_matrix.shape[0])

    ENSTs = variant_list['ENST'].unique()
    gene_totals = dict.fromkeys(variant_list['ENST'].unique(), 0)  # Quickly make a blank dictionary

    # We iterate per-ENST here, as we need to calculate both per-sample and per-gene totals, so no reason to iterate
    # twice.
    for ENST in ENSTs:
        current_genotypes = genotype_matrix[:, summary_dict[ENST]['gene_index']]

        sample_sums = current_genotypes.sum(axis=1).A1  # Get genotype totals per-sample
        gene_sums = np.where(sample_sums > 0., 1., 0.)  # Get total number of individuals with at least 1 allele
        gene_totals[ENST] = current_genotypes.shape[1]  # Get total number of variants per gene

        ac_table = np.add(ac_table, sample_sums)
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
