from pathlib import Path
from typing import Tuple, Dict, List

import numpy as np
import pandas as pd
from bgen import BgenReader
from general_utilities.mrc_logger import MRCLogger
from scipy.sparse import coo_matrix, csr_matrix

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
                                  temp_dir=None, max_genes_per_batch=10) -> csr_matrix:
    """Generate a sparse matrix of genotypes from a BGEN file."""
    search_list = variant_list.groupby('ENST').aggregate(
        CHROM=('CHROM', 'first'),
        MIN=('POS', 'min'),
        MAX=('POS', 'max'),
        VARS=('varID', list)
    )

    j_lookup = {var: idx for idx, var in enumerate(variant_list['varID'])}
    temp_dir = temp_dir or Path('/tmp')

    with BgenReader(bgen_path, sample_path=sample_path, delay_parsing=True) as bgen_reader:
        num_samples = len(bgen_reader.samples)
        num_variants = len(variant_list)

        # Pre-allocate memory-mapped arrays with estimated size
        estimated_nonzero = len(variant_list) * num_samples // 10  # Assume ~10% density
        i_array = np.memmap(temp_dir / 'i_array.memmap', dtype=np.int64, mode='w+', shape=(estimated_nonzero,))
        j_array = np.memmap(temp_dir / 'j_array.memmap', dtype=np.int64, mode='w+', shape=(estimated_nonzero,))
        d_array = np.memmap(temp_dir / 'd_array.memmap', dtype=np.float64, mode='w+', shape=(estimated_nonzero,))

        pos = 0
        for batch_start in range(0, len(search_list), max_genes_per_batch):
            batch_end = min(batch_start + max_genes_per_batch, len(search_list))
            batch_genes = search_list.iloc[batch_start:batch_end]

            for gene in batch_genes.itertuples():
                chrom = str(gene.CHROM).replace("chr", "").strip()
                chrom = int(chrom) if chrom not in ['X', 'Y'] else chrom
                variants = bgen_reader.fetch(chrom, gene.MIN, gene.MAX)

                for variant in variants:
                    if variant.rsid in gene.VARS:
                        prob = variant.probabilities
                        genotypes = np.where(prob[:, 1] == 1, 1., np.where(prob[:, 2] == 1, 2., 0.))
                        nonzero_idx = genotypes.nonzero()[0]

                        if len(nonzero_idx) > 0:
                            # Resize arrays if needed
                            if pos + len(nonzero_idx) > len(i_array):
                                new_size = int(len(i_array) * 1.5)
                                i_array.flush()
                                j_array.flush()
                                d_array.flush()
                                i_array = np.memmap(temp_dir / 'i_array.memmap', dtype=np.int64, mode='r+',
                                                   shape=(new_size,))
                                j_array = np.memmap(temp_dir / 'j_array.memmap', dtype=np.int64, mode='r+',
                                                   shape=(new_size,))
                                d_array = np.memmap(temp_dir / 'd_array.memmap', dtype=np.float64, mode='r+',
                                                   shape=(new_size,))

                            chunk_size = len(nonzero_idx)
                            i_array[pos:pos + chunk_size] = nonzero_idx
                            j_array[pos:pos + chunk_size] = j_lookup[variant.rsid]
                            d_array[pos:pos + chunk_size] = genotypes[nonzero_idx]
                            pos += chunk_size

        # Create final matrix
        matrix = csr_matrix((d_array[:pos], (i_array[:pos], j_array[:pos])),
                           shape=(num_samples, num_variants))

        # Clean up
        for name in ['i_array.memmap', 'j_array.memmap', 'd_array.memmap']:
            (temp_dir / name).unlink()

        return matrix


def check_matrix_stats(genotypes: csr_matrix, variant_list: pd.DataFrame) -> Tuple[
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

    ac_table = np.zeros(genotypes.shape[0])
    gene_ac_table = np.zeros(genotypes.shape[0])

    ENSTs = variant_list['ENST'].unique()
    gene_totals = dict.fromkeys(variant_list['ENST'].unique(), 0)  # Quickly make a blank dictionary

    # We iterate per-ENST here, as we need to calculate both per-sample and per-gene totals, so no reason to iterate
    # twice.
    for ENST in ENSTs:
        current_variant_list = variant_list[variant_list['ENST'] == ENST]
        current_genotypes = genotypes[:, current_variant_list.index]

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