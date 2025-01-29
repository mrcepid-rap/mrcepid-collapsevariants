from pathlib import Path
from typing import Tuple, Dict, List

import numpy as np
import pandas as pd
from bgen import BgenReader
from general_utilities.association_resources import download_dxfile_by_name
from general_utilities.mrc_logger import MRCLogger
from scipy.sparse import coo_matrix, csr_matrix

from collapsevariants.collapse_logger import CollapseLOGGER
from collapsevariants.ingest_data import BGENIndex

LOGGER = MRCLogger().get_logger()


def download_bgen(chrom_bgen_index: BGENIndex) -> Tuple[Path, Path, Path]:
    """Download the BGEN file from DNANexus

    This method downloads the BGEN file from DNANexus using the bgen_index dictionary and the bgen_prefix. It then
    returns the path to the downloaded BGEN file, index, and sample.

    Note that this method is the only method untested in this package as it requires a DNANexus connection to test.

    :return: A Tuple containing the paths to the BGEN file, index file, and sample file
    """

    # Download the requisite files for this chromosome according to the index dict:
    bgen_path = download_dxfile_by_name(chrom_bgen_index['bgen'], print_status=False)
    index_path = download_dxfile_by_name(chrom_bgen_index['index'], print_status=False)
    sample_path = download_dxfile_by_name(chrom_bgen_index['sample'], print_status=False)

    return bgen_path, index_path, sample_path


def get_sample_ids(sample_path: Path) -> List[str]:
    with sample_path.open('r') as sample_reader:
        sample_ids = [line.split()[0] for line in sample_reader.readlines()]
        sample_ids = sample_ids[2:]

    return sample_ids


def generate_csr_matrix_from_bgen(variant_list: pd.DataFrame, bgen_path: Path, sample_path: Path) -> csr_matrix:
    """Generate a sparse matrix of genotypes from a BGEN file.

    Independent of a specific method, the goal of this method is to convert bgen outputs into a sparse matrix format
    that is manipulable outside bgen format. In reality, this method enables the use of STAAR / GLMs in downstream
    association testing software.

    We first extract all variants and filter out individuals with a non-alternate genotype. This is converted into a
    coo_matrix for the purposes of efficient creation, and then immediately converted into a csr_matrix for efficient
    slicing operations on the columns (i.e., variants) in the matrix. We choose a csr_matrix as sample-level filtering
    is NOT done by methods that use this matrix.

    :param variant_list: A Pandas DataFrame containing the list of variants to extract from the BGEN file.
    :param bgen_path: A path to a local bgen file for conversion.
    :param sample_path: The associated sample file for the bgen file.
    :return: A csr_matrix with columns (j) representing variants and rows (i) representing samples.
    """

    # First aggregate across genes to generate a list of genes and their respective variants
    search_list = variant_list.groupby('ENST').aggregate(
        CHROM=('CHROM', 'first'),
        MIN=('POS', 'min'),
        MAX=('POS', 'max'),
        VARS=('varID', list)
    )

    j_lookup = variant_list[['varID']]
    j_lookup = j_lookup.reset_index()
    j_lookup = j_lookup.set_index('varID').to_dict(orient='index')

    with BgenReader(bgen_path, sample_path=sample_path, delay_parsing=True) as bgen_reader:
        current_samples = np.array(bgen_reader.samples)

        i_array = []
        j_array = []
        d_array = []

        for current_gene in search_list.itertuples():

            # Note to future devs: it is MUCH faster to fetch a window of variants and iterate through them, checking if
            # they are in our list of variants, then to iterate through the list of variants and fetch each one individually.
            # Important to note that this holds true for SNP and GENE masks as well, as we store the original data in
            # variant_list by gene and LATER modify the name of the ENST to be our dummy values.
            variants = bgen_reader.fetch(current_gene.CHROM, current_gene.MIN, current_gene.MAX)

            for current_variant in variants:
                if current_variant.rsid in current_gene.VARS:
                    current_probabilities = current_variant.probabilities

                    genotype_array = np.where(current_probabilities[:, 1] == 1, 1.,
                                              np.where(current_probabilities[:, 2] == 1, 2., 0.))
                    current_i = genotype_array.nonzero()[0]
                    current_j = [j_lookup[current_variant.rsid]['index']] * len(current_i)
                    current_d = genotype_array[current_i].tolist()

                    i_array.extend(current_i.tolist())
                    j_array.extend(current_j)
                    d_array.extend(current_d)

    genotypes = coo_matrix((d_array, (i_array, j_array)), shape=(len(current_samples), len(variant_list)))
    genotypes = csr_matrix(genotypes)  # Convert to csr_matrix for quick slicing / calculations of variants

    return genotypes


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
