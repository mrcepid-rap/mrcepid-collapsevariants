import os
import tempfile
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
                                  chunk_size: int = 100) -> csr_matrix:
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

    # 1) Aggregate by gene to get minimal region fetch
    search_list = variant_list.groupby('ENST').aggregate(
        CHROM=('CHROM', 'first'),
        MIN=('POS', 'min'),
        MAX=('POS', 'max'),
        VARS=('varID', list)
    ).reset_index(drop=True)  # optionally reset_index

    # 2) Build j_lookup: varID -> {'index': col_index}
    j_lookup = (
        variant_list[['varID']]
        .reset_index()
        .set_index('varID')
        .to_dict(orient='index')
    )
    n_variants = len(variant_list)

    # 3) We'll store partial results in a temporary directory
    temp_dir = tempfile.mkdtemp(prefix="bgen_sparse_")

    # 4) Open BGEN once, read the sample list
    with BgenReader(bgen_path, sample_path=sample_path, delay_parsing=True) as bgen_reader:
        samples = np.array(bgen_reader.samples)
        n_samples = len(samples)

        # 5) Split the gene list into chunks
        chunk_file_paths = []  # store the paths to partial npz files
        for start_idx in range(0, len(search_list), chunk_size):
            end_idx = start_idx + chunk_size
            chunk_df = search_list.iloc[start_idx:end_idx]

            # We'll name a temporary file for each chunk
            tmp_file_path = os.path.join(temp_dir, f"chunk_{start_idx}.npz")
            chunk_file_paths.append(tmp_file_path)

            # Process this chunk: parse, produce (i,j,d), save to disk
            _process_gene_chunk(
                bgen_reader=bgen_reader,
                chunk_df=chunk_df,
                j_lookup=j_lookup,
                tmp_file_path=tmp_file_path
            )

    # 6) Read partial files back and concatenate in a second pass
    all_i = []
    all_j = []
    all_d = []
    for file_path in chunk_file_paths:
        with np.load(file_path) as data:
            i_array = data['i']
            j_array = data['j']
            d_array = data['d']
            if len(i_array) > 0:
                all_i.append(i_array)
                all_j.append(j_array)
                all_d.append(d_array)
        # Optionally delete the chunk file as soon as we’ve used it
        os.remove(file_path)

    # Optionally remove the entire temp directory
    os.rmdir(temp_dir)

    # Edge case: no data found at all
    if len(all_i) == 0:
        return csr_matrix((n_samples, n_variants), dtype=np.float32)

    # 7) Build final arrays
    i_array = np.concatenate(all_i)
    j_array = np.concatenate(all_j)
    d_array = np.concatenate(all_d)

    # 8) Construct the final CSR matrix
    coo = coo_matrix((d_array, (i_array, j_array)), shape=(n_samples, n_variants))
    return csr_matrix(coo, dtype=np.float32)


def _process_gene_chunk(
        bgen_reader,
        chunk_df: pd.DataFrame,
        j_lookup: dict,
        tmp_file_path: str
):
    """
    Process a chunk of genes (rows of chunk_df), and write partial (i, j, d) arrays to disk (NPZ).
    """
    chunks_i = []
    chunks_j = []
    chunks_d = []

    # For each gene in the chunk, fetch variants from the BGEN
    for gene in chunk_df.itertuples():
        chrom = gene.CHROM
        if isinstance(chrom, str):
            chrom = chrom.replace("chr", "").strip()
        chrom = int(chrom)

        variants = bgen_reader.fetch(chrom, gene.MIN, gene.MAX)
        for current_variant in variants:
            if current_variant.rsid in gene.VARS:
                probs = current_variant.probabilities
                genotype_array = np.where(probs[:, 1] == 1, 1.,
                                          np.where(probs[:, 2] == 1, 2., 0.))

                current_i = genotype_array.nonzero()[0]
                if len(current_i) == 0:
                    continue

                # global column index from j_lookup
                col_idx = j_lookup[current_variant.rsid]['index']
                current_j = np.full_like(current_i, fill_value=col_idx, dtype=np.int64)
                current_d = genotype_array[current_i]

                chunks_i.append(current_i)
                chunks_j.append(current_j)
                chunks_d.append(current_d)

    # Concatenate partial results for this chunk
    if len(chunks_i) == 0:
        # Nothing to write (maybe no variants in this chunk)
        np.savez_compressed(
            tmp_file_path,
            i=np.array([], dtype=np.int64),
            j=np.array([], dtype=np.int64),
            d=np.array([], dtype=np.float32)
        )
    else:
        i_array = np.concatenate(chunks_i)
        j_array = np.concatenate(chunks_j)
        d_array = np.concatenate(chunks_d).astype(np.float32, copy=False)
        np.savez_compressed(tmp_file_path, i=i_array, j=j_array, d=d_array)


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
