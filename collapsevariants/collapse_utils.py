import os
import tempfile
from pathlib import Path
from typing import Tuple, Dict, List
from general_utilities.job_management.thread_utility import ThreadUtility

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
                                  chunk_size: int = 10) -> csr_matrix:
    """
    Generate a sparse matrix of genotypes from a BGEN file.

    This function converts BGEN file outputs into a sparse matrix format (CSR) that can be used in downstream
    association testing software such as STAAR or GLMs. The process involves extracting all variants, filtering out
    individuals with a non-alternate genotype, and converting the data into a COO matrix for efficient creation,
    which is then converted into a CSR matrix for efficient slicing operations on the columns (variants).

    Parameters:
    variant_list (pd.DataFrame): A Pandas DataFrame containing the list of variants to extract from the BGEN file.
                                 It should include columns such as 'ENST', 'CHROM', 'POS', and 'varID'.
    bgen_path (Path): The path to the local BGEN file for conversion.
    sample_path (Path): The path to the associated sample file for the BGEN file.
    chunk_size (int): The number of genes to process in each chunk. Default is 100.

    Returns:
    csr_matrix: A sparse matrix with columns (j) representing variants and rows (i) representing samples.

    The function performs the following steps:
    1. Aggregates the variant list by gene to get minimal region fetch.
    2. Builds a lookup dictionary (j_lookup) mapping variant IDs to column indices.
    3. Creates a temporary directory to store partial results.
    4. Opens the BGEN file and reads the sample list.
    5. Splits the gene list into chunks and processes each chunk to produce partial (i, j, d) arrays, which are saved to disk.
    6. Reads the partial files back and concatenates them in a second pass.
    7. Builds the final arrays for the COO matrix.
    8. Constructs and returns the final CSR matrix.

    Edge Cases:
    - If no data is found at all, the function returns an empty CSR matrix with the appropriate shape.
    """

    # Aggregate by gene to get minimal region fetch
    search_list = variant_list.groupby('ENST').aggregate(
        CHROM=('CHROM', 'first'),
        MIN=('POS', 'min'),
        MAX=('POS', 'max'),
        VARS=('varID', list)
    ).reset_index(drop=True)  # optionally reset_index

    # Build j_lookup: varID -> {'index': col_index}
    j_lookup = (
        variant_list[['varID']]
        .reset_index()
        .set_index('varID')
        .to_dict(orient='index')
    )
    n_variants = len(variant_list)

    # We'll store partial results in a temporary directory
    temp_dir = tempfile.mkdtemp(prefix="bgen_sparse_")

    # Open BGEN once, read the sample list
    with BgenReader(bgen_path, sample_path=sample_path, delay_parsing=True) as bgen_reader:
        samples = np.array(bgen_reader.samples)
        n_samples = len(samples)

        # split the gene list into chunks
        chunk_file_paths = []  # store the paths to partial npz files
        for start_idx in range(0, len(search_list), chunk_size):
            end_idx = start_idx + chunk_size
            chunk_df = search_list.iloc[start_idx:end_idx]

            # We'll name a temporary file for each chunk
            tmp_file_path = os.path.join(temp_dir, f"chunk_{start_idx}.npz")
            chunk_file_paths.append(tmp_file_path)

            # Process this chunk: parse, produce (i,j,d), save to disk
            process_gene_chunk(
                bgen_reader=bgen_reader,
                chunk_df=chunk_df,
                j_lookup=j_lookup,
                tmp_file_path=tmp_file_path
            )

    # Read partial files back and concatenate in a second pass
    print("Reading the files back in")
    # Pre-calculate total number of non-zero elements
    total_nnz = 0
    for file_path in chunk_file_paths:
        with np.load(file_path) as data:
            total_nnz += len(data['i'])

    # Pre-allocate final arrays
    i_array = np.empty(total_nnz, dtype=np.int64)
    j_array = np.empty(total_nnz, dtype=np.int64)
    d_array = np.empty(total_nnz, dtype=np.float32)

    # Fill arrays sequentially
    pos = 0
    for file_path in chunk_file_paths:
        with np.load(file_path) as data:
            chunk_size = len(data['i'])
            if chunk_size > 0:
                i_array[pos:pos + chunk_size] = data['i']
                j_array[pos:pos + chunk_size] = data['j']
                d_array[pos:pos + chunk_size] = data['d']
                pos += chunk_size
        os.remove(file_path)

    os.rmdir(temp_dir)

    # Construct the final CSR matrix
    coo = coo_matrix((d_array, (i_array, j_array)), shape=(n_samples, n_variants))
    return csr_matrix(coo, dtype=np.float32)


def process_gene_chunk(
        bgen_reader,
        chunk_df: pd.DataFrame,
        j_lookup: dict,
        tmp_file_path: str
):
    """
    Process a chunk of genes and write partial (i, j, d) arrays to disk in NPZ format.

    This function processes a subset of genes (rows of chunk_df) by fetching variants from a BGEN file,
    extracting genotype probabilities, and converting them into sparse matrix format. The results are
    saved as partial (i, j, d) arrays to a specified temporary file.

    Parameters:
    bgen_reader (BgenReader): An instance of BgenReader used to fetch variant data from the BGEN file.
    chunk_df (pd.DataFrame): A DataFrame containing a subset of genes to process. Each row should include
                             columns 'CHROM', 'MIN', 'MAX', and 'VARS'.
    j_lookup (dict): A dictionary mapping variant IDs (varID) to their corresponding column indices in the
                     final sparse matrix. The dictionary should have the structure {varID: {'index': col_index}}.
    tmp_file_path (str): The path to the temporary file where the partial results will be saved in NPZ format.

    Returns:
    None: This function does not return any value. It writes the partial results to the specified temporary file.

    The function performs the following steps:
    1. Initializes empty lists to store the row indices (i), column indices (j), and data values (d) for the sparse matrix.
    2. Iterates over each gene in the chunk DataFrame.
    3. For each gene, fetches variants from the BGEN file within the specified chromosomal range.
    4. For each variant, checks if it is in the list of variants to process (gene.VARS).
    5. Extracts genotype probabilities and converts them into a genotype array.
    6. Identifies non-zero entries in the genotype array and maps them to the appropriate row and column indices.
    7. Appends the non-zero entries to the lists of row indices, column indices, and data values.
    8. Concatenates the partial results and saves them to the specified temporary file in NPZ format.
    9. If no variants are found in the chunk, saves empty arrays to the temporary file.

    Edge Cases:
    - If no variants are found in the chunk, the function saves empty arrays to the temporary file.
    - If the chromosome information is in string format, it removes the "chr" prefix and converts it to an integer.
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
