from pathlib import Path
from typing import Tuple, Dict, Any

import numpy as np
import pandas as pd
from bgen import BgenReader
from general_utilities.mrc_logger import MRCLogger
from scipy.sparse import csr_matrix


LOGGER = MRCLogger(__name__).get_logger()


def generate_csr_matrix_from_bgen(variant_list: pd.DataFrame, bgen_path: Path, sample_path: Path,
                                  should_collapse_matrix: bool = True) -> Tuple[csr_matrix, Dict[str, Any]]:
    """
    Convert BGEN genotypes into a sparse matrix format.

    Creates a CSR matrix from BGEN file genotypes for use in STAAR/GLM association testing.
    The matrix represents samples as rows and either genes or variants as columns.
    Non-alternate genotypes are filtered out during conversion.

    Note that in order to optimise this process we are collapsing the variants here (stored in the CSR matrix part
    of our output tuple). We are also saving gene and sample level variant information in a dictionary,
    so that this information can be used downstream (for example, when generating the log file).
    There is also a functionality where if should_collapse_matrix is set to False, then we append the entire variant
    matrix in a CSR format. This is used in the associationtesting modules downstream. The rationale for doing it this
    way is that we can save memory when collapsing, but still be able to call the un-collapsed matrices when needed.

    :param variant_list: DataFrame containing variants to extract.
    :param bgen_path: Path to the BGEN file.
    :param sample_path: Path to the sample file.
    :param should_collapse_matrix: If True (default) collapse variants - sum variants per gene.
        If False, don't sum variants per gene and return the matrix instead.
    :return: A tuple containing:
        - csr_matrix: A sparse matrix with shape (n_samples, n_genes_or_variants).
        - summary_dict: A dictionary with summary information for each gene.
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
        variant_n = 0

        # iterate through each gene in our search list
        for gene_n, current_gene in enumerate(search_list.itertuples()):

            # Fetch actual data from the BGEN file
            variants = bgen_reader.fetch(current_gene.CHROM, current_gene.MIN, current_gene.MAX)

            # create a store for the variant level information
            variant_arrays = []

            gene_variant_start = variant_n
            # collect genotype arrays for each variant
            for current_variant in variants:

                if current_variant.rsid in current_gene.VARS:
                    variant_n += 1
                    # pull out the actual genotypes
                    current_probabilities = current_variant.probabilities

                    # store variant codings
                    variant_array = np.where(current_probabilities[:, 1] == 1, 1,
                                             np.where(current_probabilities[:, 2] == 1, 2, 0))

                    # store variant level information in the array we created
                    variant_arrays.append(variant_array)

            # stack the variant information for all variants in the gene
            stacked_variants = np.column_stack(variant_arrays)

            # if we are collapsing here (to save on memory), leave as False. If set to True, we won't collapse
            # and instead the uncollapsed stacked variants will be appended
            # note the vector naming convention in this small section is a bit hacky but we want the vectors naming
            # to be consistent so it works for the rest of the function
            if should_collapse_matrix is True:
                stacked_variants = stacked_variants.sum(axis=1)
                search_list = gene_n + 1
                gene_n = [gene_n]
            else:
                gene_n = [var for var in range (gene_variant_start, variant_n)]
                search_list = variant_n

            # append the variant arrays to the genotype array
            genotype_arrays.append(stacked_variants)

            # record the per-gene stats in a dict
            summary_dict[current_gene.Index] = {
                'sum': np.sum(stacked_variants),  # use np.sum to get total of all values
                'variants_per_gene': len(variant_arrays),  # get number of variants
                'gene_index': gene_n,
            }

        # stack all genotype arrays into a matrix (samples × variants)
        final_genotypes = np.column_stack(genotype_arrays)

        # convert this to a csr matrix
        final_genotypes = csr_matrix(final_genotypes, shape=(len(current_samples), search_list))

    return final_genotypes, summary_dict
