from pathlib import Path
from typing import Tuple, Dict, List, Union

import numpy as np
import pandas as pd
from bgen import BgenWriter
from general_utilities.mrc_logger import MRCLogger
from scipy.sparse import csr_matrix

from collapsevariants.tool_parsers.tool_parser import ToolParser

LOGGER = MRCLogger(__name__).get_logger()


class BOLTParser(ToolParser):
    """Generate files used by BOLT for association testing. Implements the abstract class ToolParser.

    This class is used to generate files that can be used by BOLT-LMM to perform association testing. BOLT-LMM
    requires 3 files to run: a .bgen file, a .bgen.bgi file, and a .sample file. This class generates these files
    sequentially.

    :param genes: A dictionary containing a list of genes and their respective variants.
    :param genotype_index: A dictionary containing a list of genotype matrices for each bgen file.
    :param sample_ids: A list of sample IDs to write to the sample file.
    :param output_prefix: A string representing the prefix to use for the output files.
    """

    def __init__(self, genes: Dict[str, pd.DataFrame], genotype_index: Dict[str, csr_matrix], sample_ids: List[str],
                 output_prefix: str):
        super().__init__(genes=genes, genotype_index=genotype_index, output_prefix=output_prefix, sample_ids=sample_ids,
                         tool_name='BOLT')

    def _make_output_files(self, bgen_prefix: str) -> List[Path]:
        """Wrapper method to parallelize the conversion of individual BGEN files into BOLT format genotypes

        :param bgen_prefix: A string representing the prefix of a BGEN file to convert
        :return: A list containing the paths to the BOLT-formatted BGEN file, index file, and sample file
        """

        variant_list = self._genes[bgen_prefix]
        genotypes = self._genotype_index[bgen_prefix]

        bolt_matrix = self._make_bolt_matrix(genotypes, variant_list)
        bolt_bgen, bolt_index = self._write_bolt_bgen(bgen_prefix, bolt_matrix)
        bolt_sample = self._write_bolt_sample(bgen_prefix)

        return [bolt_bgen, bolt_index, bolt_sample]

    @staticmethod
    def _make_bolt_matrix(genotypes: Tuple[csr_matrix, Dict[str, Dict[str, int]]], variant_list: pd.DataFrame) -> Dict[
        str, Dict[str, Union[np.ndarray, int, str]]]:
        """
        Create a dictionary of gene arrays for BOLT association testing.

        This method processes the provided genotype matrix and variant list to generate a dictionary of gene arrays.
        Each gene array contains a boolean array indicating the presence of variants in each sample, along with the
        minimum position and chromosome information for each gene.

        :param genotypes: A dictionary containing the genotype matrix and summary dictionary.
        :param variant_list: A Pandas DataFrame containing the list of variants to be processed.
        :return: A dictionary where keys are ENST IDs and values are dictionaries containing:
            - 'genotype_boolean': A boolean array indicating the presence of variants in each sample.
            - 'min_pos': The minimum position of the variants for the gene.
            - 'chrom': The chromosome of the gene.
        """
        # Iterate through the provided .bgen file and collect genotypes for each gene
        gene_arrays = {}

        # unpack the tuple
        genotype_matrix, summary_dict = genotypes

        search_list = variant_list.groupby('ENST').aggregate(
            CHROM=('CHROM', 'first'),
            MIN=('POS', 'min'))

        # Iterate over each gene in the search list
        for current_gene in search_list.itertuples():
            # Extract the current gene's ENST (Index), position (MIN), and chromosome (CHROM)
            current_enst = current_gene.Index
            current_pos = current_gene.MIN
            current_chrom = current_gene.CHROM

            # Subset the genotype matrix to get the current_gene as defined by the gene index (defined in
            # the generate_csr_matrix_from_bgen function
            current_genotypes = genotype_matrix[:, summary_dict[current_enst]['gene_index']]

            # Sum genotypes across samples to determine the presence of variants
            sample_sums = current_genotypes.sum(axis=1).A1

            # Create a boolean array indicating whether a variant is present in each sample
            sample_booleans = np.where(sample_sums > 0., True, False)

            # Store the processed information in the dictionary
            gene_arrays[current_enst] = {'genotype_boolean': sample_booleans,
                                         'min_pos': current_pos, 'chrom': current_chrom}

        return gene_arrays

    def _write_bolt_bgen(self, bgen_prefix: str, gene_arrays: Dict[str, dict]) -> Tuple[Path, Path]:
        """Generate input format files that can be provided to BOLT.

        BOLT only runs on individual variants and has no built-in functionality to perform burden testing. To run
        BOLT in an approximation of 'burden' mode, we have to trick it by first collapsing all variants found within
        a gene across a set of individuals into a single 'variant' that contains information for all individuals for
        that gene. Thus, we create what amounts to a binary yes/no found a variant in the listed gene for a given
        individual.

        Gene-Variants are encoded as heterozygous, with all genes having a REF allele of A and alternate allele of C.
        All variant IDs are ENST IDs for MANE transcripts (where possible). The actual output of this code is not
        returned, but is spilled to disk as a v1.2, ref-last, 8-bit .bgen format file with a partnered .sample file.

        :param bgen_prefix: A string representing the prefix of a BGEN file to convert
        :param gene_arrays: A dictionary containing a list of genes and their respective genotypes from
            :func:`_make_bolt_matrix`
        :return: A Tuple containing the paths to the BOLT-formatted BGEN file and it's accompanying index file
        """

        bolt_bgen = Path(f'{self._output_prefix}.{bgen_prefix}.BOLT.bgen')

        with BgenWriter(bolt_bgen, n_samples=len(self._sample_ids), samples=self._sample_ids, compression='zstd',
                        layout=2) as bgen_writer:

            prev_pos = 0

            # Should already be sorted by min_pos, but just in case sort when iterating:
            for ENST, gene_info in sorted(gene_arrays.items(), key=lambda value: value[1]['min_pos']):

                if prev_pos > gene_info['min_pos']:
                    raise ValueError('Genes are not sorted by position!')

                prev_pos = gene_info['min_pos']

                genotypes = np.array(
                    [[0., 1., 0.] if genotype == True else [1., 0., 0.] for genotype in gene_info['genotype_boolean']])

                bgen_writer.add_variant(varid=ENST, rsid=ENST, chrom=str(gene_info['chrom']), pos=gene_info['min_pos'],
                                        alleles=['A', 'C'], genotypes=genotypes, ploidy=2, phased=False, bit_depth=8)

        bolt_index = bolt_bgen.with_suffix('.bgen.bgi')
        return bolt_bgen, bolt_index

    def _write_bolt_sample(self, bgen_prefix: str) -> Path:
        """Generate a sample file that can be provided to BOLT.

        This method was created to ensure that this method can be tested independently of DNANexus. This ensures that
        the sample file is created directly from a list of sample IDs, rather than from the (downloaded) sample file
        itself.

        :param bgen_prefix: A string representing the prefix of a BGEN file to convert.
        :return: A Path to the BOLT-formatted sample file
        """

        bolt_sample = Path(f'{self._output_prefix}.{bgen_prefix}.BOLT.sample')

        with bolt_sample.open('w') as bolt_sample_writer:
            bolt_sample_writer.write('ID_1 ID_2 missing sex\n')
            bolt_sample_writer.write('0 0 0 D\n')
            for sample in self._sample_ids:
                bolt_sample_writer.write(f'{sample} {sample} 0 NA\n')

        return bolt_sample
