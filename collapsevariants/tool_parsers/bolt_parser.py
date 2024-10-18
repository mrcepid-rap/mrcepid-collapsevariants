import shutil
import pandas as pd
import numpy as np

from pathlib import Path
from typing import Tuple, Dict, List

from bgen import BgenReader, BgenWriter

from collapsevariants.ingest_data import BGENIndex
from general_utilities.association_resources import download_dxfile_by_name
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.mrc_logger import MRCLogger

LOGGER = MRCLogger(__name__).get_logger()


class BOLTParser:

    def __init__(self, genes: Dict[str, pd.DataFrame], bgen_index: Dict[str, BGENIndex],
                 output_prefix: str):
        """A class to convert BGEN files into BOLT-formatted BGEN files.

        This class functions in two steps, parallelizing the conversion of individual BGEN files into BOLT format for
        all BGEN files provided in the :param genes: / :param bgen_index: dictionaries:

        1. Download the BGEN files from DNANexus using :func:`_download_bgen`

        2. Convert the BGEN files into BOLT format using :func:`_make_bolt_bgen`

        :param genes: A dictionary containing the genes to collapse with keys of the BGEN file prefixes and values of
            Pandas DataFrames for the respective BGEN file containing the list of variants to collapse with columns of
            varID, ENST, CHROM, POS.
        :param bgen_index: A dictionary containing values of BGEN files (:func:`BGENIndex`) with keys of each BGEN file
            prefix.
        :param output_prefix: A string representing the prefix of the output BGEN files.
        """

        self._logger = MRCLogger(__name__).get_logger()
        self._genes = genes
        self._bgen_index = bgen_index
        self._output_prefix = output_prefix

        # Three steps:
        # 1. Next we need to filter the BGEN files according to the SNP list that we generated in snp_list_generator
        thread_utility = ThreadUtility(error_message='Error in filtering BGEN files')
        for bgen_prefix in self._genes:
            thread_utility.launch_job(self._convert_to_bolt_bgen,
                                      bgen_prefix=bgen_prefix)

        self._bolt_output = []
        for result in thread_utility:
            self._bolt_output.append(result)
        thread_utility.collect_futures()

    def get_bolt_output(self) -> List[Tuple[Path, Path, Path]]:
        """Getter for the file-based output of this class"""
        return self._bolt_output

    def _convert_to_bolt_bgen(self, bgen_prefix: str) -> Tuple[Path, Path, Path]:
        """Wrapper method to parallelize the conversion of individual BGEN files into BOLT format genotypes

        :param bgen_prefix: A string representing the prefix of a BGEN file to convert
        :return: A Tuple containing the paths to the BOLT-formatted BGEN file, index file, and sample file
        """

        chrom_bgen_index = self._bgen_index[bgen_prefix],
        variant_list = self._genes[bgen_prefix]

        # Note that index is required but is not explicitly taken as input by BgenReader. It MUST have the same
        # name as the bgen file, but with a .bgi suffix.
        bgen_path, index_path, sample_path = self._download_bgen(chrom_bgen_index)
        bolt_bgen, bolt_index, bolt_sample = self._make_bolt_bgen(bgen_prefix, bgen_path, sample_path, variant_list)
        bgen_path.unlink()
        index_path.unlink()
        sample_path.unlink()

        return bolt_bgen, bolt_index, bolt_sample

    @staticmethod
    def _download_bgen(chrom_bgen_index: BGENIndex) -> Tuple[Path, Path, Path]:
        """Download the BGEN file from DNANexus

        This method downloads the BGEN file from DNANexus using the bgen_index dictionary and the bgen_prefix. It then
        returns the path to the downloaded BGEN file, index, and sample.

        :return: A Tuple containing the paths to the BGEN file, index file, and sample file
        """

        # Download the requisite files for this chromosome according to the index dict:
        bgen_path = download_dxfile_by_name(chrom_bgen_index['bgen'], print_status=False)
        index_path = download_dxfile_by_name(chrom_bgen_index['index'], print_status=False)
        sample_path = download_dxfile_by_name(chrom_bgen_index['sample'], print_status=False)

        return bgen_path, index_path, sample_path

    def _make_bolt_bgen(self, bgen_prefix: str, bgen_path: Path, sample_path: Path,
                        variant_list: pd.DataFrame) -> Tuple[Path, Path, Path]:
        """Generate input format files that can be provided to BOLT

        BOLT only runs on individual variants and has no built-in functionality to perform burden testing. To run BOLT
        in an approximation of 'burden' mode, we have to trick it by first collapsing all variants found within a gene across
        a set of individuals into a single 'variant' that contains information for all individuals for that gene. Thus, we
        create what amounts to a binary yes/no found a variant in the listed gene for a given individual.

        Gene-Variants are encoded as heterozygous, with all genes having a REF allele of A and alternate allele of C. All
        variant IDs are ENST IDs for MANE transcripts (where possible). The actual output of this code is not returned,
        but is spilled to disk as a v1.2, ref-last, 8-bit .bgen format file with a partnered .sample file.

        :param bgen_prefix: A string representing the prefix of a BGEN file to convert
        :param bgen_path: A Path object pointing to the BGEN file to convert
        :param sample_path: A Path object pointing to the sample file for the BGEN file
        :param variant_list: A Pandas DataFrame containing the list of variants to collapse with columns of varID,
            ENST, CHROM, POS
        :return: A Tuple containing the paths to the BOLT-formatted BGEN file, index file, and sample file
        """

        # First aggregate across genes to generate a list of genes and their respective variants
        search_list = variant_list.groupby('ENST').aggregate(
            CHROM=('CHROM', 'first'),
            MIN=('POS', 'min'),
            MAX=('POS', 'max'),
            VARS=('varID', set)
        )

        # Iterate through the provided .bgen file and collect genotypes for each gene
        with BgenReader(bgen_path, sample_path=sample_path, delay_parsing=True) as bgen_reader:
            current_samples = np.array(bgen_reader.samples)

            gene_arrays = {}

            for current_Gene in search_list.itertuples():

                # Note to future devs: it is MUCH faster to fetch a window of variants and iterate through them, checking if
                # they are in our list of variants, then to iterate through the list of variants and fetch each one individually.
                variants = bgen_reader.fetch(current_Gene.CHROM.replace('chr', ''), current_Gene.MIN, current_Gene.MAX)

                current_ENST = current_Gene.Index
                current_POS = current_Gene.MIN
                current_CHROM = current_Gene.CHROM

                gene_arrays[current_ENST] = {'genotype_boolean': np.zeros(len(current_samples), dtype=bool),
                                             'min_pos': current_POS, 'chrom': current_CHROM}

                for current_variant in variants:
                    if current_variant.rsid in current_Gene.VARS:
                        current_probabilities = current_variant.probabilities

                        # AFAIK genotype probabilities can ONLY be 0 / 1
                        genotype_boolean = np.logical_or(current_probabilities[:, 1] == 1,
                                                         current_probabilities[:, 2] == 1)
                        gene_arrays[current_ENST]['genotype_boolean'] = np.logical_or(
                            gene_arrays[current_ENST]['genotype_boolean'], genotype_boolean)

        bolt_bgen = Path(f'{self._output_prefix}.{bgen_prefix}.BOLT.bgen')

        with BgenWriter(bolt_bgen, n_samples=len(current_samples), samples=current_samples, compression='zstd',
                        layout=2) as bgen_writer:

            prev_pos = 0

            # Should already be sorted by min_pos, but just in case sort when iterating:
            for ENST, gene_info in sorted(gene_arrays.items(), key=lambda value: value[1]['min_pos']):

                if prev_pos > gene_info['min_pos']:
                    raise ValueError('Genes are not sorted by position!')

                prev_pos = gene_info['min_pos']

                genotypes = np.array(
                    [[0., 1., 0.] if genotype == True else [1., 0., 0.] for genotype in gene_info['genotype_boolean']])

                bgen_writer.add_variant(varid=ENST, rsid=ENST, chrom=gene_info['chrom'], pos=gene_info['min_pos'],
                                        alleles=['A', 'C'], genotypes=genotypes, ploidy=2, phased=False, bit_depth=8)

        bolt_index = bolt_bgen.with_suffix('.bgen.bgi')
        bolt_sample = shutil.copy(sample_path, bolt_bgen.with_suffix('.sample'))
        return bolt_bgen, bolt_index, bolt_sample
