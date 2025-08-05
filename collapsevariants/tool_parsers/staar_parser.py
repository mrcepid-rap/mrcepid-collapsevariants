import csv
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
from scipy.sparse import csr_matrix

from collapsevariants.tool_parsers.tool_parser import ToolParser
from utilities.collapse_utils import GenotypeInfo


class STAARParser(ToolParser):

    def __init__(self, genes: Dict[str, pd.DataFrame], genotype_index: Dict[str, Tuple[csr_matrix, Dict[str, GenotypeInfo]]],
                 sample_ids: List[str], output_prefix: str):
        """Generate files used by STAAR for association testing. Implements the abstract class ToolParser.

        This class is used to generate files that can be used by STAAR to perform association testing. STAAR
        requires 2 files to run: a samples table and a variants table. This class generates these files sequentially.

        :param genes: A dictionary containing a list of genes and their respective variants.
        :param genotype_index: A dictionary containing a list of genotype matrices for each bgen file.
        :param sample_ids: A list of sample IDs to write to the sample file.
        :param output_prefix: A string representing the prefix to use for the output files.
        """
        super().__init__(genes=genes, genotype_index=genotype_index, output_prefix=output_prefix, sample_ids=sample_ids,
                         tool_name='STAAR')

    def _make_output_files(self, bgen_prefix: str) -> List[Path]:
        """Wrapper method to parallelize the conversion of individual variant data into STAAR compatible files.

        :param bgen_prefix: A string representing the prefix of a BGEN file to convert
        :return: A list containing the paths to the STAAR samples and variants tables
        """

        variant_list = self._genes[bgen_prefix]

        samples_table_file = self.make_samples_dict(self._output_prefix, bgen_prefix, self._sample_ids)
        variants_table_file = self.make_variants_dict(self._output_prefix, bgen_prefix, variant_list)

        return [samples_table_file, variants_table_file]

    @staticmethod
    def make_samples_dict(output_prefix: str, bgen_prefix: str, sample_ids: List[str]) -> Path:
        """Make a sample 'dictionary' in text format

        Make a text file of sample ID / rowNum for writing the final matrix.

        Note that this is a staticmethod to enable use by GENE / SNP mask generation methods in the main applet. This
        is not my favourite way to do this, but it is the most expedient way to do so.

        :param output_prefix: The prefix to use for the output files
        :param bgen_prefix: The prefix of the BGEN file being processed
        :param sample_ids: A list of sample IDs to write to the file
        :return: A PathLike for the sample file in STAAR-ready format
        """

        samples_table_path = Path(f'{output_prefix}.{bgen_prefix}.STAAR.samples_table.tsv')
        with samples_table_path.open('w') as samples_dict_writer:
            samples_dict_writer.write('sampID\trow\n')

            for col_num, sample in enumerate(sample_ids):  # Don't use enumerate due to the 2nd header row!
                samples_dict_writer.write(f'{sample}\t{col_num}\n')

            samples_dict_writer.close()

        return samples_table_path

    @staticmethod
    def make_variants_dict(output_prefix: str, bgen_prefix: str, variant_list: pd.DataFrame) -> Path:
        """Make a variant 'dictionary' in text format

        Makes a text file of variant ID / chrom / pos / ENST / column for use in making the final matrix .tsv but also
        as a dictionary of possible variants for downstream processing both by STAAR and other tools.

        Note that this is a staticmethod to enable use by GENE / SNP mask generation methods in the main applet. This
        is not my favourite way to do this, but it is the most expedient way to do so.

        :param output_prefix: The prefix to use for the output files
        :param bgen_prefix: The prefix of the BGEN file being processed
        :param variant_list: A pandas DataFrame containing the variant information
        :return: A Tuple containing a PathLike for the variant dictionary and a python dict() representation of this
            data
        """

        variants_table_file = Path(f'{output_prefix}.{bgen_prefix}.STAAR.variants_table.tsv')
        with variants_table_file.open('w') as variants_dict_writer:
            variants_dict_csv = csv.DictWriter(variants_dict_writer, delimiter='\t', extrasaction='ignore',
                                               fieldnames=['varID', 'CHROM', 'POS', 'ENST', 'column'])
            variants_dict_csv.writeheader()

            for column_num, var in enumerate(variant_list.sort_values(by='POS').itertuples()):
                out_dict = {'varID': var.varID, 'CHROM': var.CHROM, 'POS': var.POS, 'ENST': var.ENST,
                            'column': column_num}
                variants_dict_csv.writerow(out_dict)

        return variants_table_file
