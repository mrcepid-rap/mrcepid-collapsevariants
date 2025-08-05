import csv
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
from scipy.sparse import csr_matrix

from collapsevariants.tool_parsers.tool_parser import ToolParser
from collapsevariants.utilities.collapse_utils import GenotypeInfo


class REGENIEParser(ToolParser):
    """Generate files used by REGENIE for association testing. Implements the abstract class ToolParser.

    This class is used to generate files that can be used by REGENIE to perform association testing. REGENIE
    requires 3 files to run: a setlist file, an annotation file, and a mask file. This class generates these files
    sequentially.

    :param genes: A dictionary containing a list of genes and their respective variants.
    :param genotype_index: A dictionary containing a list of genotype matrices for each bgen file.
    :param sample_ids: A list of sample IDs to write to the sample file.
    :param output_prefix: A string representing the prefix to use for the output files.
    """

    def __init__(self, genes: Dict[str, pd.DataFrame], genotype_index: Dict[str, Tuple[csr_matrix, Dict[str, GenotypeInfo]]],
                 sample_ids: List[str], output_prefix: str):
        super().__init__(genes=genes, genotype_index=genotype_index, output_prefix=output_prefix, sample_ids=sample_ids,
                         tool_name='REGENIE')

    def _make_output_files(self, bgen_prefix: str) -> List[Path]:
        """Wrapper method to parallelize the conversion of individual variant data into REGENIE-compatible files.

        :param bgen_prefix: A string representing the prefix of a BGEN file to convert
        :return: A list containing the paths to the REGENIE-formatted setlist, annotation, and mask files
        """

        variant_list = self._genes[bgen_prefix]
        setlist_path = self._make_regenie_setlist_file(bgen_prefix, variant_list)
        annotation_path = self._make_regenie_annotation_file(bgen_prefix, variant_list)
        mask_path = self._make_regenie_mask_file(bgen_prefix)

        return [setlist_path, annotation_path, mask_path]

    def _make_regenie_annotation_file(self, bgen_prefix: str, variant_list: pd.DataFrame) -> Path:
        """Generate the annotation file for REGENIE (*.REGENIE.annotationFile.txt).

        The output of this method is in a tab-delimited format with the following columns:
        - VAR: The variant ID
        - ENST: The ENST ID
        - MASK: The mask name

        :param bgen_prefix: A name to append to beginning of output files.
        :param variant_list: A DataFrame containing the variants to be included in the annotation file.
        :return: A Path object representing the path to the annotation file.
        """

        output_annotation_path = Path(f'{self._output_prefix}.{bgen_prefix}.REGENIE.annotationFile.txt')

        with output_annotation_path.open('w') as output_annotation_writer:
            output_setfile_csv = csv.DictWriter(output_annotation_writer, delimiter='\t',
                                                fieldnames=['VAR', 'ENST', 'MASK'])

            for var_info in variant_list.itertuples():
                out_dict = {'VAR': var_info.varID, 'ENST': var_info.ENST, 'MASK': self._output_prefix}
                output_setfile_csv.writerow(out_dict)

        return output_annotation_path

    @staticmethod
    def _reformat_regenie_variants(var_id_list: List[str]) -> str:
        """Format variant IDs according to SAIGE input specifications

        This method serves as an input function to :func:`pandas.DataFrame.aggregate` to both reformat variant IDs
        and join them into a single (tab-delimited) string for output to a file.

        :param var_id_list: A list of N variant IDs from a single ENST group
        :return: A single tab-delimited string consisting of N variant IDs reformatted for SAIGE input
        """

        joined_reformatted_ids = ','.join(var_id_list)
        return joined_reformatted_ids

    def _make_regenie_setlist_file(self, bgen_prefix: str, variant_list: pd.DataFrame) -> Path:
        """Generate input format file that can be provided to REGENIE (*.REGENIE.setListFile.txt)

        For the way SAIGE is implemented downstream of this applet, SAIGE requires a standard .bcf file and a 'groupFile'.
        This groupFile is a simple tab-delimited file with individual genes (here defined as ENSTs) as the first column,
        followed by all the variant IDs that should be included in that gene per our mask definition.

        Variant IDs are slightly different from that included in other files produced by this applet, in that they must
        follow the format of CHR:POS_REF/ALT rather than CHR_POS_REF_ALT as defined in the original .bgen files produced
        prior to running this applet.

        :param bgen_prefix: A name to append to beginning of output files.
        :return: A tuple containing two dictionaries of ENSTs mapped to variants and variants mapped to ENSTs, respectively
        """

        output_setlist_path = Path(f'{self._output_prefix}.{bgen_prefix}.REGENIE.setListFile.txt')

        with output_setlist_path.open('w') as output_setlist_writer:

            output_group_csv = csv.DictWriter(output_setlist_writer, delimiter='\t',
                                              fieldnames=['ENST', 'CHROM', 'MIN', 'VARS'],
                                              extrasaction='ignore')

            variants_grouped = variant_list.groupby('ENST').aggregate(
                CHROM=('CHROM', 'first'),
                MIN=('POS', 'min'),
                VARS=('varID', self._reformat_regenie_variants)
            )
            variants_grouped = variants_grouped.reset_index()

            prev_pos = 0

            for gene_info in variants_grouped.sort_values(by='MIN').to_dict(orient='index').values():

                if prev_pos > gene_info['MIN']:
                    raise ValueError('Genes are not sorted by position!')

                prev_pos = gene_info['MIN']
                output_group_csv.writerow(gene_info)

        return output_setlist_path

    def _make_regenie_mask_file(self, bgen_prefix: str) -> Path:
        """Create the required mask file for REGENIE (*.REGENIE.maskfile.txt).

        This method creates a mask file for REGENIE, which is a simple tab-delimited file with one row and two columns,
        both with the same value: the output prefix. This is a dummy file that is not used in the downstream analysis,
        but is required by REGENIE to run.

        :param bgen_prefix: A name to append to beginning of output files.
        :return: A Path object representing the path to the mask file.
        """

        output_mask_path = Path(f'{self._output_prefix}.{bgen_prefix}.REGENIE.maskfile.txt')

        with output_mask_path.open('w') as output_mask_writer:
            output_mask_writer.write(f'{self._output_prefix}\t{self._output_prefix}\n')

        return output_mask_path
