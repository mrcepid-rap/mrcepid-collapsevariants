import pandas as pd

from pathlib import Path
from typing import Dict, List
from scipy.sparse import csr_matrix

from collapsevariants.tool_parsers.tool_parser import ToolParser


class SAIGEParser(ToolParser):

    def __init__(self, genes: Dict[str, pd.DataFrame], genotype_index: Dict[str, csr_matrix], sample_ids: List[str],
                 output_prefix: str):
        super().__init__(genes=genes, genotype_index=genotype_index, output_prefix=output_prefix, sample_ids=sample_ids,
                         tool_name='SAIGE')

    def _make_output_files(self, bgen_prefix: str) -> List[Path]:

        variant_list = self._genes[bgen_prefix]

        return [self._make_saige_group_file(bgen_prefix, variant_list)]

    @staticmethod
    def _reformat_saige_variants(var_id_list: List[str]) -> str:
        """Format variant IDs according to SAIGE input specifications

        This method serves as an input function to :func:`pandas.DataFrame.aggregate` to both reformat variant IDs
        and join them into a single (tab-delimited) string for output to a file.

        :param var_id_list: A list of N variant IDs from a single ENST group
        :return: A single tab-delimited string consisting of N variant IDs reformatted for SAIGE input
        """

        reformatted_ids = []
        for varID in var_id_list:

            split_id = varID.split("_")
            reformatted_id = "{0}:{1}_{2}/{3}".format(*split_id)
            reformatted_ids.append(reformatted_id)

        joined_reformatted_ids = "\t".join(reformatted_ids)
        return joined_reformatted_ids

    def _make_saige_group_file(self, bgen_prefix: str, variant_list: pd.DataFrame) -> Path:
        """Generate input format file that can be provided to SAIGE

        For the way SAIGE is implemented downstream of this applet, SAIGE requires a standard .bcf file and a 'groupFile'.
        This groupFile is a simple tab-delimited file with individual genes (here defined as ENSTs) as the first column,
        followed by all the variant IDs that should be included in that gene per our mask definition.

        Variant IDs are slightly different from that included in other files produced by this applet, in that they must
        follow the format of CHR:POS_REF/ALT rather than CHR_POS_REF_ALT as defined in the original .bgen files produced
        prior to running this applet.

        :param bgen_prefix: A name to append to beginning of output files.
        :return: A tuple containing two dictionaries of ENSTs mapped to variants and variants mapped to ENSTs, respectively
        """

        output_group_file = Path(f'{self._output_prefix}.{bgen_prefix}.SAIGE.groupFile.txt')

        with output_group_file.open('w') as output_setfile_SAIGE:

            variants_grouped = variant_list.groupby('ENST').aggregate(
                MIN=('POS', 'min'),
                VARS=('varID', self._reformat_saige_variants)
            )
            variants_grouped = variants_grouped.reset_index()

            prev_pos = 0

            for gene_info in variants_grouped.sort_values(by='MIN').to_dict(orient='index').values():

                if prev_pos > gene_info['MIN']:
                    raise ValueError('Genes are not sorted by position!')

                prev_pos = gene_info['MIN']
                output_setfile_SAIGE.write(f'{gene_info["ENST"]}\t{gene_info["VARS"]}\n')

        return output_group_file
