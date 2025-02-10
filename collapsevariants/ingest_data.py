import csv
from pathlib import Path
from typing import TypedDict, Dict, Optional, IO, Tuple

import dxpy
from general_utilities.association_resources import download_dxfile_by_name
from general_utilities.mrc_logger import MRCLogger

from collapsevariants.collapse_utils import get_sample_ids


class BGENIndex(TypedDict):
    """A TypedDict containing information on a downloaded bgen file to make for more obvious return types

    :cvar index: A DXFile representation of the .bgen.bgi index file
    :cvar sample: A DXFile representation of the bgen .sample file
    :cvar bgen: A DXFile representation of the .bgen genetic data file
    :cvar vep: A DXFile representation of the annotation .tsv.gz file generated by VEP
    """

    index: str
    sample: str
    bgen: str
    vep: str


class IngestData:
    """A helper class to ingest data to a DNANexus AWS instance

    This class will download and (in some cases) perform pre-processing of the downloaded files to conform to the
    requirements of downstream processing.

    Note that this class is untested as it requires a DNANexus connection to test.

    :param filtering_expression: A string filtering expression to filter variants (must be compatible with pandas.query())
    :param bgen_index: A file containing information of bgen files to collapse on
    :param snp_list: A DXFile containing a list of varIDs to use as a custom mask
    :param gene_list: A DXFile containing a list of gene symbols to collapse into a custom mask
    """

    def __init__(self, bgen_index: dict, filtering_expression: str, snp_list: Optional[dict],
                 gene_list: Optional[dict]):

        # Instantiate the MRC logger
        self._logger = MRCLogger(__name__).get_logger()

        self.filtering_expression = filtering_expression

        self.sample_ids = []
        self.bgen_index = self._ingest_bgen_index(bgen_index)
        self.vep_dict = self._open_vep_filepaths()
        self.snp_list_path = self._define_filter_list(snp_list)
        self.gene_list_path = self._define_filter_list(gene_list)

    def _ingest_bgen_index(self, bgen_index: dict) -> Dict[str, BGENIndex]:
        """Index filtered bgen files from the mrc filtering & annotation workflow

        This class will NOT download the larger genetic data but only the vep information. bgen download is handled
        later in this workflow so that it can be parallelized. This information is stored in a TypedDict (BGENIndex)
        for easier key access.

        :param bgen_index: A file containing information of bgen files to collapse on
        :return: A dictionary with keys equal to chromosome (with 'chr' prefix) and keys of BGENIndex
        """

        bgen_dict = {}
        bgen_index = download_dxfile_by_name(bgen_index)
        samples_collected = False

        # and load it into a dict:
        with bgen_index.open('r') as bgen_reader:
            bgen_index_csv = csv.DictReader(bgen_reader, delimiter='\t')

            for batch in bgen_index_csv:

                bgen_dict[batch['prefix']] = {'index': batch['bgen_index_dxid'], 'sample': batch['sample_dxid'],
                                              'bgen': batch['bgen_dxid'], 'vep': batch['vep_dxid']}

                # Collect sample IDs from the 1st bgen encountered. This is done on the premise that all bgen files
                # have the same sample IDs...
                if samples_collected is False:
                    sample_file = download_dxfile_by_name(batch['sample_dxid'], print_status=False)
                    self.sample_ids = get_sample_ids(sample_file)
                    samples_collected = True
                    sample_file.unlink()  # Make sure to delete to avoid conflict with later download

        if len(self.sample_ids) == 0:
            raise ValueError('No sample IDs found in the bgen files. Please check the sample files and try again.')

        return bgen_dict

    def _open_vep_filepaths(self) -> Dict[str, IO]:

        vep_dict = {}
        for bgen_prefix, bgen_info in self.bgen_index.items():
            vep_dict[bgen_prefix] = dxpy.open_dxfile(bgen_info['vep'], mode='rb')

        return vep_dict

    @staticmethod
    def _define_filter_list(filtering_list: dict) -> Optional[Path]:
        """Download a filtering list file (if provided)

        :param filtering_list: A DXFile ID pointing to some filtering file (likely Gene / SNP-based) on the RAP
        :return: A Path object pointing to the downloaded file, if found, otherwise None
        """

        if filtering_list:
            return download_dxfile_by_name(filtering_list)
        else:
            return None


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
