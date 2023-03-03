import csv
import dxpy

from pathlib import Path
from typing import TypedDict, Dict

from general_utilities.association_resources import run_cmd
from general_utilities.mrc_logger import MRCLogger


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

    :param filtering_expression: A string filtering expression to filter variants (must be compatible with
    pandas.query())
    :param bgen_index: A file containing information of bgen files to collapse on
    :param snplist: A DXFile containing a list of varIDs to use as a custom mask
    :param genelist: A DXFile containing a list of gene symbols to collapse into a custom mask
    """

    def __init__(self, filtering_expression: str, bgen_index: dict, snplist: dict, genelist: dict):

        # Instantiate the MRC logger
        self._logger = MRCLogger(__name__).get_logger()

        self.filtering_expression = filtering_expression

        self._ingest_docker()
        self.bgen_index = self._ingest_bgen_index(bgen_index)
        self.found_snps = self._define_snplist(snplist)
        self.found_genes = self._define_genelist(genelist)

    # Grab our docker image
    def _ingest_docker(self):
        """Download the standard mrcepid docker instance to this instance

        This method just runs `docker pull egardner413/mrcepid-burdentesting:latest`
        """

        cmd = 'docker pull egardner413/mrcepid-burdentesting:latest'
        run_cmd(cmd, is_docker=False)
        self._logger.info('Docker loaded')

    # Ingest the INDEX of bgen files and download VEP indices
    @staticmethod
    def _ingest_bgen_index(bgen_index: dict) -> Dict[str, BGENIndex]:
        """Index filtered bgen files from the mrc filtering & annotation workflow

        This class will NOT download the larger genetic data but only the vep information. bgen download is handled
        later in this workflow so that it can be parallelized. This information is stored in a TypedDict (BGENIndex)
        for easier key access.

        :param bgen_index: A file containing information of bgen files to collapse on
        :return: A dictionary with keys equal to chromosome (with 'chr' prefix) and keys of BGENIndex
        """

        bgen_dict = {}
        bgen_index = dxpy.DXFile(bgen_index)
        dxpy.download_dxfile(bgen_index.get_id(), 'bgen_locs.tsv')

        # and load it into a dict:
        Path('filtered_bgen/').mkdir()  # For downloading later...
        with Path('bgen_locs.tsv').open('r') as bgen_reader:
            bgen_index_csv = csv.DictReader(bgen_reader, delimiter='\t')

            for chrom in bgen_index_csv:
                bgen_dict[chrom['chrom']] = {'index': chrom['bgen_index_dxid'], 'sample': chrom['sample_dxid'],
                                             'bgen': chrom['bgen_dxid'], 'vep': chrom['vep_dxid']}

                # but download the vep index because of how we generate the SNP list:
                vep = dxpy.DXFile(chrom['vep_dxid'])
                dxpy.download_dxfile(vep.get_id(), f'filtered_bgen/chr{chrom["chrom"]}.filtered.vep.tsv.gz')

        return bgen_dict

    @staticmethod
    def _define_snplist(snplist: dict) -> bool:
        """Download the SNPList file (if provided)

        :param snplist: A DXFile ID pointing to the SNPList file on the RAP
        :return: A boolean defining if a SNPList file was found
        """

        found_snps = False
        if snplist is not None:
            snplist = dxpy.DXFile(snplist)
            dxpy.download_dxfile(snplist, 'snp_list.snps')
            found_snps = True

        return found_snps

    @staticmethod
    def _define_genelist(genelist: dict) -> bool:
        """Download the SNPList file (if provided)

        :param genelist: A DXFile ID pointing to the GeneList file on the RAP
        :return: A boolean defining if a GeneList file was found
        """

        found_genes = False
        if genelist is not None:
            genelist = dxpy.DXFile(genelist)
            dxpy.download_dxfile(genelist, 'gene_list.genes')
            found_genes = True

        return found_genes

    def _check_filtering_expression(self) -> None:
        """Check to make sure the filtering expression, SNPFile, and GeneFile is provided in the correct combination

        We check to make sure that we fit one of three filtering styles:

        1. Filtering expression alone
        2. Filtering expression with a Gene List (for making a gene collapsing mask)
        3. A SNP List alone (for making a SNP collapsed mask)
        """
        if self.filtering_expression is not None and self.found_snps is False and self.found_genes is False:
            # For logging purposes output the filtering expression provided by the user
            self._logger.info('Current Filtering Expression:')
            self._logger.info(self.filtering_expression)
        elif self.filtering_expression is not None and self.found_snps is False and self.found_genes:
            self._logger.info('Gene list provided together with filtering expression - will filter for variant '
                              'classes of interest in gene set')
            self._logger.info('Current Filtering Expression:')
            self._logger.info(self.filtering_expression)
        elif self.filtering_expression is None and self.found_snps and self.found_genes is False:
            self._logger.info('Using provided SNPlist to generate a mask...')
        else:
            raise ValueError('Incorrect input for snp/gene/filtering expression provided... exiting!')
