import os
import csv
import dxpy

from pathlib import Path
from typing import TypedDict, Dict

from general_utilities.association_resources import run_cmd


# A TypedDict to make for more obvious return types
class BGENIndex(TypedDict):
    index: str
    sample: str
    bgen: str
    vep: str


class IngestData:

    def __init__(self, filtering_expression: str, bgen_index: dict, snplist: dict, genelist: dict):

        self.filtering_expression = filtering_expression

        self._ingest_docker()
        self._set_num_threads()
        self.bgen_index = self._ingest_bgen_index(bgen_index)
        self.found_snps = self._define_snplist(snplist)
        self.found_genes = self._define_genelist(genelist)

    # Grab our docker image
    @staticmethod
    def _ingest_docker():

        cmd = 'docker pull egardner413/mrcepid-burdentesting:latest'
        run_cmd(cmd, is_docker=False)
        print('Docker loaded')

    def _set_num_threads(self):
        self.threads = os.cpu_count()
        print(f'Number of threads available: {self.threads}')

    # Ingest the INDEX of bgen files and download VEP indices
    @staticmethod
    def _ingest_bgen_index(bgen_index: dict) -> Dict[str, BGENIndex]:

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

        found_snps = False
        if snplist is not None:
            snplist = dxpy.DXFile(snplist)
            dxpy.download_dxfile(snplist, 'snp_list.snps')
            found_snps = True

        return found_snps

    @staticmethod
    def _define_genelist(genelist: dict) -> bool:

        found_genes = False
        if genelist is not None:
            genelist = dxpy.DXFile(genelist)
            dxpy.download_dxfile(genelist, 'gene_list.genes')
            found_genes = True

        return found_genes

    def _check_filtering_expression(self):
        if self.filtering_expression is not None and self.found_snps is False and self.found_genes is False:
            # For logging purposes output the filtering expression provided by the user
            print('Current Filtering Expression:')
            print(self.filtering_expression)
        elif self.filtering_expression is not None and self.found_snps is False and self.found_genes:
            print('Gene list provided together with filtering expression - will filter for variant classes of interest in gene set')
            print('Current Filtering Expression:')
            print(self.filtering_expression)
        elif self.filtering_expression is None and self.found_snps and self.found_genes is False:
            print('Using provided SNPlist to generate a mask...')
        else:
            raise dxpy.AppError('Incorrect input provided...quitting!')


