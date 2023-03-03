from pathlib import Path
from typing import List, Dict

from general_utilities.association_resources import run_cmd

from tool_parsers.bolt_parser import parse_filters_BOLT
from tool_parsers.saige_parser import GeneDict
from tool_parsers.staar_parser import STAARParser


class SNPMerger:
    """A class to merged SNPs into either a single SNP or Gene-specific set of burden testing files

    :param valid_chromosomes: A list containing chromosomes we actually found data for
    :param file_prefix: A name to append to beginning of output files.
    :param found_genes: Did we find a GeneList to process?
    """

    def __init__(self, valid_chromosomes: List[str], file_prefix: str, found_genes: bool):

        self._valid_chromosomes = valid_chromosomes
        self._file_prefix = file_prefix
        self._found_genes = found_genes

        # file prefix to distinguish output from gene list or snp list and prepping for downstream processing
        # self._genes is full of infomation to fake out the BOLT methods to write a single file
        if self._found_genes:
            self._tar_type = 'GENE'
            self._genes: Dict[str, GeneDict] = {'ENST99999999999': {'CHROM': '1', 'min_poss': 1,
                                                                    'varIDs': [], 'poss': []}}
        else:
            self._tar_type = 'SNP'
            self._genes: Dict[str, GeneDict] = {'ENST00000000000': {'CHROM': '1', 'min_poss': 1,
                                                                    'varIDs': [], 'poss': []}}

        self._merge_across_snplist()

    def _merge_across_snplist(self):
        """Merges all chromosomes into a single set of files

        This is used regardless if either a SNP or Gene list is found (apologies for SNP naming...)
        Note that the default is to merge a SNP list. An extra flag for found genes can be added to do a gene list
        """

        # First generate a list of vcfs that we are going to mash together and get variant names:
        cmd = f'bcftools concat -Ob -o /test/{self._file_prefix}.{self._tar_type}.SAIGE.pre.bcf '
        variant_ids = ['ENST99999999999' if self._found_genes else 'ENST00000000000']  # 1st ENST replacement
        snp_gene_map = {}
        for chrom in self._valid_chromosomes:
            cmd += f'/test/{self._file_prefix}.{chrom}.SAIGE.bcf '
            with Path(f'{self._file_prefix}.{chrom}.SAIGE.groupFile.txt').open('r') as group_file:
                for line in group_file:
                    line = line.rstrip()
                    variants = line.split('\t')
                    variant_ids.extend(variants[1:len(variants)])
                    for variant in variants[1:len(variants)]:
                        bolt_format_id = variant.replace('_', ':').replace('/', ':')
                        snp_gene_map[bolt_format_id] = 'ENST00000000000'
                        if self._found_genes:  # 2nd ENST replacement
                            snp_gene_map[bolt_format_id] = 'ENST99999999999'

        # Combine with bcftools concat
        run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting:latest')

        # Make sure sorted properly...
        cmd = f'bcftools sort -Ob -o /test/{self._file_prefix}.{self._tar_type}.SAIGE.bcf ' \
              f'/test/{self._file_prefix}.{self._tar_type}.SAIGE.pre.bcf'
        run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting:latest')
        Path(f'{self._file_prefix}.{self._tar_type}.SAIGE.pre.bcf').unlink()  # Delete old file

        # And index:
        cmd = f'bcftools index /test/{self._file_prefix}.{self._tar_type}.SAIGE.bcf'
        run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting:latest')

        # Write new groupFile:
        with Path(f'{self._file_prefix}.{self._tar_type}.SAIGE.groupFile.txt').open('w') as snp_groupfile:
            snp_groupfile.write("\t".join(variant_ids))

        # Trick the already made BOLT code above to build a new merged BOLT file:
        # This copy is slightly dodgy, as it assumes at least one chromosome has come through in the variable
        # 'chrom' from the loop above
        Path(f'{self._file_prefix}.{chrom}.sample').rename(f'{self._file_prefix}.{self._tar_type}.sample')
        parse_filters_BOLT(self._file_prefix, self._tar_type, self._genes, snp_gene_map)

        # Trick the already made STAAR code above to build a new merged set of STAAR files
        STAARParser(self._file_prefix, self._tar_type).parse_filters_STAAR()

        # Delete old files to avoid confusion:
        for chrom in self._valid_chromosomes:
            Path(f'{self._file_prefix}.{chrom}.SAIGE.bcf').unlink()
            Path(f'{self._file_prefix}.{chrom}.SAIGE.bcf.csi').unlink()
            Path(f'{self._file_prefix}.{chrom}.SAIGE.groupFile.txt').unlink()
            Path(f'{self._file_prefix}.{chrom}.BOLT.bgen').unlink()
            Path(f'{self._file_prefix}.{chrom}.BOLT.sample').unlink()
            Path(f'{self._file_prefix}.{chrom}.STAAR.matrix.rds').unlink()
            Path(f'{self._file_prefix}.{chrom}.variants_table.STAAR.tsv').unlink()
