import shutil
from pathlib import Path

from general_utilities.association_resources import run_cmd

from tool_parsers.bolt_parser import BOLTParser
from tool_parsers.staar_parser import STAARParser


class SNPMerger:

    def __init__(self, valid_chromosomes: list, file_prefix: str, found_genes: bool):

        self._merge_across_snplist(valid_chromosomes, file_prefix, found_genes)

    @staticmethod
    # This is used regardless of if a SNP of Gene list (apologies for naming...)
    # Note that the default is to merge a SNP list. An extra flag for found genes can be added to do a gene list
    def _merge_across_snplist(valid_chromosomes: list, file_prefix: str, found_genes: bool = False):

        # file prefix to distinguish output from gene list or snp list
        if found_genes:
            tar_type = 'GENE'
        else:
            tar_type = 'SNP'

        # First generate a list of vcfs that we are going to mash together and get variant names:
        cmd = f'bcftools concat -Ob -o /test/{file_prefix}.{tar_type}.SAIGE.pre.bcf '
        variant_ids = ['ENST99999999999' if found_genes else 'ENST00000000000']  # 1st ENST replacement
        snp_gene_map = {}
        for chrom in valid_chromosomes:
            cmd += f'/test/{file_prefix}.{chrom}.SAIGE.bcf '
            with Path(f'{file_prefix}.{chrom}.SAIGE.groupFile.txt').open('r') as group_file:
                for line in group_file:
                    line = line.rstrip()
                    variants = line.split('\t')
                    variant_ids.extend(variants[1:len(variants)])
                    for variant in variants[1:len(variants)]:
                        bolt_format_id = variant.replace('_', ':').replace('/', ':')
                        snp_gene_map[bolt_format_id] = 'ENST00000000000'
                        if found_genes:  # 2nd ENST replacement
                            snp_gene_map[bolt_format_id] = 'ENST99999999999'

        # Combine with bcftools concat
        run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting:latest')
        # Make sure sorted properly...
        cmd = f'bcftools sort -Ob -o /test/{file_prefix}.{tar_type}.SAIGE.bcf ' \
              f'/test/{file_prefix}.{tar_type}.SAIGE.pre.bcf'
        run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting:latest')
        Path(f'{file_prefix}.{tar_type}.SAIGE.pre.bcf').unlink()
        # And index:
        cmd = f'bcftools index /test/{file_prefix}.{tar_type}.SAIGE.bcf'
        run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting:latest')

        # Write new groupFile:
        with Path(f'{file_prefix}.{tar_type}.SAIGE.groupFile.txt').open('w') as snp_groupfile:
            snp_groupfile.write("\t".join(variant_ids))

        # Trick the already made BOLT code above to build a new merged BOLT file:
        genes = {}
        if found_genes:
            genes['ENST99999999999'] = {'CHROM': 1, 'min_poss': 1}
        else:
            genes['ENST00000000000'] = {'CHROM': 1, 'min_poss': 1}

        # This copy is slightly dodgy, as it assumes at least one chromosome has come through in the variable
        # 'chrom' from the loop above
        Path(f'{file_prefix}.{chrom}.sample').rename(f'{file_prefix}.{tar_type}.sample')
        BOLTParser(file_prefix, tar_type, genes, snp_gene_map)

        # Trick the already made STAAR code above to build a new merged set of STAAR files
        STAARParser(file_prefix, tar_type)

        # Delete old files to avoid confusion:
        for chrom in valid_chromosomes:
            Path(f'{file_prefix}.{chrom}.SAIGE.bcf').unlink()
            Path(f'{file_prefix}.{chrom}.SAIGE.bcf.csi').unlink()
            Path(f'{file_prefix}.{chrom}.SAIGE.groupFile.txt').unlink()
            Path(f'{file_prefix}.{chrom}.BOLT.bgen').unlink()
            Path(f'{file_prefix}.{chrom}.BOLT.sample').unlink()
            Path(f'{file_prefix}.{chrom}.STAAR.matrix.rds').unlink()
            Path(f'{file_prefix}.{chrom}.variants_table.STAAR.tsv').unlink()
