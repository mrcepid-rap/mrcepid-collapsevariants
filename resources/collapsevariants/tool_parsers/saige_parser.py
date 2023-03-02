import csv
from pathlib import Path
from typing import Tuple

from general_utilities.association_resources import run_cmd


class SAIGEParser:

    def __init__(self, file_prefix: str, chromosome: str):

        self.genes, self.snp_gene_map = self.parse_filters_SAIGE(file_prefix, chromosome)

    # Generate input format files that can be provided to SAIGE
    @staticmethod
    def parse_filters_SAIGE(file_prefix: str, chromosome: str) -> Tuple[dict, dict]:

        # Easier to run SAIGE with a BCF file as I already have that pipeline set up
        cmd = f'plink2 --memory 9000 --threads 1 --bgen /test/{file_prefix}.{chromosome}.bgen \'ref-last\' ' \
              f'--sample /test/{file_prefix}.{chromosome}.sample ' \
              f'--export bcf ' \
              f'--out /test/{file_prefix}.{chromosome}.SAIGE'
        run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting:latest')

        # and index...
        cmd = f'bcftools index /test/{file_prefix}.{chromosome}.SAIGE.bcf'
        run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting:latest')

        # Need to make the SAIGE groupFile. I can use the file 'snp_ENST.txt' created above to generate this...
        with Path('snp_ENST.txt').open('r') as snp_reader,\
                Path(f'{file_prefix}.{chromosome}.SAIGE.groupFile.txt').open('w') as output_setfile_SAIGE:

            genes = {}
            snp_gene_map = {}
            snp_csv = csv.DictReader(snp_reader, delimiter='\t')
            for snp in snp_csv:
                if snp['CHROM'] == chromosome:
                    snp_gene_map[snp['varID']] = snp['ENST']
                    if snp['ENST'] in genes:
                        genes[snp['ENST']]['varIDs'].append(snp['varID'])
                        genes[snp['ENST']]['poss'].append(int(snp['POS']))
                    else:
                        genes[snp['ENST']] = {'CHROM': snp['CHROM'], 'poss': [int(snp['POS'])], 'varIDs': [snp['varID']]}

            for gene in genes:
                min_pos = min(genes[gene]['poss'])
                genes[gene]['min_poss'] = min_pos

            for gene in genes:
                id_string = "\t".join(["{0}:{1}_{2}/{3}".format(*item2) for item2 in [item.split(":") for item in genes[gene]['varIDs']]])
                output_setfile_SAIGE.write(f'{gene}\t{id_string}\n')

            Path(f'{file_prefix}.{chromosome}.SAIGE.log').unlink()

        print(snp_gene_map)
        return genes, snp_gene_map
