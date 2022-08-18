import os
import csv
from typing import Tuple

from ..collapse_resources import *

class SAIGEParser:

    def __init__(self, file_prefix: str, chromosome: str):

        self.genes, self.snp_gene_map = self.parse_filters_SAIGE(file_prefix, chromosome)


    # Generate input format files that can be provided to SAIGE
    @staticmethod
    def parse_filters_SAIGE(file_prefix: str, chromosome: str) -> Tuple[dict, dict]:

        # Easier to run SAIGE with a BCF file as I already have that pipeline set up
        cmd = "plink2 --memory 9000 --threads 1 --bgen /test/" + file_prefix + "." + chromosome + ".bgen 'ref-last' " \
                     "--sample /test/" + file_prefix + "." + chromosome + ".sample " \
                     "--export bcf " \
                     "--out /test/" + file_prefix + "." + chromosome + ".SAIGE"
        run_cmd(cmd, True)

        # and index...
        cmd = "bcftools index /test/" + file_prefix + "." + chromosome + ".SAIGE.bcf"
        run_cmd(cmd, True)

        # Need to make the SAIGE groupFile. I can use the file 'snp_ENST.txt' created above to generate this...
        genes = {}
        snp_gene_map = {}
        snp_reader = csv.DictReader(open('snp_ENST.txt', 'r'), delimiter = "\t")
        for snp in snp_reader:
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

        with open(file_prefix + "." +  chromosome + '.SAIGE.groupFile.txt', 'w') as output_setfile_SAIGE:
            for gene in genes:
                id_string = "\t".join(["{0}:{1}_{2}/{3}".format(*item2) for item2 in [item.split(":") for item in genes[gene]['varIDs']]])
                output_setfile_SAIGE.write("%s\t%s\n" % (gene, id_string))
            output_setfile_SAIGE.close()

        os.remove(file_prefix + "." + chromosome + ".SAIGE.log")

        return genes, snp_gene_map
