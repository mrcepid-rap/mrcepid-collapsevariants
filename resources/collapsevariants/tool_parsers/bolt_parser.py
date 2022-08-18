import os
import csv
from typing import Tuple
import pandas as pd
import numpy as np

from ..collapse_resources import *


class BOLTParser:

    def __init__(self, file_prefix: str, chromosome: str, genes: dict, snp_gene_map: dict):

        self.poss_indiv, self.samples = self._parse_filters_BOLT(file_prefix, chromosome, genes, snp_gene_map)
        self.sample_table = self._check_vcf_stats(self.poss_indv, self.samples)

    # Generate input format files that can be provided to BOLT
    @staticmethod
    def _parse_filters_BOLT(file_prefix: str, chromosome: str, genes: dict, snp_gene_map: dict) -> Tuple[list, dict]:

        # Get out BCF file into a tab-delimited format that we can parse for BOLT.
        # We ONLY want alternate alleles here (-i 'GT="alt") and then for each row print:
        # 1. Sample ID: UKBB eid format
        # 2. The actual genotype (0/1 or 1/1 in this case)
        # 3. The ENST ID so we know what gene this value is derived for
        cmd = "bcftools query -i \'GT=\"alt\"'  -f \'[%SAMPLE\\t%ID\\t%GT\\n]\' -o /test/" + file_prefix + "." + chromosome + ".parsed.txt /test/" + file_prefix + "." + chromosome + ".SAIGE.bcf"
        run_cmd(cmd, True)
        # This is just a list-format of the above file's header so we can read it in below with index-able columns
        header = ['sample', 'varID', 'genotype']

        samples = {}
        # Now we are going to read this file in and parse the genotype information into the dictionary
        # we created above (samples). This dictionary has a structure like:
        # {'sample1': {'gene1': 1}, 'sample3': {'gene2': 1}}
        # Note that only samples with a qualifying variant are listed
        filtered_variants = csv.DictReader(open(file_prefix + "." + chromosome + '.parsed.txt', 'r', newline='\n'),
                                           fieldnames=header,
                                           delimiter="\t",
                                           quoting=csv.QUOTE_NONE)
        for var in filtered_variants:
            # IF the gene is already present for this individual, increment based on genotype
            # ELSE the gene is NOT already present for this individual, instantiate a new
            # level in the samples dict and set according to current genotype
            ENST = snp_gene_map[var['varID']]
            if var['sample'] in samples:
                if ENST in samples[var['sample']]:
                    samples[var['sample']][ENST] += 1 if (var['genotype'] == '0/1') else 2
                else:
                    samples[var['sample']][ENST] = 1 if (var['genotype'] == '0/1') else 2
            else:
                samples[var['sample']] = {ENST: 1 if (var['genotype'] == '0/1') else 2}

        # We have to write this first into plink .ped format and then convert to bgen for input into BOLT
        # We are tricking BOLT here by setting the individual "variants" within bolt to genes. So our map file
        # will be a set of genes, and if an individual has a qualifying variant within that gene, setting it
        # to that value
        output_map = open(file_prefix + "." + chromosome + '.BOLT.map', 'w')
        output_ped = open(file_prefix + "." + chromosome + '.BOLT.ped', 'w')
        output_fam = open(file_prefix + "." + chromosome + '.BOLT.fam', 'w')

        # Make map file (just list of genes with the chromosome and start position of that gene):
        for gene in genes:
            output_map.write("%s %s 0 %i\n" % (genes[gene]['CHROM'],
                                               gene,
                                               genes[gene]['min_poss']))

        # Make ped / fam files:
        # ped files are coded with dummy genotypes of A A as a ref individual and A C as a carrier

        # We first need a list of samples we expect to be in this file. We can get this from the REGENIE .psam
        poss_indv = []  # This is just to help us make sure we have the right numbers later
        sample_file = csv.DictReader(open(file_prefix + "." + chromosome + '.sample', 'r', newline='\n'),
                                     delimiter=" ",
                                     quoting=csv.QUOTE_NONE)
        for sample in sample_file:
            if sample['ID'] != "0":  # This gets rid of the weird header row in bgen sample files...
                sample = sample['ID']
                poss_indv.append(sample)
                output_ped.write("%s %s 0 0 0 -9 " % (sample, sample))
                output_fam.write("%s %s 0 0 0 -9\n" % (sample, sample))
                genes_processed = 0
                for gene in genes:
                    genes_processed += 1  # This is a helper value to make sure we end rows on a carriage return (\n)
                    if sample in samples:
                        if gene in samples[sample]:
                            if genes_processed == len(genes):
                                output_ped.write("A C\n")
                            else:
                                output_ped.write("A C ")
                        else:
                            if genes_processed == len(genes):
                                output_ped.write("A A\n")
                            else:
                                output_ped.write("A A ")
                    else:
                        if genes_processed == len(genes):
                            output_ped.write("A A\n")
                        else:
                            output_ped.write("A A ")

        output_ped.close()
        output_map.close()
        output_fam.close()

        # And convert to bgen
        # Have to use OG plink to get into .bed format first
        cmd = 'plink --threads 1 --memory 9000 --make-bed --file /test/' + file_prefix + "." + chromosome + '.BOLT --out /test/' + file_prefix + "." + chromosome + '.BOLT'
        run_cmd(cmd, True)
        # And then use plink2 to make a bgen file
        cmd = "plink2 --threads 1 --memory 9000 --export bgen-1.2 'bits='8 --bfile /test/" + file_prefix + "." + chromosome + ".BOLT --out /test/" + file_prefix + "." + chromosome + ".BOLT"
        run_cmd(cmd, True)

        # Purge unecessary intermediate files to save space on the AWS instance:
        os.remove(file_prefix + "." + chromosome + '.BOLT.ped')
        os.remove(file_prefix + "." + chromosome + '.BOLT.map')
        os.remove(file_prefix + "." + chromosome + '.BOLT.fam')
        os.remove(file_prefix + "." + chromosome + '.BOLT.bed')
        os.remove(file_prefix + "." + chromosome + '.BOLT.bim')
        os.remove(file_prefix + "." + chromosome + '.BOLT.log')
        os.remove(file_prefix + "." + chromosome + '.BOLT.nosex')

        return poss_indv, samples

    # This gets information relating to included variants in bcf format files (per-variant)
    @staticmethod
    def _check_vcf_stats(poss_indv: list, genotypes: dict) -> pd.DataFrame:

        sample_table = pd.DataFrame(data={'sample_id': poss_indv})
        geno_dict = {'sample_id': [], 'ac': [], 'ac_gene': []}
        for sample in genotypes:
            samp_total = 0
            gene_total = 0
            for gene in genotypes[sample]:
                samp_total += genotypes[sample][gene]
                gene_total += 1
            geno_dict['sample_id'].append(sample)
            geno_dict['ac'].append(samp_total)
            geno_dict['ac_gene'].append(gene_total)

        found_genotypes = pd.DataFrame.from_dict(geno_dict)
        sample_table = pd.merge(sample_table, found_genotypes, on='sample_id', how="left")
        sample_table['ac'] = np.where(sample_table['ac'] >= 1, sample_table['ac'], 0)
        sample_table['ac_gene'] = np.where(sample_table['ac_gene'] >= 1, sample_table['ac_gene'], 0)

        return sample_table
