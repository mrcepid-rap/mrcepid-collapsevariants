import csv
import sys

import pandas as pd
import numpy as np

from pathlib import Path
from typing import Tuple, List, Dict

from collapsevariants.tool_parsers.saige_parser import GeneDict
from general_utilities.job_management.command_executor import CommandExecutor
from general_utilities.mrc_logger import MRCLogger

LOGGER = MRCLogger(__name__).get_logger()


def parse_filters_BOLT(file_prefix: str, chromosome: str, genes: Dict[str, GeneDict],
                       snp_gene_map: dict, cmd_exec: CommandExecutor) -> Tuple[List[str], Dict[str, Dict[str, int]]]:
    """Generate input format files that can be provided to BOLT

    BOLT only runs on individual variants and has no built-in functionality to perform burden testing. Thus, to run BOLT
    in an approximation of 'burden' mode we have to trick it by first collapsing all variants found within a gene across
    a set of individuals into a single 'variant' that contains information for all individuals for that gene. Thus we
    create what amounts to a binary yes/no found a variant in the listed gene for a given individual.

    Gene-Variants are encoded as heterozygous, with all genes having a REF allele of A and alternate allele of C. All
    variant IDs are ENST IDs for MANE transcripts (where possible). The actual output of this code is not returned,
    but is spilled to disk as a v1.2, ref-last, 8-bit .bgen format file with a partnered .sample file.

    :param file_prefix: A name to append to beginning of output files.
    :param chromosome: The chromosome currently being processed. This must be the short form of the chromosome name
        (e.g., '1' not 'chr1').
    :param genes: A Dictionary of genes identified by saige_parser.parse_filters_SAIGE of ENSTs mapped to variant IDs
    :param snp_gene_map: A Dictionary of variants identified by saige_parser.parse_filters_SAGE of variantIDs
        mapped to ENSTs
    :param cmd_exec: An instance of CommandExecutor to run commands on a docker container
    :return: A Tuple containing a List of samples that were found and a Dictionary with keys of sample IDs and values of
        a dictionary of ENST / Genotype pairs for that given individual
    """

    # Get out BCF file into a tab-delimited format that we can parse for BOLT.
    # We ONLY want alternate alleles here (-i 'GT="alt") and then for each row print:
    # 1. Sample ID: UKBB eid format
    # 2. The actual genotype (0/1 or 1/1 in this case)
    # 3. The ENST ID so we know what gene this value is derived for
    cmd = f'bcftools query -i \'GT="alt"\' -f \'[%SAMPLE\\t%ID\\t%GT\\n]\' ' \
          f'-o /test/{file_prefix}.{chromosome}.parsed.txt ' \
          f'/test/{file_prefix}.{chromosome}.SAIGE.bcf'
    cmd_exec.run_cmd_on_docker(cmd)
    # This is just a list-format of the above file's header so we can read it in below with index-able columns
    header = ['sample', 'varID', 'genotype']

    samples: Dict[str, Dict[str, int]] = dict()
    # Now we are going to read this file in and parse the genotype information into the dictionary
    # we created above (samples). This dictionary has a structure like:
    # {'sample1': {'gene1': 1}, 'sample3': {'gene2': 1}}
    # Note that only samples with a qualifying variant are listed
    with Path(f'{file_prefix}.{chromosome}.parsed.txt').open('r') as filtered_reader:
        filtered_variants = csv.DictReader(filtered_reader,
                                           fieldnames=header,
                                           delimiter='\t',
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

    # We have to write this first into .vcf format and then convert to .bgen for input into BOLT
    # We are tricking BOLT here by setting the individual "variants" within bolt to genes. So our .vcf file
    # will be a set of genes, and if an individual has a qualifying variant within that gene, setting genotype
    # to 0/1

    vcf_path = Path(f'{file_prefix}.{chromosome}.BOLT.vcf')

    with vcf_path.open('w') as output_vcf, \
            Path(f'{file_prefix}.{chromosome}.sample').open('r') as sample_reader:

        # First extract samples from the sample file
        sample_file = csv.DictReader(sample_reader, delimiter=' ', quoting=csv.QUOTE_NONE)
        poss_indv = []  # This is just to help us make sure we have the right numbers later
        for sample in sample_file:
            if sample['ID'] != "0":  # This gets rid of the weird header row in bgen sample files...
                sample = sample['ID']
                poss_indv.append(sample)

        # Write the vcf header
        output_vcf.write('##fileformat=VCFv4.2\n')
        output_vcf.write(f'##contig=<ID={chromosome},length=999999999>\n') # fake len as it doesn't matter to .bgen
        output_vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        sample_string = '\t'.join(poss_indv)
        output_vcf.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_string}\n')

        # Write the vcf body
        for gene in genes:
            gt_array = []
            for sample in poss_indv:
                if sample in samples:
                    if gene in samples[sample]:
                        gt_array.append('0/1')
                    else:
                        gt_array.append('0/0')
                else:
                    gt_array.append('0/0')
            gt_string = '\t'.join(gt_array)
            output_vcf.write(f'{genes[gene]["CHROM"]}\t{genes[gene]["min_poss"]}\t{gene}\tA\tC\t.\tPASS\t.\tGT\t{gt_string}\n')

    # And convert to .bgen
    cmd = f'qctool -filetype "vcf" -ofiletype "bgen_v1.2" -sort ' \
          f'-g /test/{vcf_path} ' \
          f'-og /test/{file_prefix}.{chromosome}.BOLT.bgen ' \
          f'-os /test/{file_prefix}.{chromosome}.BOLT.sample'
    cmd_exec.run_cmd_on_docker(cmd)

    vcf_path.unlink()

    return poss_indv, samples


def check_vcf_stats(poss_indv: List[str], genotypes: Dict[str, Dict[str, int]]) -> pd.DataFrame:
    """Get information relating to included variants in bcf format files (per-variant)

    This method calculates per-sample and per-gene totals for this chromosome

    :param poss_indv: Identical to return 1 from parse_filters_BOLT – a List of samples that were found
    :param genotypes: Identical to return 2 from parse_filters_BOLT – a Dictionary with keys of sample IDs and values of
        a dictionary of ENST / Genotype pairs for that given individual
    :return: A pandas.DataFrame containing per-sample and per-gene totals for this chromosome
    """

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
