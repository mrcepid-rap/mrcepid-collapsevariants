#!/usr/bin/env python3
import os
import csv
import sys
import random
import shutil
import pandas as pd

from pathlib import Path
from typing import Tuple
from pysam import VariantFile

# simuPOP has to be imported this way to enable optionsetting
import simuOpt

from general_utilities.job_management.command_executor import CommandExecutor

simuOpt.setOptions(numThreads=os.cpu_count())
import simuPOP as sim

# Set RNGs to ensure consistent results:
random.seed(1234)
sim.getRNG().set('mt19937', 1234)

CMD_EXEC = CommandExecutor()


def get_recomb_rate(position: int) -> Tuple[float, float]:
    """ Get the cM position and recombination rate for a given locus.

    This calculation was adapted from simuGWAS by Bo Peng with modifications to use an IntervalTree rather than an
    ordered list search.

    :param position: bp position of a variant on a given chromosome
    :return: A tuple of cM and recombination rate (in cM / Mb)
    """

    found_intervals = gm.query('`Position(bp)` == @position')
    if position < start_interval['Position(bp)']:
        cM = start_interval['Map(cM)'] - (start_interval['Position(bp)'] - position) * 1e-8
        cM = 0 if cM < 0 else cM  # cM doesn't make sense less than 0.
        rate = start_interval['Rate(cM/Mb)']

    elif position > end_interval['Position(bp)']:
        cM = end_interval['Map(cM)'] + (position - end_interval['Position(bp)']) * 1e-8
        rate = end_interval['Rate(cM/Mb)']

    elif len(found_intervals) == 1:
        found_interval = gm.iloc[found_intervals.index[0]]
        cM = found_interval['Map(cM)']
        rate = found_interval['Rate(cM/Mb)']

    else:
        found_intervals = gm.query('`Position(bp)` > @position')
        next_interval = gm.iloc[found_intervals.index[0]]
        cM = next_interval['Map(cM)'] - (next_interval['Position(bp)'] - position) * 1e-8
        cM = 0 if cM < 0 else cM
        rate = next_interval['Rate(cM/Mb)']

    return cM, rate


def exp_demo_func(N0, N1, expandGen):
    '''
    Return an exponential population expansion demographic function.

    This function is adapted from simuGWAS by Bo Peng.

    N0: a list of initial subpopulation sizes.
    N1: final population size.

    gen: generations to evolve.
    '''
    N1 = [int(N1 * 1.0 * x / sum(N0)) for x in N0]
    step = [float(x - y) / expandGen for x, y in zip(N1, N0)]

    def func(gen):
        # try to make the exact population size
        if gen == expandGen - 1:
            return N1
        else:
            return [int(x + (gen + 1) * y) for x, y in zip(N0, step)]

    return func


chromosome = sys.argv[1]
vcf_file = Path(sys.argv[2])
bim_path = Path(sys.argv[3])
recomb_map = Path('genetic_map_hg19.formatted.txt.gz')
print(f'Starting simulation of chromosome {chromosome}.')

# Iterate through the appropriate .bim file from the real data, so we can have approx. the same number of variants:
with bim_path.open('r') as bim_reader:
    n_vars = 0
    for _ in bim_reader:
        n_vars += 1
print(f'Number of variants to be simulated: {n_vars}')

# Make sure the genetic Map has only the correct chromosomes on it:
gm = pd.read_csv(recomb_map, sep='\t', dtype={'Chromosome': str})
gm = gm[gm['Chromosome'] == chromosome]
gm.reset_index(inplace=True)

# Get rsID values (required by BOLT). bigBedToBed is already available within the Docker this script runs on:
get_snps_cmd = f'bigBedToBed https://hgdownload.soe.ucsc.edu/gbdb/hg19/snp/dbSnp153Common.bb ' \
               f'-chrom=chr{chromosome} -header rsID.bed'
CMD_EXEC.run_cmd(get_snps_cmd)
dbSNP = pd.read_csv('rsID.bed', sep='\t')
dbSNP = dbSNP[['#chrom', 'chromEnd', 'name']]  # Only chromEnd required since table is in 0-based bed format

# Start first and last intervals for easy calculations:
start_interval = gm.iloc[0]
end_interval = gm.iloc[len(gm) - 1]

# Decide a random subset of variants to include based on total number of variants:
# Also have to remove multiallelic sites
with VariantFile(vcf_file) as loaded_vcf_file:
    num_multiallelics = 0
    last_pos = 0
    variant_ranges = []
    for var_no, variant in enumerate(loaded_vcf_file.fetch()):
        if variant.pos != last_pos:
            variant_ranges.append(var_no)
        else:
            num_multiallelics += 1
        last_pos = variant.pos

    kept_loci_index = random.sample(variant_ranges, n_vars)
    kept_loci_index.sort()
print(f'Total variants found in parent VCF             : {var_no}')
print(f'Total multiallelic variants found in parent VCF: {num_multiallelics}')

# Import the VCF file:
print(f'Loading {n_vars} variants into simuPOP.')
with VariantFile(vcf_file) as loaded_vcf_file:

    num_indvs = len(loaded_vcf_file.header.samples)
    pop = sim.Population(size=num_indvs, ploidy=2, chromNames=[chromosome])

    variants = []

    total_variant_no = 0  # This tracks the total number of variants in pop, NOT in the VCF
    for variant_no, variant in enumerate(loaded_vcf_file.fetch()):
        if variant_no in kept_loci_index:
            ref_allele = variant.ref
            alt_allele = variant.alts[0]  # I have filtered out all multi-allelics
            # Decide varID based on presence / absence in dbSNP database:
            foundSNPs = dbSNP.query('chromEnd == @variant.pos')
            if len(foundSNPs) == 1:
                varID = foundSNPs['name'].item()
            else:
                varID = ':'.join(map(str, [variant.chrom, variant.pos, ref_allele, alt_allele]))
            cM, rate = get_recomb_rate(variant.pos)
            variants.append({'CHR': chromosome, 'SNP': varID, 'POS': variant.pos, 'COUNTED': ref_allele,
                             'ALT': alt_allele, 'CM': cM, 'rate': rate})
            pop.addLoci(chrom=pop.chromByName(chromosome), pos=variant.pos, lociNames=[varID],
                        alleleNames=[ref_allele, alt_allele])

            for sample_no, genotype in enumerate(variant.samples.values()):
                for allele_no, allele in enumerate(genotype.allele_indices):
                    pop.individual(sample_no).setAllele(allele, total_variant_no, allele_no)
            total_variant_no += 1

    # And get all rates into an ordered list to create a recombinator
    rates = [var['rate'] for var in variants]

print(f'Initiating simuPOP evolution simulation to generate 10,000 individuals.')
demo_func = exp_demo_func(pop.subPopSizes(), 10000, 500)
simulator = sim.Simulator(pop, rep=1)
rec_op = sim.Recombinator(rates=rates, loci=[x for x in range(0, n_vars)])
simulator.evolve(
    initOps=[
        sim.InitSex()
    ],
    preOps=[
        sim.SNPMutator(u=0, v=0, loci=[x for x in range(0, n_vars - 1)]),
        sim.NoneOp()
    ],
    matingScheme=sim.ControlledRandomMating(ops=rec_op, subPopSize=demo_func),
    gen=500
)

pop = simulator.extract(0)

# Now parse the data and generate a plink-compatible ped/bim/fam
print(f'Simulation complete, writing .traw format file (sim_chromosome_{chromosome}.traw) for input to plink2.')
with Path(f'sim_chromosome_{chromosome}.traw').open('w') as traw_writer,\
        Path(f'sim_chromosome_{chromosome}.fam').open('w') as fam_writer:

    sample_ids = [f'{id + 1000000}' for id in range(10000)]
    for id in sample_ids:
        fam_writer.write(f'{id} {id} 0 0 0 -9\n')

    sample_ids = [f'{id}_{id}' for id in sample_ids]
    traw_csv = csv.DictWriter(traw_writer, delimiter='\t',
                              fieldnames=['CHR', 'SNP', 'CM', 'POS', 'COUNTED', 'ALT'] + sample_ids,
                              extrasaction='ignore')
    traw_csv.writeheader()

    for variant_no in range(n_vars):
        variant = variants[variant_no].copy()
        for sample_no in range(10000):
            # I think .traw format is a bit weird in that the 'count' is of the REF allele, so 2 - alt. Regardless, it
            # imports into plink in the way you would expect
            variant[sample_ids[sample_no]] = \
                2 - sum([pop.individual(sample_no).allele(variant_no, ploidy) for ploidy in range(2)])

        traw_csv.writerow(variant)

print(f'Converting .traw to plink .bed format and writing final output.')
plink_cmd = f'plink2 --import-dosage sim_chromosome_{chromosome}.traw skip0=1 skip1=2 id-delim=_ chr-col-num=1 ' \
            f'pos-col-num=4 ref-first ' \
            f'--fam sim_chromosome_{chromosome}.fam ' \
            f'--make-bed ' \
            f'--out sim_chromosome_{chromosome}'
CMD_EXEC.run_cmd(plink_cmd, print_cmd=True)

# Delete files we don't want uploaded
Path(f'sim_chromosome_{chromosome}.traw').unlink()
Path(f'sim_chromosome_{chromosome}.log').unlink()
Path('rsID.bed').unlink()
if Path('udcCache/').exists():
    shutil.rmtree('udcCache/')
