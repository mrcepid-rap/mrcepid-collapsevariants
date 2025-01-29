# mrcepid-collapsevariants Developer Readme

## Testing

Tests for this applet are implemented using the [mrcepid-testing](https://github.com/mrcepid-rap/mrcepid-testing) suite.
Please see this repository for more specifics on how testing works in the DNANexus environment.

Note for development: the test data to run unit tests for this applet are living on GCloud.
If you would like to gain access, please contact Eugene Gardner.
If you already have access to GCloud, run this command in the `/test_data/` directory:

```
gcloud storage rsync gs://iiuk-human-genetics/Syris/test_data/collapsevariants .
```

### System Requirements

The following language-specific packages must be installed on your system for tests to run:

* R
  * data.table
  * jsonlite
  * tidyverse
  * broom
* python
  * dxpy
* command-line tools
  * bcftools
  * plink
  * plink2
  * bgenix
  * docker

### Simulated Data

To enable tests for this applet, a large amount of data needs to be simulated. This is to enable testing of:

1. The collapsing analyses that this applet performs
2. Downstream associationtesting modules with the same data and without re-simulating different data in multiple locations.

The exact data that needs to be simulated here is:

1. A corresponding .bed + .fam + .map of common alleles that can be fed into runassociationtesting
2. WES BGEN files for chromosomes 1 : 22 & X, with:
   * SNPs with specific annotations that vary based on
     + MAF
     + Consequence
     + Some reasonable number of additional variant categories
   * SNPs within the same GENE that add up to some effect
   * Individual SNPs that add up to some effect
   * Collective GENEs that add up to some effect
3. A set of base covariates (PC1..PC40, age, age_squared, sex, wes_batch)
4. Simulated effects for all of the above

Simulation of data proceeds in two parts:

1. Simulation of a binary plink dataset that represents the 'genotyped' SNPs available within UKBiobank and DNANexus.
2. Simulation of a phenotype consisting of technical / phenotypic covariates and genetic covariates, of which genetic covariates include:
   * Genetic Markers
   * WES Data

**Note:** All of these simulations rely on fixing the random number generator for each step to ensure consistent results 
when these data need to be regenerated. This number is fixed to '1234' in the two scripts that use random generation 
outlined in the sections below. 

#### 1. Generating Binary Plink Files 

To simulate 'genotype' data, we use the [simuPOP](http://bopeng.github.io/simuPOP/index.html) tool developed for python3. 
Use of this tool requires a set of known variants to draw from to generate realistic haplotypes. This data can be 
acquired from the original [1000 Genomes project hg19 version](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/).
simuPOP requires a relatively long time to install (~30m) and we have created a docker image to run the simulations 
outlined below. This image is in a public [Docker Repo](https://hub.docker.com/repository/docker/egardner413/simupop/general), 
with a `Dockerfile` also available in this testing suite at `test/simuPop.Dockerfile`. Using this image should also ensure
reproducibility due to potential differences in random number generators used by different operating systems. A workflow to use 
simuPOP (and other tools) to simulate data required for running this test suite is outlined below:

For the more detail-oriented: The below commands use a mix of cloud-instances and local workstation. Generally, when 
downloading files, a cloud instance was used and data subsequently uploaded to the DNANexus platform. When running
the `simulate_chromosome.py` script, the 'swiss army knife' applet provided by DNANexus was used.

```commandline
# Download phased, autosomal hg19 chromosomal data:
perl -e 'for (1..22) {print "wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$_" . ".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz;\n";}' | bash
dx upload --destination project_resources/simulation_data/1KGP_phase3_hg19/ *.genotypes.vcf.gz

# Download sample manifest and get a list of GBR samples (note that I already checked, and none of the GBR samples are part of a trio):
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
perl -ane 'chomp $_; if ($F[1] eq 'GBR') {print "$F[0]\n";}' integrated_call_samples_v3.20130502.ALL.panel > GBR.txt

# Download a genetic map (using the same one for BOLT from the Alkes group) and format for sim1000G:
wget https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg19.txt.gz
gzcat genetic_map_hg19.txt.gz | perl -ane 'chomp $_; if ($F[0] eq "chr") {print "Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)\n";} else {$F[0] = "chr" . $F[0]; print join("\t", @F) . "\n";}' | gzip - > genetic_map_hg19.formatted.txt.gz 

# Filter the VCF to only GBR individuals and exclude multiallelics and SVs:
perl -e 'for (1..22) {print "nohup bcftools view -S GBR.txt -i '\''MULTI_ALLELIC=0 & VT != \"SV\"'\'' \
    -Ob ALL.chr$_.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | \
    bcftools view --min-af 0.1 --max-af 0.9 -Ob -o GBR.chr$_.vcf.gz - > $_.out &\n";}' | bash

# Upload the resulting GBR-only VCF data to DNANexus (I use a folder called 'project_resources/simulation_data/1KGP_phase3_hg19/', but it can be any named folder):
dx upload --brief --destination project_resources/simulation_data/1KGP_phase3_hg19/ GBR.chr*.vcf.gz

# Download rsID information and upload to DNANexus

# Run the chromosome simulation script for each chromosome:
# Can either run using 
# 1. swiss-army-knife on DNANexus (Make sure the script file is uploaded to DNANexus PRIOR to running this command).
# The perl one-liner wrapper just runs a seperate command for chr 1 -> 22 to parallelise tasks
perl -e 'for (1..22) {print "dx run swiss-army-knife --instance-type mem1_ssd1_v2_x36 -iin={'\''simulate_chromosome.py'\'','\''/project_resources/simulation_data/simulated_genotypes/genetic_map_hg19.formatted.txt.gz'\'','\''/project_resources/simulation_data/1KGP_phase3_hg19/GBR.chr$_.vcf.gz'\'','\''/Bulk/Genotype Results/Genotype calls/ukb22418_c$_" . "_b0_v2.bim'\''} -icmd='\''python3 simulate_chromosome.py $_ GBR.chr$_.vcf.gz ukb22418_c$_" . "_b0_v2.bim'\'' -iimage='\''egardner413/simupop'\'';\n";}' | bash
# 2. Locally using docker (example for chromosome 22):
docker run -v /path/to/current/directory/:/home/ -w /home/ egardner413/simupop bash -c 'python3 simulate_chromosome.py 22 GBR.chr22.vcf.gz ukb22418_c22_b0_v2.bim'

# Merge the resulting data, get low MAC sites, calculate a rough GRM, and calculate first 40 PCs from the resulting files using plink
ls sim_chromosome_*.bed | sed 's_.bed__g' > bgen_list.txt
plink2 --pmerge-list bgen_list.txt bfile --make-bed --out sim_chromosome
plink2 --bfile sim_chromosome --max-mac 100 --write-snplist --out sim_chromosome
plink2 --maf 0.01 --bfile sim_chromosome --make-rel --out sim_chromosome
plink --pca 40 'header' --bfile sim_chromosome --out sim_chromosome

dx upload --brief --destination project_resources/simulation_data/simulated_genotypes/ sim_chromosome.bed
dx upload --brief --destination project_resources/simulation_data/simulated_genotypes/ sim_chromosome.bim
dx upload --brief --destination project_resources/simulation_data/simulated_genotypes/ sim_chromosome.fam
dx upload --brief --destination project_resources/simulation_data/simulated_genotypes/ sim_chromosome.snplist
dx upload --brief --destination project_resources/simulation_data/simulated_genotypes/ sim_chromosome.eigenvec

# See below how to create these two files using python...
dx upload --brief --destination project_resources/simulation_data/simulated_genotypes/ sim_chromosome.sparseGRM.mtx
dx upload --brief --destination project_resources/simulation_data/simulated_genotypes/ sim_chromosome.sparseGRM.mtx.sampleIDs.txt
```

The following code will generate a relatedness matrix for the genetic data from the output of `plink2 make-rel`:

```python
from pathlib import Path

with Path('sim_chromosome.rel').open('r') as rel,\
    Path('sim_chromosome.sparseGRM.mtx').open('w') as mtx,\
        Path('sim_chromosome.sparseGRM.mtx.sampleIDs.txt').open('w') as mtx_sample:

    to_print = list()
    i = 1
    for line in rel:
        line = line.rstrip()
        data = map(float, line.split('\t'))

        for j, king in enumerate(data):
            j_mod = j + 1
            if king >= 0.02:
                if i == j_mod:
                    king = 0.5
                else:
                    king /= 2
                to_print.append({'i': i, 'j': j_mod, 'kin': king})

        i += 1

    mtx.write(f'%%MatrixMarket matrix coordinate real symmetric\n')
    mtx.write(f'10000 10000 {len(to_print)}\n')
    for rec in to_print:
        mtx.write(f'{rec["i"]} {rec["j"]} {rec["kin"]}\n')

    for samp in range(1000000, 1010000):
        mtx_sample.write(f'{samp}\n')
```

#### 2. Generating Phenotype and WES Data

Here we generate simulated WES data. While the generation of phenotypic data is not specifically necessary for running 
tests for this module, at this time we also generate genes and variants that have a simulated (random) effect on our
simulated phenotype. This ensures that we can use the data generated here to test associationtesting modules downstream.

```commandline
# Run the phenotype & WES simulation script: – this should take ~1hr on a 4 core machine. This script requires the 
# data.table, tidyverse, broom, and jsonlite packages for R to be installed.
./simulate_data.R

# Convert to .bgen & index
perl -e 'for (1..22,"X") {print "plink2 --import-dosage sim_data/chr" . $_ . ".filtered.traw skip0=1 skip1=2 noheader chr-col-num=1 pos-col-num=4 ref-first --fam sim_data/filtered.fam --export bgen-1.2 'bits='8 'sample-v2' --out sim_data/chr" . $_ . ".filtered; bgenix -g sim_data/chr" . $_ . ".filtered.bgen -index;\n"}' | bash

# Sort, bgzip, tabix annotations:
perl -e 'for (1..22,"X") {print "(head -n 1 sim_data/chr" . $_ . ".filtered.vep.tsv && tail -n +2 sim_data/chr" . $_ . ".filtered.vep.tsv | sort -k2,2n) | bgzip -c > sim_data/chr" . $_ . ".filtered.vep.tsv.gz; tabix -c C -s 1 -b 2 -e 2 sim_data/chr" . $_ . ".filtered.vep.tsv.gz\n";}' | bash

# Upload the data to DNANexus (I use a folder called 'project_resources/simulation_data/simulated_wes/', but it can be any named folder)
dx upload --brief --destination project_resources/simulation_data/simulated_wes/ sim_data/*.bgen*
dx upload --brief --destination project_resources/simulation_data/simulated_wes/ sim_data/*.sample
dx upload --brief --destination project_resources/simulation_data/simulated_wes/ sim_data/*.vep.tsv.gz*

#And generate a file list and upload to DNANexus (this script always generates a file named 'sim_data/bgen_locs.tsv'):
python3 build_bgen_list.py /project_resources/simulation_data/simulated_wes/
dx upload --destination project_resources/simulation_data/simulated_wes/ sim_data/bgen_locs.tsv
```

## Running Tests

To run tests, we use the mrcepid-testing suite. Test data required for tests should be located in the `/test_data/` 
folder included in this repository. If that folder is empty, then please follow the instructions to simulate data as 
outlined above. The DNANexus file ID (`file-GQz5Z78J0zVXK6633vP3GV8K`) is a place-holder. This will need to be replaced 
with the file-ID for the bgen location file uploaded in the above section.

An example command-line to run tests is provided below:

```commandline
# --script is the pytest compatible script
# --files are the test data required for testing
# --root_dir is the path to the root directory containing the source code for mrcepid-collapsevariants
# --modules are modules required for the current test. A branch (e.g., v1.1.0) of a given module can be requested using syntax like: general_utilities:v1.1.0 
# --add_opts are additional options to pass to the applet itself. bgen_index MUST be set to the actual bgen_index generated above
./test_launch.py --script /path/to/mrcepid-collapsevariants/test/collapsevariants_test.py \ 
   --files /path/to/mrcepid-collapsevariants/test/test_data/ \
   --root_dir /path/to/mrcepid-collapsevariants/ \
   --modules general_utilities mrcepid-collapsevariants \
   --add_opts bgen_index:file-GQz5Z78J0zVXK6633vP3GV8K
```