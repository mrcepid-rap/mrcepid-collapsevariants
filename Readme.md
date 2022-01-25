# CollapseVariants (DNAnexus Platform App)

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

### Table of Contents

- [Introduction](#introduction)
    * [Background](#background)
    * [Dependencies](#dependencies)
        + [Docker](#docker)
        + [Resource Files](#resource-files)
- [Methodology](#methodology)
    * [1. Filtering With BCFTools Filtering Expressions](#1-filtering-with-bcftools-filtering-expressions)
    * [2. Generating Outputs](#2-generating-outputs)
- [Running on DNANexus](#running-on-dnanexus)
    * [Inputs](#inputs)
    * [Outputs](#outputs)
    * [Command line example](#command-line-example)
        + [Batch Running](#batch-running)

## Introduction

This applet generates raw data necessary to perform rare variant burden testing using [bcftools](https://samtools.github.io/bcftools/bcftools.html)
or [plink2](https://www.cog-genomics.org/plink/2.0/). Please see these two tool's respective documentation for more
information on how individual commands used in this applet work.

This README makes use of DNANexus file and project naming conventions. Where applicable, an object available on the DNANexus
platform has a hash ID like:

* file – `file-1234567890ABCDEFGHIJKLMN`
* project – `project-1234567890ABCDEFGHIJKLMN`

Information about files and projects can be queried using the `dx describe` tool native to the DNANexus SDK:

```commandline
dx describe file-1234567890ABCDEFGHIJKLMN
```

**Note:** This README pertains to data included as part of the DNANexus project "MRC - Variant Filtering" (project-G2XK5zjJXk83yZ598Z7BpGPk)

### Background

Downstream of this applet, we have implemented four tools / methods for rare variant burden testing:

* [BOLT](https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html)
* [SAIGE-GENE](https://github.com/weizhouUMICH/SAIGE/wiki/Genetic-association-tests-using-SAIGE)
* [STAAR](https://github.com/xihaoli/STAAR)
* GLMs – vanilla linear/logistic models implemented with python's [statsmodels module](https://www.statsmodels.org/stable/index.html)

These four tools / methods require very different input files to run. The purpose of this applet is to generate inputs 
that are compatible with each of these tools input requirements. This tool is part (1) of a two-step process (in bold):

1. **Generate initial files from each VCF filtered/annotated by [mrcepid-filterbcf](https://github.com/mrcepid-rap/mrcepid-filterbcf) 
   and [mrcepid-annotatecadd](https://github.com/mrcepid-rap/mrcepid-annotatecadd)**
2. Merge these resulting files into a single set of inputs for the four tools that we have implemented 

For more information on the format of these files, please see the [mrcepid-runassociationtesting](https://github.com/mrcepid-rap/mrcepid-runassociationtesting)
documentation.

### Dependencies

#### Docker

This applet uses [Docker](https://www.docker.com/) to supply dependencies to the underlying AWS instance
launched by DNANexus. The Dockerfile used to build dependencies is available as part of the MRCEpid organisation at:

https://github.com/mrcepid-rap/dockerimages/blob/main/filterbcf.Dockerfile

This Docker image is built off of the primary 20.04 Ubuntu distribution available via [dockerhub](https://hub.docker.com/layers/ubuntu/library/ubuntu/20.04/images/sha256-644e9b64bee38964c4d39b8f9f241b894c00d71a932b5a20e1e8ee8e06ca0fbd?context=explore).
This image is very light-weight and only provides basic OS installation. Other basic software (e.g. wget, make, and gcc) need
to be installed manually. For more details on how to build a Docker image for use on the UKBiobank RAP, please see:

https://github.com/mrcepid-rap#docker-images

In brief, the primary **bioinformatics software** dependencies required by this Applet (and provided in the associated Docker image)
are:

* [htslib and samtools](http://www.htslib.org/)
* [bcftools](https://samtools.github.io/bcftools/bcftools.html)
* [plink2](https://www.cog-genomics.org/plink/2.0/)

This list is not exhaustive and does not include dependencies of dependencies and software needed
to acquire other resources (e.g. wget). See the referenced Dockerfile for more information.

#### Resource Files

This applet does not have any external dependencies.

## Methodology

This applet is step 4 (mrc-collapsevariants) of the rare variant testing pipeline developed by Eugene Gardner for the UKBiobank
RAP at the MRC Epidemiology Unit:

![](https://github.com/mrcepid-rap/.github/blob/main/images/RAPPipeline.png)

This applet has two major steps:

1. Select variants from a filtered and annotated VCF file using a BCF tools filtering expression.
2. Generate various output files in varying formats that can be merged later as part of 
   [mrcepid-mergecollapsevariants](https://github.com/mrcepid-rap/mrcepid-mergecollapsevariants).

For more details please see the commented source code available at `src/mrcepid-collapsevariants.py` of this repository.

### 1. Filtering With BCFTools Filtering Expressions

The user of this applet must provide a filtering expression that is compatible with bcftools filtering expressions. For 
extensive details and tutorials on how to construct such expressions, please see the [bcftools EXPRESSIONS documentation](https://samtools.github.io/bcftools/bcftools.html#expressions).
Briefly, one can construct various filtering expressions to generate variants that they want to test during rare variant burden tests.
These expressions **MUST** be based on INFO fields generated by the annotation and filtering parts of this workflow. For possible
fields that can be filtered on, please see:

https://github.com/mrcepid-rap/mrcepid-annotatecadd#outputs

For possible consequences (PARSED_CSQ) to filter on, please see:

https://github.com/mrcepid-rap/mrcepid-filterbcf#4-parsing-vep-consequences

**Note:** One can also filter based on raw [VEP consequence](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences)
using the "CSQ" field if desired (e.g. stop_gained).

For example, lets say that I want to test the effects of rare (MAF < 1x10<sup>-3</sup>) Protein Truncating Variants (PTVs) 
on some phenotype. Thus, I could construct a filtering expression like:

`AF<0.001 & PARSED_CSQ="PTV"`

Note that this expression would retain variants that FAILED quality control and did not pass basic PTV QC filters (e.g. 
[LoFTEE](https://github.com/konradjk/loftee)). So, we could further modify this expression like so:

`FILTER="PASS" & AF<0.001 & LOFTEE="HC" & PARSED_CSQ="PTV"`

* FILTER="PASS" means we only retain variants that passed quality control
* LOFTEE="HC" retains only PTVs that are LoFTEE high-confidence LoFs

This filtering expression can be increased in complexity to generate a more stringent set of variants:

`FILTER="PASS" & AF<0.001 & LOFTEE="HC" & PARSED_CSQ="PTV" & CADD>25 & gnomAD_AF < 0.001`

And so forth... 

These expressions can be mixed and modified according to user need. In summary, the above expression is run by the applet
as part of the a command like:

```commandline
bcftools view -i 'FILTER="PASS" & AF<0.001 & LOFTEE="HC" & PARSED_CSQ="PTV" & CADD>25 & gnomAD_AF < 0.001' \ 
        -Ob -o variants.filtered.bcf variants.vcf.gz
```

The file `variants.filtered.bcf` is used for the next step.

### 2. Generating Outputs

This applet then performs a series of formatting steps to generate various output files that are compatible with the
different tools listed in [Background](#background). I am not going to go into detail on the format of these files as 
they are mainly intermediate and not used by any other applets / tools.

## Running on DNANexus

### Inputs

|input|description             |
|---- |------------------------|
|input_vcfs  | Input vcf file(s) from [mrcepid-annotatecadd](https://github.com/mrcepid-rap/mrcepid-annotatecadd) to filter for variants according to `filtering_expression` |
|filtering_expression | [bcftools](https://samtools.github.io/bcftools/bcftools.html) compatible filtering expression. See [above](#1-filtering-with-bcftools-filtering-expressions) |
|file_prefix | descriptive file prefix for output name |

**BIG Note:** The value provided to `file_prefix` **MUST** be identical for all VCF files that you wish to merge and test during
rare variant burden testing.

`input_vcfs` is a file list that **MUST** contain DNANexus file hash keys (e.g. like file-1234567890ABCDEFGHIJ). A simple
way to generate such a list is with the following bash/perl one-liner:

```commandline
dx ls -l filtered_vcfs/ukb23148_c7_b*_v1_chunk*cadd.bcf | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "$1\n";}' > collapse_list.txt
```

This command will:

1. Find all filtered vcf files on chromosome 7 and print in dna nexus "long" format which includes a column for file hash (column 7)
2. Extract the file hash using a perl one-liner and print one file hash per line

The final input file will look something like:

```text
file-1234567890ABCDEFGHIJ
file-2345678901ABCDEFGHIJ
file-3456789012ABCDEFGHIJ
file-4567890123ABCDEFGHIJ
```

This file then needs to be uploaded to the DNANexus platform, so it can be provided as input:

```commandline
dx upload collapse_list.txt
```

#### File Coordinate Reference

We also provide a list of all bcf chunks in the file: `coordinates.files.tsv.gz (file-G7YGYPjJJv8kz6QvP41q5KYg)` in `project-G6BJF50JJv8p4PjGB9yy7YQ2`. 
This file has the following rows:

1. chrom – chunk chromosome
2. start - coordinate of the first variant in this chunk
3. end - coordinate of the last variant in this chunk
4. chunk_prefix – chunk prefix for all files for this chunk (regex `/ukb23148_c[0-9XY]{1,2}_b\d+_v1_chunk\d/`)
5. bcf_dxpy – DNA Nexus file hash for the bcf file
6. bcf_indx_dxpy – DNA Nexus file hash for the bcf index (.csi) file
7. vep_dxpy – DNA Nexus file hash for the variant VEP annotation file
8. vep_index_dxpy – DNA Nexus file hash for the variant VEP annotation index (.tbi) file

This file is tab-delimited and tabix-indexed (.tbi) to enable searching for desired file-chunk(s). One can query this 
coordinates file to generate a list as described above for specific genetic regions.

### Outputs

|output                 | description       |
|-----------------------|-------------------|
|output_tarball         |  Output tarball containing filtered and processed variant counts  |

output_tarball is named based on the name of the `input_vcf` combined with `file_prefix` like:

`ukb23156_c1_b0_v1.norm.filtered.tagged.missingness_filtered.annotated.cadd.PTV.tar.gz`

While I am not going into detail about the format of the files contained in this tar file, I list here the files for 
record-keeping purposes. All files have a standard prefix identical to that of the tarball with an extra descriptor:

* <prefix>.BOLT.json – BOLT-ready .json file of sample - gene pairs.
* <prefix>.REGENIE.annotation.txt – Per variant annotation information in tsv format with coordinate and ID
* <prefix>.REGENIE.pgen – plink pgen format-file of filtered genotypes
* <prefix>.REGENIE.psam - plink psam format-file of filtered genotypes
* <prefix>.REGENIE.pvar - plink pvar format-file of filtered genotypes
* <prefix>.REGENIE.log - log file from creating pgen files
* <prefix>.STAAR.matrix.txt - STAAR-ready tsv file of sample - variant - genotype sets

### Command line example

If this is your first time running this applet within a project other than "MRC - Variant Filtering", please see our
organisational documentation on how to download and build this app on the DNANexus Research Access Platform:

https://github.com/mrcepid-rap

Running this command is fairly straightforward using the DNANexus SDK toolkit. For the input vcf (provided with the flag
`-iinput_vcf`) one can use a file hash from the VCF output of `mrcepid-annotatecadd`:

```commandline
dx run mrcepid-collapsevariants --priority low --destination filtered_vcfs/ -iinput_vcfs=file-A12345 \
        -ifiltering_expression='FILTER="PASS" & AF<0.001 & LOFTEE="HC" & PARSED_CSQ="PTV"' \
        -ifile_prefix="PTV"
```

Brief I/O information can also be retrieved on the command line:

```commandline
dx run mrcepid-collapsevariants --help
```

I have set a sensible (and tested) default for compute resources on DNANexus that is baked into the json used for building 
the app (at `dxapp.json`) so setting an instance type is unnecessary. This current default is for a mem1_ssd2_v2_x2 instance
(2 CPUs, 4 Gb RAM, 50Gb storage). If necessary to adjust compute resources, one can provide a flag like `--instance-type mem1_ssd1_v2_x4`.

#### Batch Running

It is easier to implement batch running manually, rather than use built-in DNANexus batch functionality. In brief, first
generate a list of all files that need to be run through the process as outlined [above](#inputs):

```commandline
dx ls -l filtered_vcfs/ukb23148_c7_b*_v1_chunk*cadd.bcf | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "$1\n";}' > collapse_list.txt
```

Then, simply use the *NIX default split command to generate a set of individual files that can work through all the files
found above on individual instances:

```commandline
split -a 1 -l 35 collapse_list.txt collapse_list_
```

A few important notes on the above:
1. We set the number of files per-list to 35 (using `-l 35`) because our default instance has 35 cores and requires 2 
   cores per file. This means we should be able to run a total of 35 files at a time, but we need a core to be able to 
   monitor these processes, thus why we do 35 files.
2. We CAN set the number of files to greater than 35, but this means other files need to finish processing before others
   can start, meaning runtime will be longer than expected.
3. This will create files named bcf_list_a, bcf_list_b, bcf_list_c, etc.

Then we upload to dna nexus, and generate a set of commands that will then run this applet:

```commandline
dx upload collapse_list_* --destination batch_lists/
dx ls -l batch_lists/collapse_list_* | \ 
    perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {print "dx run mrcepid-collapsevariants --priority low --yes --brief --destination filtered_vcfs/ -iinput_vcfs=$1 -ifiltering_expression='\''FILTER=\"PASS\" \& AF<0.0001 \& LOFTEE=\"HC\" \& PARSED_CSQ=\"PTV\"'\'' -ifile_prefix=\"PTV\";\n";}' | bash
```
