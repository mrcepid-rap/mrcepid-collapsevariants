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

This applet is step 3 (mrc-collapsevariants) of the rare variant testing pipeline developed by Eugene Gardner for the UKBiobank
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
|input_vcf  | Input vcf file from [mrcepid-annotatecadd](https://github.com/mrcepid-rap/mrcepid-filterbcf) to annotate with CADD |
|filtering_expression | [bcftools](https://samtools.github.io/bcftools/bcftools.html) compatible filtering expression. See [above](#1-filtering-with-bcftools-filtering-expressions) |
|file_prefix | descriptive file prefix for output name |

**BIG Note:** The value provided to `file_prefix` **MUST** be identical for all VCF files that you wish to merge and test during
rare variant burden testing.

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
dx run mrcepid-collapsevariants --priority low --destination filtered_vcfs/ -iinput_vcf=file-A12345 \
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

t.b.d. A fast way of performing filtering for a large number of VCF files
