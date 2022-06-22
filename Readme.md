# CollapseVariants (DNAnexus Platform App)

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

### Table of Contents

- [Introduction](#introduction)
    * [Background](#background)
    * [Changelog](#changelog)
    * [Dependencies](#dependencies)
        + [Docker](#docker)
        + [Resource Files](#resource-files)
- [Methodology](#methodology)
    * [1. Filtering With Pandas Query Expressions](#1-filtering-with-pandas-query-expressions)
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

### Changelog

* v1.0.0
  * Initial numbered release. Changes going forward will be tracked in this section of the documentation

### Background

Downstream of this applet, we have implemented four tools / methods for rare variant burden testing:

* [BOLT](https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html)
* [SAIGE-GENE](https://github.com/weizhouUMICH/SAIGE/wiki/Genetic-association-tests-using-SAIGE)
* [STAAR](https://github.com/xihaoli/STAAR)
* GLMs – vanilla linear/logistic models implemented with python's [statsmodels module](https://www.statsmodels.org/stable/index.html)

These four tools / methods require very different input files to run. The purpose of this applet is to generate inputs 
that are compatible with each of these tools input requirements. For more information on the format of these input files, 
please see the [mrcepid-runassociationtesting](https://github.com/mrcepid-rap/mrcepid-runassociationtesting) documentation.

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
* [qctool](https://www.well.ox.ac.uk/~gav/qctool_v2/index.html)
* [bgenix](https://enkre.net/cgi-bin/code/bgen/doc/trunk/doc/wiki/bgenix.md)

This list is not exhaustive and does not include dependencies of dependencies and software needed
to acquire other resources (e.g. wget). See the referenced Dockerfile for more information.

#### Resource Files

This applet does not have any external dependencies.

## Methodology

This applet is step 5 (mrc-collapsevariants) of the rare variant testing pipeline developed by Eugene Gardner for the UKBiobank
RAP at the MRC Epidemiology Unit:

![](https://github.com/mrcepid-rap/.github/blob/main/images/RAPPipeline.png)

This applet has two major steps:

1. Select variants from a filtered and annotated VCF file using either:
   1. a [pandas query](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html) filtering expression 
   2. A list of variant IDs
   3. A list of HGNC gene symbols (one per line) combined with a filtering expression (see i.)
2. Generate various output files in varying formats per-chromosome to be fed into
   [mrcepid-runassociationtesting](https://github.com/mrcepid-rap/mrcepid-runassociationtesting).

For more details please see the commented source code available at `src/mrcepid-collapsevariants.py` of this repository.

### 1. Filtering Variants

#### i. Filtering With Pandas Query Expressions

The user of this applet can provide a filtering expression that is compatible with the pandas query function to select 
variants for association testing. For details and tutorials on how to construct such queries, please see the 
[pandas query documentation](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html). Briefly, one can
construct various expressions to generate variants that they want to test during rare variant burden tests.
These expressions **MUST** be based on fields within the annotation's tsv file provided as output of the [mrcepid-makebgen](https://github.com/mrcepid-rap/mrcepid-makebgen)
step of this pipeline. This file is stored on the RAP as file `file-G857Z4QJJv8x7GXfJ3y5v1qV` in project `project-G6BJF50JJv8p4PjGB9yy7YQ2`.
For possible fields and values that can be used for filtering, please see:

https://github.com/mrcepid-rap/mrcepid-annotatecadd#outputs

For possible consequences (PARSED_CSQ) to filter on, please see:

https://github.com/mrcepid-rap/mrcepid-filterbcf#4-parsing-vep-consequences

**Note:** One can also filter based on raw [VEP consequence](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences)
using the "CSQ" field, if desired (e.g. stop_gained).

For example, lets say that I want to test the effects of rare (MAF < 1x10<sup>-3</sup>) Protein Truncating Variants (PTVs) 
on some phenotype. Thus, I could construct a filtering expression like:

`AF<=0.001 & PARSED_CSQ=="PTV"`

Note that this expression would retain variants that FAILED quality control and did not pass basic PTV QC filters (e.g. 
[LoFTEE](https://github.com/konradjk/loftee)). So, we could further modify this expression like so:

`FILTER=="PASS" & AF<=0.001 & LOFTEE=="HC" & PARSED_CSQ=="PTV"`

* FILTER=="PASS" means we only retain variants that passed quality control
* LOFTEE=="HC" retains only PTVs that are LoFTEE high-confidence LoFs

This filtering expression can be increased in complexity to generate a more stringent set of variants:

`FILTER=="PASS" & AF<=0.001 & LOFTEE=="HC" & PARSED_CSQ=="PTV" & CADD>25 & gnomAD_AF < 0.001`

And so forth... 

These expressions can be mixed and modified according to user need. In summary, the above expression is run by the applet
using the pseudo-code:

```python
import pandas as pd
import gzip

variant_index = pd.read_csv(gzip.open('450k_vep.sorted.tsv.gz', 'rt'), sep = "\t")
variant_index = variant_index.query('FILTER=="PASS" & AF<=0.001 & LOFTEE=="HC" & PARSED_CSQ=="PTV" & CADD>25 & gnomAD_AF < 0.001')
```

The file `variants.filtered.bcf` is used for the next step.

#### ii. Filtering with Variant IDs

The user can also opt to provide a list of variant IDs that they want to use to build a single mask using the `snplist` 
input. This mode will create a single 'Gene' with the ID ENST00000000000 that includes information for all variants. Variants
listed can be from multiple chromosomes or a single gene. variant IDs **MUST** match the format found in the VEP annotation
files generated by [mrcepid-makebgen](https://github.com/mrcepid-rap/mrcepid-makebgen) applet:

```text
1:1234:A:G
2:4567:ATA:A
```

Any number of variants can be included in this mask but we have not tested the limits in terms of compute on generating 
a mask with a large number of variants. This list can also include variants NOT found in the VEP files (e.g. from a CRISPR
screen). This applet will report the variants NOT found in the log file. See [outputs](#outputs) for more information on 
this file.

Use of a SNP list will generate similar files to those generated when using a pandas-compatible filtering expression. The 
only major difference is that instead of a per-chromosome set of files, a single set of files with the prefix 'SNP' will 
be generated. This is to allow [mrcepid-runassociationtesting](https://github.com/mrcepid-rap/mrcepid-runassociationtesting)
to recognise that we are running collapsed SNPs rather than per-GENE tests.

**Big Note** – Masks generated with a SNP-list are currently only compatible with the phewas and extract modes of 
[mrcepid-runassociationtesting](https://github.com/mrcepid-rap/mrcepid-runassociationtesting).

#### iii. Filtering with HGNC gene symbols + query expression

Another option is to provide a list HGNC gene symbols **combined** with a pandas query expression to create a custom mask. The relevant parameter to specify the gene list is `genelist`. As with the SNP list input we create a single 'GENE' which here has the ID ENST99999999999. Many considerations for the SNP list input also apply to gene list input (e.g. maximum number of underlying variants etc.).
Gene IDs are checked against the HGNC gene symbols reported in the VEP annotation files, example input would look like this:

```text
ATM 
ATR 
BRCA1 
```
The rules for the generation of the pandas query expression that is passed with the `filtering_expression` parameter are the same as outlined in i. If a gene list is provided without a filtering expression, the applet/app will throw an error. Very briefly, this is because testing for any variants in a gene would include many common variants and would not answer a biologically meaningful question.

Another note of caution is that you generate a list of gene symbols via excel or a similar tool and then save it as txt/csv it is reasonable to check the resulting file for special characters e.g. with `cat -A` or `less` in the Unix command-line, as these are not automatically removed by the underlying Python scripts.

The applet will also provide information on which gene symbols were not found in the VEP files and flag cases where no variants remained after the application of the filtering expression.

The output from the gene list input will be very similar to the files generated when using a SNP list. The key difference is that prefix 'SNP' will be replaced by 'GENE'. This is to allow [mrcepid-runassociationtesting](https://github.com/mrcepid-rap/mrcepid-runassociationtesting)
to recognise that we are running a collapse gene list ather than per-GENE tests.

**Big Note** – As with a SNP list mask generated with a gene list combined with a filtering expression are currently only compatible with the phewas and extract modes of 
[mrcepid-runassociationtesting](https://github.com/mrcepid-rap/mrcepid-runassociationtesting).

### 2. Generating Outputs

This applet then performs a series of formatting steps to generate various output files that are compatible with the
different tools listed in [Background](#background). This README does not go into detail on the format of these files.
Full descriptions of these files as input to individual association tests can be found as part of the repository for 
mrcepid-runassociationtesting:

https://github.com/mrcepid-rap/mrcepid-runassociationtesting

## Running on DNANexus

### Inputs

| input                | description                                                                                                                   |
|----------------------|-------------------------------------------------------------------------------------------------------------------------------|
| filtering_expression | [pandas query](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html) compatible filtering expression.     |
| snplist              | A list of SNPids to generate a mask.                                                                                          |
| genelist             | A list of HGNC gene symbols to generate a mask.                                                                                          |
| file_prefix          | descriptive file prefix for output name                                                                                       |
| bgen_index           | index of bgen information and corresponding annotations **[project-G6BJF50JJv8p4PjGB9yy7YQ2:file-G86GJ3jJJv8fbXVB9PQ2pjz6]**. |

**BIG NOTE**: The name given to 'file_prefix' will be used by the next step in this analysis pipeline as a column in the 
final tab-delimited output files provided for each tool. These columns are derived by splitting `file_prefix` on "-". For 
example, if `file_prefix` is "PTV-HC", where "PTV" is the type of variant and "HC" is the LoFTee setting used to filter PTVs,
two columns in the final associationtesting output will be PTV and HC. As many "dashes" as desired can be included to derive
multiple columns of information (e.g. "MISSENSE-CADD_25-REVEL_07-MAF_01-POPALL" will result in 5 additional columns in the
output of associationtesting). 

**BIG NOTE**: As part of the naming process, associationtesting searches for a special case: where the second column includes
either the keywords "MAF" or "AC" (e.g. "HC_PTV-MAF_01"). This will result in associationtesting naming additional columns
"MASK" and "MAF" rather than generic names.

### Outputs

| output                   | description                                                     |
|--------------------------|-----------------------------------------------------------------|
| output_tarball           | Output tarball containing filtered and processed variant counts |
| log_file                 | Output logfile containing statistics for each                   |

The output_tarball and logfile are named based on the value of `file_prefix` like:

`PTV.tar.gz` and `PTV.log`

While I am not going into detail about the format of the files contained in this tar file, I list here the files for 
record-keeping purposes. All files have a standard prefix identical to that of the tarball with an extra descriptor, with
one set of files for each chromosome:

* <file_prefix>.<chr>.BOLT.bgen – BOLT-ready .bgen file of per-gene 'genotypes'.
* <file_prefix>.<chr>.BOLT.bgen.sample – BOLT-ready bgen sample file of per-gene 'genotypes'.
* <file_prefix>.<chr>.SAIGE.bcf – BCF file for SAIGE of genotypes per-variant.
* <file_prefix>.<chr>.SAIGE.bcf.csi – BCF file index for SAIGE of genotypes per-variant.
* <file_prefix>.<chr>.SAIGE.groupFile.txt – TSV file of variants assigned to each gene. See the [SAIGE documentation](https://github.com/weizhouUMICH/SAIGE/wiki/Genetic-association-tests-using-SAIGE) for more details.
* <file_prefix>.<chr>.STAAR.matrix.rds - STAAR-ready matrix sample - variant - genotype sets in R .rds format
* <file_prefix>.<chr>.variants_table.STAAR.tsv - per-variant annotations for STAAR

The `.log` file contains information on total number of variants, number of variants per-participant, annotations, and a 
histogram of allele counts. Please see this file for more information. 

If providing a SNP list, the only difference for the above is that a single set of files with the prefix <file_prefix>.SNP
is generated. The log file is identical except for a section that documents variant IDs that were NOT found during 
filtering.

If we provide a gene list combined with a pandas query expression the only relevant difference is that a single set of files 
with the prefix <file_prefix>.GENE is generated. The log file is again identical except for a section in the beginning that lists how many variants that fulfil the filtering criteria were found in each gene. Also, in the log file that is automatically generated by DNAnexus (accessible from the "Monitor" tab in your project) information is given on genes that for a range of reasons did not contribute any variants to the analysis (see section 1 iii above).

### Command line example

If this is your first time running this applet within a project other than "MRC - Variant Filtering", please see our
organisational documentation on how to download and build this app on the DNANexus Research Access Platform:

https://github.com/mrcepid-rap

Running this command is fairly straightforward using the DNANexus SDK toolkit:

1. If using a pandas filtering expression:

```commandline
dx run mrcepid-collapsevariants --priority low --destination collapsed_variants/
        -ifiltering_expression='FILTER=="PASS" & AF<=0.001 & LOFTEE=="HC" & PARSED_CSQ=="PTV"' \
        -ifile_prefix="HC_PTV-MAF_01"
```

2. If using a SNP list:

```commandline
dx run mrcepid-collapsevariants --priority low --destination collapsed_variants/
        -isnplist=file-G9B297QJJv8Q3KVP3yJ2v945 \
        -ifile_prefix="Xchr_SNPs"
```

3. If using a gene list combined with a filtering expression:

```commandline
dx run mrcepid-collapsevariants --priority low --destination collapsed_variants/ \
        -ifiltering_expression='FILTER=="PASS" & AF<=0.001 & LOFTEE=="HC" & PARSED_CSQ=="PTV" & CADD>25 & gnomAD_AF < 0.001' \
        -igenelist=file-GBk2pf8JZGX713bjBZ61v72B \
        -ifile_prefix="brca1_brca2_atr"

```


Brief I/O information can also be retrieved on the command line:

```commandline
dx run mrcepid-collapsevariants --help
```

I have set a sensible (and tested) default for compute resources on DNANexus that is baked into the json used for building 
the app (at `dxapp.json`) so setting an instance type is unnecessary. This current default is for a mem3_ssd1_v2_x32 instance
(32 CPUs, 256 Gb RAM, 1200Gb storage). If necessary to adjust compute resources, one can provide a flag like `--instance-type mem1_ssd1_v2_x4`.

#### Batch Running

This tool is not compatible with batch running. All processes are parallelised internally.
