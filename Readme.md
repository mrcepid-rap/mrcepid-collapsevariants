# CollapseVariants (DNAnexus Platform App)

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

### Table of Contents

- [Introduction](#introduction)
  * [Changelog](#changelog)
  * [Background](#background)
  * [Dependencies](#dependencies)
- [Methodology](#methodology)
  * [1. Filtering Variants](#1-filtering-variants)
    + [Filtering With Pandas Query Expressions](#filtering-with-pandas-query-expressions)
    + [Filtering with Variant IDs](#filtering-with-variant-ids)
    + [Filtering with HGNC gene symbols and query expression](#filtering-with-hgnc-gene-symbols-and-query-expression)
  * [2. Generating Outputs](#2-generating-outputs)
- [Running on DNANexus](#running-on-dnanexus)
  * [Inputs](#inputs)
    + [BGEN Index Format](#bgen-index-format)
  * [Outputs](#outputs)
    + [Output Log](#output-log)
    + [Output Tarball](#output-tarball)
      - [Tool-Specific Outputs](#tool-specific-outputs)
        * [BOLT](#bolt)
        * [SAIGE](#saige)
        * [REGENIE](#regenie)
        * [STAAR](#staar)
  * [Command line example](#command-line-example)

## Introduction

This applet generates filtered variant data necessary to perform rare variant burden testing. bcf / vcf and bgen
manipulation is done using either native python or the [bgen](https://github.com/jeremymcrae/bgen) package.

This README makes use of DNANexus file and project naming conventions. Where applicable, an object available on the
DNANexus platform has a hash ID like:

* file – `file-1234567890ABCDEFGHIJKLMN`
* project – `project-1234567890ABCDEFGHIJKLMN`

Information about files and projects can be queried using the `dx describe` tool native to the DNANexus SDK:

```commandline
dx describe file-1234567890ABCDEFGHIJKLMN
```

### Changelog

* v2.0.0
  * Major overhaul of the applet to support WGS. A list of major changes is included below; for more information, please 
    see various pull requests. **These changes are breaking** and will require v2.0.0 of the [mrcepid-runassociationtesting](https://github.com/mrcepid-rap/mrcepid-runassociationtesting) 
    applet to use outputs from this applet for associationtests. 
    * Implemented testing across the tool using simulated data
    * Added support for WGS data
      * .bgen files are no longer tied to chromosomes; they can be 'chunks' of the genome with multiple files representing variant data for a single chromosome
        * **Note**: This means that genes CANNOT cross chunk boundaries! Please see the documentation in [makebgen](https://github.com/mrcepid-rap/mrcepid-makebgen) for information on how to chunk the genome
      * Modified the flow of information through the applet; parallelized the majority of the applet to function on individual .bgen files rather than chromosomes
      * To avoid excessively large outputs, individual variant-level data is no longer stored in output tarballs, except for those used by BOLT.
      * REGENIE-specific outputs are now included in outputs
    * Refactor of most of the codebase to be more functional
    * Several output file formats from v1.* are no longer present. Please refer to the documentation for changes.
    * Removal of all external system calls to 3rd party software. All file manipulation is now done in native python with additional library support.
  * Updated package manager to uv

* v1.2.2
    * Bugfix release to solve issues with new version of plink2
        * This removes plink2 from most aspects of the applet due to issues with silent crashes and formatting problems
        * This may slow down the processing of data due to the need to use qctool rather than plink2 – please report
          excessive slowdown!

* v1.2.1
    * Removed run_cmd from all code due to deprecation
      in [general_utilities](https://github.com/mrcepid-rap/general_utilities)

* v1.2.0
    * The vast majority of changes in this release are invisible to the day-to-day user but should greatly improve
      maintainability of the code
    * Deleted the collapse_resources library and instead implemented general_utilities
        * Docker commands are now run in a consistent way to other applets in this project
    * Modernised the codebase to be more in line with the mrcepid project code style
        * Cleaned up package imports
        * Added python docstrings to all methods
        * Implemented the MRCLOGGER to print information in a more synergistic way with the DNANexus platform rather
          than using `print()` commands
    * Implemented tests using the [mrcepid-testing](https://github.com/mrcepid-rap/mrcepid-testing) suite
        * Please see the developer README (Readme.developer.md) for more details
        * This adds two new command-line options: `testing_script` and `testing_directory`. **DO NOT** use these options
          unless you know what you are doing!
    * Added a `write_string()` method to CollapseLOGGER
    * Ensured that items written to log_file output would be truncated rather than print very long strings
    * Refactored the `filter_bgen()` method into filtering.py to enable easier testing
    * Filter expressions / gene lists / and SNP lists are now actually checked for proper input combinations
    * Created a new method that performs final logging capabilities
    * Modified code to be more in-line with `pandas` v2.0 to allow eventual version change

* v1.1.0
    * Did a major refactor of the codebase to implement object-oriented style for code maintainability.
    * Code is functionally identical from the user's perspective

* v1.0.0
    * Initial numbered release. Changes going forward will be tracked in this section of the documentation

### Background

Downstream of this applet, we have implemented four tools / methods for rare variant burden testing:

* [BOLT](https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html)
* [REGENIE](https://rgcgithub.github.io/regenie/)
* [SAIGE-GENE](https://github.com/weizhouUMICH/SAIGE/wiki/Genetic-association-tests-using-SAIGE)
* [STAAR](https://github.com/xihaoli/STAAR)
* GLMs – vanilla linear/logistic models implemented with python's [statsmodels module](https://www.statsmodels.org/stable/index.html)

These five tools / methods require very different input files to run. The purpose of this applet is to generate inputs
that are compatible with each of these tools input requirements. For more information on the format of these input
files, please see the [outputs](#outputs) section of this README.

### Dependencies

Due to how the DNANexus platform works, this applet is dependent on itself. In short, this means that at launch,
the applet will automatically download the latest version of itself from GitHub and install itself, and other required 
python dependencies, via uv. This allows the module subpackages (in `collapsevariants`) to be imported by the main class.

## Methodology

This applet is step 4 (mrcepid-collapsevariants) of the rare variant testing pipeline developed by Eugene Gardner for the
UKBiobank
RAP at the MRC Epidemiology Unit:

![](https://github.com/mrcepid-rap/.github/blob/main/images/RAPPipeline.v3.png)

This applet has two major steps:

1. Select variants from a filtered and annotated VCF file using either:
    1. a [pandas query](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html) filtering expression
    2. A list of variant IDs
    3. A list of HGNC gene symbols (one per line) combined with a filtering expression (see i.)
2. Generate various output files in varying formats per-chromosome to be fed into
   [mrcepid-runassociationtesting](https://github.com/mrcepid-rap/mrcepid-runassociationtesting).

For more details please see the commented source code available in the `collapsevariants/` python package in this repository.

### 1. Filtering Variants

#### Filtering With Pandas Query Expressions

The user of this applet can provide a filtering expression that is compatible with the pandas query function to select
variants for association testing. For details and tutorials on how to construct such queries, please see the
[pandas query documentation](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html). Briefly, one can
construct various expressions to generate variants that they want to test during rare variant burden tests.
These expressions **MUST** be based on fields within the annotation's tsv file provided as output of
the [mrcepid-makebgen](https://github.com/mrcepid-rap/mrcepid-makebgen) step of this pipeline. 
For possible consequences (PARSED_CSQ) to filter on, please see:

https://github.com/mrcepid-rap/mrcepid-filterbcf#4-parsing-vep-consequences

**Note:** One can also filter based on raw [VEP consequence](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences) using the "CSQ" field, if desired (e.g. stop_gained).

For example, lets say that I want to test the effects of rare (MAF < 1x10<sup>-3</sup>) Protein Truncating Variants (
PTVs) on some phenotype. Thus, I could construct a filtering expression like:

`AF<=0.001 & PARSED_CSQ=="PTV"`

Note that this expression would retain variants that FAILED quality control and did not pass basic PTV QC filters (e.g.
[LoFTEE](https://github.com/konradjk/loftee)). So, we could further modify this expression like so:

`FILTER=="PASS" & AF<=0.001 & LOFTEE=="HC" & PARSED_CSQ=="PTV"`

* FILTER=="PASS" means we only retain variants that passed quality control
* LOFTEE=="HC" retains only PTVs that are LoFTEE high-confidence LoFs

This filtering expression can be increased in complexity to generate a more stringent set of variants:

`FILTER=="PASS" & AF<=0.001 & LOFTEE=="HC" & PARSED_CSQ=="PTV" & REVEL>0.5 & gnomAD_AF < 0.001`

And so forth...

These expressions can be mixed and modified according to user need. In summary, the above expression is run by the
applet using the pseudocode:

```python
import pandas as pd
import gzip

variant_index = pd.read_csv(gzip.open('470k_vep.sorted.tsv.gz', 'rt'), sep="\t")
variant_index = variant_index.query(
    'FILTER=="PASS" & AF<=0.001 & LOFTEE=="HC" & PARSED_CSQ=="PTV" & REVEL>0.5 & gnomAD_AF < 0.001')
```

The table `variant_index` is used for the next step.

#### Filtering with Variant IDs

The user can also opt to provide a list of variant IDs that they want to use to build a single mask using the `snplist`
input. This mode will create a single 'Gene' with the ID ENST00000000000 that includes information for all variants.
Variants listed can be from multiple chromosomes or a single gene. variant IDs **MUST** match the format found in the VEP
annotation files generated by [mrcepid-makebgen](https://github.com/mrcepid-rap/mrcepid-makebgen) applet:

```text
1_1234_A_G
2_4567_ATA_A
```

Any number of variants can be included in this mask but we have not tested the limits in terms of compute on generating
a mask with a large number of variants. This list can also include variants NOT found in the VEP files (e.g. from a
CRISPR screen). This applet will report the variants NOT found in the log file. See [outputs](#outputs) for more information on
this file.

Use of a SNP list will generate similar files to those generated when using a pandas-compatible filtering expression.
The only major difference is that instead of a per-chromosome set of files, a single set of files with the prefix 'SNP' will
be generated. This is to allow [mrcepid-runassociationtesting](https://github.com/mrcepid-rap/mrcepid-runassociationtesting) to 
recognise that we are running collapsed SNPs rather than per-GENE tests.

**Big Note** – Masks generated with a SNP-list are currently only compatible with the phewas and extract modes of
[mrcepid-runassociationtesting](https://github.com/mrcepid-rap/mrcepid-runassociationtesting).

#### Filtering with HGNC gene symbols and query expression

Another option is to provide a list HGNC gene symbols **combined** with a pandas query expression to create a custom
mask. The relevant parameter to specify the gene list is `genelist`. As with the SNP list input we create a single 'GENE'
which here has the ID ENST99999999999. Many considerations for the SNP list input also apply to gene list input (e.g. maximum number of underlying variants
etc.). Gene IDs are checked against the HGNC gene symbols reported in the VEP annotation files, example input would look like this:

```text
ATM 
ATR 
BRCA1 
```

The rules for the generation of the pandas query expression that is passed with the `filtering_expression` parameter are
the same as outlined in i. If a gene list is provided without a filtering expression, the applet/app will throw an
error. Very briefly, this is because testing for any variants in a gene would include many common variants and would not answer
a biologically meaningful question.

Another note of caution is that if you generate a list of gene symbols via excel or a similar tool and then save it as
txt/csv it is reasonable to check the resulting file for special characters e.g. with `cat -A` or `less` in the Unix
command-line, as these are not automatically removed by the underlying Python scripts.

The applet will also provide information on which gene symbols were not found in the VEP files and flag cases where
no variants remained after the application of the filtering expression.

The output from the gene list input will be very similar to the files generated when using a SNP list. The key
difference is that prefix 'SNP' will be replaced by 'GENE'. This is to
allow [mrcepid-runassociationtesting](https://github.com/mrcepid-rap/mrcepid-runassociationtesting)
to recognise that we are running a collapse gene list ather than per-GENE tests.

**Big Note** – As with a SNP list mask generated with a gene list combined with a filtering expression are currently
only compatible with the phewas and extract modes of [mrcepid-runassociationtesting](https://github.com/mrcepid-rap/mrcepid-runassociationtesting).

### 2. Generating Outputs

This applet then performs a series of formatting steps to generate various output files that are compatible with the
different tools listed in [Background](#background). This README does not go into detail on the format of these files.
Full descriptions of these files as input to individual association tests can be found as part of the repository for
mrcepid-runassociationtesting:

https://github.com/mrcepid-rap/mrcepid-runassociationtesting

## Running on DNANexus

### Inputs

| input                | description                                                                                                                                                                                                       |
|----------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| filtering_expression | [pandas query](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html) compatible filtering expression.                                                                                         |
| snplist              | A file containing a list of SNPids to generate a mask.                                                                                                                                                            |
| genelist             | A file containing a list of HGNC gene symbols to generate a mask.                                                                                                                                                 |
| output_prefix        | descriptive file prefix for output name.                                                                                                                                                                          |
| bgen_index           | index of bgen information and corresponding annotations.                                                                                                                                                          | 

**BIG NOTE**: The name given to 'output_prefix' will be used by the next step in this analysis pipeline as a column in the
final tab-delimited output files provided for each tool. These columns are derived by splitting `output_prefix` on "-".
For example, if `output_prefix` is "PTV-HC", where "PTV" is the type of variant and "HC" is the LoFTee setting used to filter
PTVs, two columns in the final associationtesting output will be PTV and HC. As many "dashes" as desired can be included to
derive multiple columns of information (e.g. "MISSENSE-CADD_25-REVEL_07-MAF_01-POPALL" will result in 5 additional columns in
the output of associationtesting).

**BIG NOTE**: As part of the naming process, associationtesting searches for a special case: where the second column
includes either the keywords "MAF" or "AC" (e.g. "HC_PTV-MAF_01"). This will result in associationtesting naming additional
columns "MASK" and "MAF" rather than generic names.

#### BGEN Index Format

The BGEN index must follow a strict format and is a required option:

```
prefix   vep_dxid   bgen_dxid    bgen_index_dxid   sample_dxid
1_1_100000    file-1234567890ABCDEFGH   file-0987654321ABCDEFGH   file-1234567890HGFEDCBA   file-0987654321HGFEDCBA
```

Where prefix can be any alphanumeric string that provides a unique identifier to a given bgen. A more detailed 
description of what this file looks like and how to make it is included in the following repository:

https://github.com/mrcepid-rap/QC_workflow

### Outputs

| output         | description                                                     |
|----------------|-----------------------------------------------------------------|
| output_tarball | Output tarball containing filtered and processed variant counts |
| log_file       | Output logfile containing summary statistics                    |

The output_tarball and logfile are named based on the value of `output_prefix` like:

`<output_prefix>.tar.gz` and `<output_prefix>.log`

#### Output Log

The `.log` file contains information on total number of variants, number of variants per-participant, annotations, and a
histogram of allele counts. Please see this file for more information.

For SNP masks, the log file is identical except for a section that documents variant IDs that were NOT found during
filtering.

For GENE masks, The log file is again identical except for a section in the beginning that lists how many variants that
fulfil the filtering criteria were found in each gene. Also, in the log file that is automatically generated by DNAnexus
(accessible from the "Monitor" tab in your project) information is given on genes that for a range of reasons did not
contribute any variants to the analysis (see [above](#filtering-with-hgnc-gene-symbols--query-expression) for more information).

#### Output Tarball

Within `<output_prefix>.tar.gz`, all output files have a standard name that composed of two parts: 1) the requested `output_prefix` and 2) for each file 
provided to `bgen_index`, a file-specific prefix. For each .bgen file provided to the applet, the following files are generated:

* <output_prefix>.<prefix>.BOLT.bgen – BOLT-ready .bgen file of per-gene 'genotypes'.
* <output_prefix>.<prefix>.BOLT.bgen.bgi – BOLT-ready .bgen.bgi bgen index of per-gene 'genotypes'.
* <output_prefix>.<prefix>.BOLT.sample – BOLT-ready bgen sample file of per-gene 'genotypes'.
* <output_prefix>.<prefix>.SAIGE.groupFile.txt – TSV file of variants assigned to each gene.
* <output_prefix>.<prefix>.REGENIE.annotationFile.txt – TSV of variant assigned to gene and mask.
* <output_prefix>.<prefix>.REGENIE.setListFile.txt – TSV of variants assigned to genes with coordinate information
* <output_prefix>.<prefix>.REGENIE.maskfile.txt – TSV of masks in the given output
* <output_prefix>.<prefix>.STAAR.samples_table.tsv - STAAR ready sample table
* <output_prefix>.<prefix>.STAAR.variants_table.tsv - per-variant annotations for STAAR

If providing a SNP or GENE list, the following differences for the above apply:

* a single set of files with the prefix `<output_prefix>.SNP.*` or `<output_prefix>.GENE.*`, respectively, is generated.
* An additional [matrixmarket format file](https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.mmwrite.html) is generated for STAAR called <output_prefix>.<prefix>.STAAR.mtx 

##### Tool-Specific Outputs

###### BOLT

* <output_prefix>.<prefix>.BOLT.bgen – BOLT-ready .bgen file of per-gene 'genotypes'.

BOLT does not allow for gene collapsing tests by default. To enable burden tests with BOLT we generate a dummy .bgen
file that collapses all variants in a given gene-mask pair into a set of biallelic variants. In brief, individuals 
with a variant specified by the current mask are assigned a heterozygous genotype (0/1) and individuals without the 
variant are assigned a homozygous reference genotype (0/0). Resulting genotypes are then encoded into bgen v1.2 
format. As .bgen can only encode nucleotide information as REF / ALT, we encode each gene mask with REF = A and ALT = C.
*This bgen file is REF FIRST*. To be clear, this was an arbitrary choice to ensure consistency throughout this 
pipeline. ID for each variant is encoded as the *ENST* for the current transcript. The text representation of this 
variant is the following (roughly in bgen specification format):

```text
varID   rsID    chr pos ref alt genotype_block
ENST00000000001 ENST00000000001 1 100000 A C <genotypes>
```

* <output_prefix>.<prefix>.BOLT.bgen.bgi

An index file to enable random access of the .bgen file. Please see [bgenix](https://enkre.net/cgi-bin/code/bgen/doc/trunk/doc/wiki/bgenix.md)
for more information.

* <output_prefix>.<prefix>.BOLT.sample – BOLT-ready bgen sample file of per-gene 'genotypes'.

A standard .bgen .sample file. To be clear, this is *NOT* the sample-v2 file. An example of the format is

```text
ID_1 ID_2 missing sex
0   0   0   D
1000000 1000000 0 NA
1000001 1000001 0 NA
```

###### SAIGE

* <output_prefix>.<prefix>.SAIGE.groupFile.txt – TSV file of variants assigned to each gene.

SAIGE only requires a single group file that assigns variants to genes and then variants to masks. This file is generated
as a TSV file with the following no-header format:

```text
ENST00000000001 var 1:100000:C:A    1:100010:T:G    1:100020:G:A
ENST00000000001 anno    foo foo foo
ENST00000000002 var 1:100100:T:A    1:100110:G:C    1:100200:C:A
ENST00000000002 anno    foo foo foo
```

Note that the variant ID is the same as the one used in other variant files. SAIGE requires that the variant ID is encoded
as `chr:pos:REF:ALT`. This is opposed to that stored in the vep.tsv files of `chr_pos_REF_ALT`. The value of 'foo' 
in the annotation line is a dummy placeholder and is used to by SAIGE when running burden tests. Do not change `foo` unless 
you also change the burden module. In theory, it can be the name of the mask, but does not make a difference to the downstream
workflow.

###### REGENIE

* <output_prefix>.<prefix>.REGENIE.annotationFile.txt – TSV of variant assigned to gene and mask.

This file assigns individual variants to genes and a mask. The format is a TSV file with the following no-header format:

```text
1_100000_C_A    ENST00000000001    HC_PTV-MAF_01
1_100010_T_G    ENST00000000001    HC_PTV-MAF_01
1_100020_G_A    ENST00000000001    HC_PTV-MAF_01
1_100100_T_A    ENST00000000002    HC_PTV-MAF_01
1_100110_G_C    ENST00000000002    HC_PTV-MAF_01
1_100200_C_A    ENST00000000002    HC_PTV-MAF_01
```

* <output_prefix>.<prefix>.REGENIE.setListFile.txt – TSV of variants assigned to genes with coordinate information

This file assigns individuals variants to a gene and provides gene-level coordinate information. This file is 
similar to that generated for SAIGE, except it also encodes chromosome and position and separates variant IDs by comma 
(`,`) rather than space. The format is a TSV file with the following no-header format:

```text
ENST00000000001 1   100000  1_100000_C_A,1_100010_T_G,1_100020_G_A
ENST00000000002 1   100100  1_100100_T_A,1_100110_G_C,1_100200_C_A
```

Where rows are ENST, chromosome, position, and a comma-separated list of variant IDs, respectively.

* <output_prefix>.<prefix>.REGENIE.maskfile.txt – TSV of masks in the given output

This file just lists possible masks in a given run. It is a TSV format file, but should only ever contain ONE line per
chromosome-mask pair. The format is:

```text
HC_PTV-MAF_01   HC_PTV-MAF_01
```

###### STAAR

Raw variant-level genetic data is NOT stored in the output files; however, variant-level data _is_ required for STAAR
to perform association testing. Therefore, we create the column (variants_table.tsv) and rows (samples_table.tsv) files
of a sparse matrix to describe what data should be filled into the matrix. This data is then generated on-the-fly at 
the time of running the association test. The format of these files is as follows:

* <output_prefix>.<prefix>.STAAR.variants_table.tsv – TSV of variant assigned to columns

```text
1_100000_C_A    1    100000  ENST00000000001 1
1_100010_T_G    1    100010  ENST00000000001 2
1_100020_G_A    1    100020  ENST00000000001 3
1_100100_T_A    1    100100  ENST00000000002 4
1_100110_G_C    1    100110  ENST00000000002 5
1_100200_C_A    1    100200  ENST00000000002 6
```

Where columns are variant ID, chromosome, position, gene ID, and column number [i] of the matrix, respectively.

* <output_prefix>.<prefix>.STAAR.samples_table.tsv - TSV of sample IDs assigned to rows

```text
1000000 1
1000001 2
1000002 3
1000003 4
```

Where columns are sample ID and row number [j] of the matrix, respectively.

* <output_prefix>.<prefix>.STAAR.mtx - MatrixMarket format file of the sparse matrix

**Note:** Only for SNP & GENE masks!

```text
%%MatrixMarket matrix coordinate integer general
%
10000 34 128
211 26 1
243 16 1
...
956 31 1
```

### Command line example

There are two ways to acquire this applet:

1. As an **applet** – clone the repository from github and `dx build` an APPLET into your own workspace. If this is your
   first time doing this within a project other than "MRC - Variant Filtering", please see our organisational documentation on how to
   download and build this app on the DNANexus Research Access Platform:

https://github.com/mrcepid-rap

2. As an **app** – use the app that has been provided in the DNANexus global namespace. This will ensure you are always
   using the latest version and keeps you from having to manually update your local version.

**Note:** All commands below have been provided as if using option (2) above!

Running this command is fairly straightforward using the DNANexus SDK toolkit:

1. If using a pandas filtering expression:

```commandline
dx run app-mrcepid-collapsevariants --priority low --destination collapsed_variants/
        -ifiltering_expression='FILTER=="PASS" & AF<=0.001 & LOFTEE=="HC" & PARSED_CSQ=="PTV"' \
        -ioutput_prefix="HC_PTV-MAF_01" \
        -ibgen_index=file-1234567890ABCDEFGH
```

2. If using a SNP list:

```commandline
dx run app-mrcepid-collapsevariants --priority low --destination collapsed_variants/
        -isnplist=file-G9B297QJJv8Q3KVP3yJ2v945 \
        -ioutput_prefix="Xchr_SNPs" \
        -ibgen_index=file-1234567890ABCDEFGH
```

3. If using a gene list combined with a filtering expression:

```commandline
dx run app-mrcepid-collapsevariants --priority low --destination collapsed_variants/ \
        -ifiltering_expression='FILTER=="PASS" & AF<=0.001 & LOFTEE=="HC" & PARSED_CSQ=="PTV" & CADD>25 & gnomAD_AF < 0.001' \
        -igenelist=file-GBk2pf8JZGX713bjBZ61v72B \
        -ioutput_prefix="brca1_brca2_atr" \
        -ibgen_index=file-1234567890ABCDEFGH
```

Brief I/O information can also be retrieved on the command line:

```commandline
dx run mrcepid-collapsevariants --help
```

I have set a sensible (and tested) default for compute resources on DNANexus that is baked into the json used for
building the app (at `dxapp.json`) so setting an instance type is unnecessary. This current default is for a mem3_ssd1_v2_x32
instance (32 CPUs, 256 Gb RAM, 1200Gb storage). If necessary to adjust compute resources, one can provide a flag like
`--instance-type mem1_ssd1_v2_x4`.