#
# Author: Eugene Gardner (eugene.gardner at mrc.epid.cam.ac.uk)
#
# Prior to using this script PLEASE follow the instructions in the developer readme (Readme.developer.md) carefully.
# This Readme provides instructions on how to regenerate testing data necessary to run these tests.

import os
import time
import json
from typing import Optional, Dict

import dxpy
import pytest
import filecmp
import pandas as pd

from pathlib import Path

from general_utilities.association_resources import fix_plink_bgen_sample_sex
from general_utilities.job_management.command_executor import build_default_command_executor

from collapsevariants.ingest_data import IngestData
from collapsevariants.collapse_logger import CollapseLOGGER
from collapsevariants.snp_list_generator import SNPListGenerator
from collapsevariants.snp_list_merger import SNPMerger
from general_utilities.mrc_logger import MRCLogger

test_folder = Path(os.getenv('TEST_DIR'))
CMD_EXEC = build_default_command_executor()
LOGGER = MRCLogger(__name__).get_logger()


def get_testing_files(file_name: str, download: bool = False) -> dict:
    """A function to pull named testing files from the provided test folder onto the local storage of this instance

    Some testing needs to involve downloading / uploading of files from the DNANexus platform. However, testing does
    not allow us to provide static DNANexus file IDs. Thus, this method is just a wrapper around
    dxpy.find_one_data_object to enable us to find files by name in the folder of test data used for running tests.
    The location of all files using this method is hard-coded to the test folder provided to testing at startup.

    :param file_name: String name of the file we want
    :param download: Download the file
    :return: A 'mock' DNANexus link in the format of a dictionary with key always as '$dnanexus_link' and value as a
        DNANexus file-id
    """
    found_file = dxpy.find_one_data_object(classname='file',
                                           project=dxpy.PROJECT_CONTEXT_ID,
                                           name_mode='exact',
                                           name=file_name,
                                           folder=f'{test_folder}',
                                           zero_ok=False)
    if download:
        dxpy.download_dxfile(dxid=found_file['id'], project=found_file['project'], filename=file_name)

    return {'$dnanexus_link': found_file['id']}


expression = 'AF<0.001 & PARSED_CSQ=="PTV" & LOFTEE=="HC" & FILTER=="PASS"'
bgen_index = {'$dnanexus_link': os.getenv('BGEN_INDEX')}
snplist_dxlink = get_testing_files('snp_list.pheno_1_MISS.txt')
genelist_dxlink = get_testing_files('gene_list.pheno_1_MISS.txt')
@pytest.mark.parametrize(
    argnames=['testing_expression', 'snplist', 'genelist', 'expected_exception'],
    argvalues=zip([expression, None, expression, expression, None, None],
                  [None, snplist_dxlink, None, snplist_dxlink, None, snplist_dxlink],
                  [None, None, genelist_dxlink, None, genelist_dxlink, genelist_dxlink],
                  [None, None, None, ValueError, ValueError, ValueError])
)
def test_ingest_data(testing_expression: Optional[str], snplist: Optional[Dict], genelist: Optional[Dict],
                     expected_exception: Exception):
    """Test IngestData in CollapseVariants

    Test in order:
    1. filtering expression
    2. SNPList
    3. gene list + filtering expression

    Errors:
    4. SNPList + expression
    5. No expression + gene list
    6. snp list + gene list

    :param testing_expression: A filtering expression to use for this test
    :param snplist: An optional dxpy.DXLink to a SNP list
    :param genelist: An optional dxpy.DXLink to a Gene list
    :param expected_exception: Exception that should be thrown by this combination of elements
    """

    # When running this code successive times, 'docker pull' functionality will just assert that the image is already
    # on the instance, and no need to test...
    if expected_exception is None:
        ingested_data = IngestData(bgen_index,
                                   testing_expression,
                                   snplist,
                                   genelist)
        assert sorted([f'{x}' for x in range(1, 23)] + ['X']) == sorted(ingested_data.bgen_index.keys())
        assert ingested_data.filtering_expression == testing_expression
        if snplist:
            assert ingested_data.snp_list_path.exists()
        else:
            assert ingested_data.snp_list_path is None
        if genelist:
            assert ingested_data.gene_list_path.exists()
        else:
            assert ingested_data.gene_list_path is None
    else:
        with pytest.raises(expected_exception):
            IngestData(bgen_index,
                       testing_expression,
                       snplist,
                       genelist)

get_testing_files('expression_test_data.json', download=True)
var_test = json.load(Path('expression_test_data.json').open('r'))
@pytest.mark.parametrize(
    argnames=['var_info'],
    argvalues=zip(var_test)
)
def test_filtering(var_info: dict):
    """Test both the SNPListGenerator and filter_bgen classes / methods for a filtering expression

    This test essentially performs an end-to-end test for collapsing using just a filtering expression we 're-ingest'
    data and use the resulting files to make sure all files are as they should be for chromosome 1 with the assumption
    that all tests would pass for all chromosomes.

    :param var_info: Dictionary of information necessary for testing. See below for individual parameters
    """

    # Set test parameters from the dictionary of known results provided to this method
    variant_type = var_info['variant_type']  # Specifically one of three PARSED_CSQ fields: PTV, MISS, SYN
    test_type = var_info['test_type']  # type of test we are running
    expected_vars = var_info['expected_vars']  # Total number of expected vars on ALL chromosomes after running a filtering expression
    test_gene = var_info['test_gene']  # Name of the gene that we are using to run tests
    test_var = var_info['test_var']  # Name of the variant that we are using to run tests
    gene_het_count = var_info['gene_het_count']  # Total number of het carriers of `test_gene`
    var_het_count = var_info['var_het_count']  # Total number of het carriers of `test_var`
    tot_gene_count = var_info['tot_gene_count']  # Total number of genes on chromosome 1 for given `variant_type`
    test_var_count = var_info['test_var_count']  # Total number of variants in `test_gene`
    tot_var_count = var_info['tot_var_count']  # Total number of variants on chromosome 1 for given `variant_type`
    snp_list = get_testing_files(var_info['snp_list'], download=False) if var_info['snp_list'] else None # Get a SNP list as a DXID (if requested)
    gene_list = get_testing_files(var_info['gene_list'], download=False) if var_info['gene_list'] else None # Get a gene list as a DXID (if requested)

    # Running the actual processing / filtering first and then doing all tests below
    if not snp_list:
        filtering_expression = f'MAF<0.001 & PARSED_CSQ=="{variant_type}"'
    else:
        filtering_expression = None
    ingested_data = IngestData(bgen_index, filtering_expression, snp_list, gene_list)
    LOG_FILE = CollapseLOGGER(f'{variant_type}_{test_type}_test')
    snp_list_generator = SNPListGenerator(ingested_data, LOG_FILE)

    filtering_total = 0
    # Don't want to waste time doing every chromosome unless a SNP or GENE list
    for chromosome in (snp_list_generator.chromosomes if (gene_list or snp_list) else ['1']):
        LOGGER.info(f'Running filtering for chromosome {chromosome}...')
        chrom_bgen_index = ingested_data.bgen_index[chromosome]
        per_chromosome_total, _, sample_table = filter_bgen(file_prefix=f'{variant_type}_{test_type}_test',
                                                            chromosome=chromosome,
                                                            chrom_bgen_index=chrom_bgen_index,
                                                            cmd_exec=CMD_EXEC)
        filtering_total += per_chromosome_total

    LOGGER.info(f'Filtering done...')

    if snp_list or gene_list:
        SNPMerger(snp_list_generator.chromosomes, f'{variant_type}_{test_type}_test',
                  ingested_data.found_genes, cmd_exec=CMD_EXEC)
    # This is the end of actual processing. All tests follow

    # Test overall site count as determined by snp_list_generator
    assert expected_vars == snp_list_generator.total_sites

    # Test composition of the pandas DataFrame of variants
    expected_cols = ['CHROM', 'POS', 'REF', 'ALT', 'ogVarID', 'FILTER', 'AF', 'F_MISSING', 'AN', 'AC', 'MANE', 'ENST',
                     'ENSG', 'BIOTYPE', 'SYMBOL', 'CSQ', 'gnomAD_AF', 'CADD', 'REVEL', 'SIFT', 'POLYPHEN', 'LOFTEE',
                     'AA', 'AApos', 'PARSED_CSQ', 'MULTI', 'INDEL', 'MINOR', 'MAJOR', 'MAF', 'MAC']
    assert snp_list_generator.variant_index.index.name == 'varID'
    assert snp_list_generator.variant_index.columns.to_list() == expected_cols
    assert expected_vars == len(snp_list_generator.variant_index)

    # Test contents of the SNP lists
    with Path('include_snps.txt').open('r') as snp_list,\
            Path('snp_ENST.txt').open('r') as gene_list:
        total_snps = 0
        for _ in snp_list:
            total_snps += 1

        total_genes = 0
        total_col_errors = 0
        for data in gene_list:
            data = data.rstrip()
            data = data.split('\t')
            if len(data) != 4:
                total_col_errors += 1
            total_genes += 1

        assert total_snps == expected_vars
        assert total_genes == expected_vars + 1  # Gene list has a header
        assert total_col_errors == 0

    # Test contents of the log file:
    get_testing_files(f'{variant_type}_{test_type}_actual.log', download=True)
    assert filecmp.cmp(f'{variant_type}_{test_type}_actual.log', f'{variant_type}_{test_type}_test.log')

    # Now test variant filtering. I think we only have to do one chromosome to prove to ourselves that everything
    # worked? Ultimately just want to make sure all the individual files written by this process do what we think
    # they do. No need for unit tests in this case.
    assert snp_list_generator.total_sites == total_snps  # Ensure the different code parts do the same thing
    assert filtering_total == tot_var_count

    # Using one randomly selected variant / gene pair to ensure proper information exists in the file
    with Path('extract_gene.txt').open('w') as extract_gene:
        extract_gene.write(f'{test_gene}\n')

    # Check if all required files exist
    required_suffixes = ['BOLT.bgen', 'BOLT.sample', 'SAIGE.bcf', 'SAIGE.bcf.csi', 'SAIGE.groupFile.txt',
                         'STAAR.matrix.rds', 'sample', 'variants_table.STAAR.tsv']
    # Suffix of final files differs depending on processing performed... Need to account for that here
    if var_info['snp_list']:
        chrom_suffix = 'SNP'
    elif var_info['gene_list']:
        chrom_suffix = 'GENE'
    else:
        chrom_suffix = '1'
    assert sorted([f'{file}' for file in Path('./').glob(f'{variant_type}_{test_type}_test.{chrom_suffix}.*')]) == \
           sorted([f'{variant_type}_{test_type}_test.{chrom_suffix}.{suffix}' for suffix in required_suffixes])

    # Also spot check the primary genotype / data files
    # BOLT.bgen
    # NOTE: For GENE/SNP lists the number autogenerated by simulate_data.R for gene_het_count will be wrong. I have
    # manually calculated this value and added it to the json.
    fix_sample = fix_plink_bgen_sample_sex(Path(f'{variant_type}_{test_type}_test.{chrom_suffix}.BOLT.sample'))
    test_cmd = f'plink2 --bgen /test/{variant_type}_{test_type}_test.{chrom_suffix}.BOLT.bgen \'ref-last\' ' \
               f'--sample /test/{fix_sample} ' \
               f'--extract /test/extract_gene.txt --export AD --out /test/{variant_type}_{test_type}_bolt ' \
               f'--split-par hg38'
    CMD_EXEC.run_cmd_on_docker(test_cmd)
    bolt_test = pd.read_csv(f'{variant_type}_{test_type}_bolt.raw', sep='\t')
    assert bolt_test[f'{test_gene}_HET'].sum() == gene_het_count  # Alt count correct
    assert bolt_test[f'{test_gene}_A'].sum() == (20000 - gene_het_count)  # Ref count correct
    assert len(bolt_test) == 10000  # Total samples correct

    # SAIGE.bcf
    test_cmd = f'bcftools query -i \'ID == "{test_var}"\' -f \'[%CHROM\\t%POS\\t%REF\\t%ALT\\t%SAMPLE\\t%GT\\n]\' ' \
               f'/test/{variant_type}_{test_type}_test.{chrom_suffix}.SAIGE.bcf > {variant_type}_{test_type}_saige.tsv'
    CMD_EXEC.run_cmd_on_docker(test_cmd)
    saige_test = pd.read_csv(f'{variant_type}_{test_type}_saige.tsv', sep="\t",
                             names=['chrom', 'pos', 'ref', 'alt', 'FID', 'gt'])
    assert saige_test['gt'].value_counts()['0/1'] == var_het_count
    assert saige_test['gt'].value_counts()['0/0'] == 10000 - var_het_count
    assert len(saige_test) == 10000

    # SAIGE.groupFile.txt
    with Path(f'{variant_type}_{test_type}_test.{chrom_suffix}.SAIGE.groupFile.txt').open('r') as saige_group:
        total_genes = 0
        for line in saige_group:
            total_genes += 1
            line = line.rstrip().split('\t')
            if line[0] == test_var:
                assert (len(line) - 1) == test_var_count

        assert total_genes == tot_gene_count

    # variants_table.STAAR.tsv
    staar_variants = pd.read_csv(f'{variant_type}_{test_type}_test.{chrom_suffix}.variants_table.STAAR.tsv', sep='\t')
    assert len(staar_variants.query(f'varID == @test_var')) == 1
    assert len(staar_variants.query('ENST == @test_gene')) == test_var_count
    assert len(staar_variants) == tot_var_count
    assert staar_variants['column'].max() == tot_var_count

# Code to regenerate logs (if necessary)
# Drop into section just above if not snp_list
# logs = [{'variant_type': 'PTV', 'test_type': 'expression', 'snp_list': None, 'gene_list': None},
#         {'variant_type': 'MISS', 'test_type': 'expression', 'snp_list': None, 'gene_list': None},
#         {'variant_type': 'SYN', 'test_type': 'expression', 'snp_list': None, 'gene_list': None},
#         {'variant_type': 'PTV', 'test_type': 'gene', 'snp_list': None, 'gene_list': {'$dnanexus_link': 'file-GQz9330J0zVyXjxgqFzB5Y09'}},
#         {'variant_type': 'MISS', 'test_type': 'gene', 'snp_list': None, 'gene_list': {'$dnanexus_link': 'file-GQz991QJ0zVjkZfV4QyJxjvg'}},
#         {'variant_type': 'PTV', 'test_type': 'snp', 'snp_list': {'$dnanexus_link': 'file-GQz9YzjJ0zVV8X3F63Yyjjzf'}, 'gene_list': None},
#         {'variant_type': 'MISS', 'test_type': 'snp', 'snp_list': {'$dnanexus_link': 'file-GQz9Yv0J0zVfz88KX1p1vgKX'}, 'gene_list': None}]
# for log in logs:
#     variant_type = log['variant_type']
#     test_type = log['test_type']
#     snp_list = log['snp_list']
#     gene_list = log['gene_list']
#     if not snp_list:
#         filtering_expression = f'MAF<0.001 & PARSED_CSQ=="{variant_type}"'
#     else:
#         filtering_expression = None
#     ingested_data = IngestData(bgen_index, filtering_expression, snp_list, gene_list)
#     LOG_FILE = CollapseLOGGER(f'{variant_type}_{test_type}_actual')
#     snp_list_generator = SNPListGenerator(ingested_data, LOG_FILE)
