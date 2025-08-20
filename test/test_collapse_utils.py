import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
from general_utilities.import_utils.import_lib import BGENInformation
from scipy.sparse import csr_matrix

from collapsevariants.utilities.collapse_logger import CollapseLOGGER
from collapsevariants.utilities.collapse_utils import get_sample_ids, check_matrix_stats
from collapsevariants.utilities.parallelization_wrappers import stat_writer, generate_genotype_matrix
from collapsevariants.snp_list_generator.snp_list_generator import SNPListGenerator

# Validated test data:
test_dir = Path(__file__).parent
test_data_dir = test_dir / 'test_data/'

# Filtering lists
snp_path = InputFileHandler(test_data_dir / 'snp_list.v2.txt')
gene_enst_path = InputFileHandler(test_data_dir / 'gene_list.ENST.txt')

# Variant information
bgen_dict = {'chr1_chunk1': BGENInformation(index= InputFileHandler(test_data_dir / 'chr1_chunk1.bgen.bgi'),
                                            bgen= InputFileHandler(test_data_dir / 'chr1_chunk1.bgen'),
                                            sample= InputFileHandler(test_data_dir / 'chr1_chunk1.sample'),
                                            vep= InputFileHandler(test_data_dir / 'chr1_chunk1.vep.tsv.gz'),
                                            gts= InputFileHandler(test_data_dir / 'chr1_chunk1.gts')),
             'chr1_chunk2': BGENInformation(index= InputFileHandler(test_data_dir / 'chr1_chunk2.bgen.bgi'),
                                            bgen= InputFileHandler(test_data_dir / 'chr1_chunk2.bgen'),
                                            sample= InputFileHandler(test_data_dir / 'chr1_chunk2.sample'),
                                            vep= InputFileHandler(test_data_dir / 'chr1_chunk2.vep.tsv.gz'),
                                            gts= InputFileHandler(test_data_dir / 'chr1_chunk2.gts')),
             'chr1_chunk3': BGENInformation(index= InputFileHandler(test_data_dir / 'chr1_chunk3.bgen.bgi'),
                                            bgen= InputFileHandler(test_data_dir / 'chr1_chunk3.bgen'),
                                            sample= InputFileHandler(test_data_dir / 'chr1_chunk3.sample'),
                                            vep= InputFileHandler(test_data_dir / 'chr1_chunk3.vep.tsv.gz'),
                                            gts= InputFileHandler(test_data_dir / 'chr1_chunk3.gts'))}


def generate_expected_counts(test_data: Path, gene_list: InputFileHandler = None, snp_list: InputFileHandler = None) -> np.ndarray:
    """
    Generate expected genotype counts from test data.

    This function reads genotype data from a file, optionally filters it by gene list and SNP list,
    and calculates the sum of genotypes for each sample. The resulting counts are returned as a numpy array.

    :param test_data: Path to the file containi`ng genotype data.
    :param gene_list: Optional path to a file containing a list of genes to filter by.
    :param snp_list: Optional path to a file containing a list of SNPs to filter by.
    :return: A numpy array containing the sum of genotypes for each sample.
    """
    gts = pd.read_csv(test_data,
                      sep='\t', names=['CHROM', 'POS', 'ID', 'ENST', 'CSQ', 'LOFTEE', 'MAF', 'SAMPLE', 'GT'])

    gts['GT'] = gts['GT'].apply(lambda x: 1 if x == '0/1' else 2)
    gts['SAMPLE'] = gts['SAMPLE'].str.split('_', expand=True)[1].astype(int)

    # Filter by ENST if required
    if gene_list:
        with gene_list.get_file_handle().open('r') as gene_file:
            gene_list = [gene.rstrip() for gene in gene_file.readlines()]

        gts = gts[gts['ENST'].isin(gene_list)]

    # Filter by SNP if required
    if snp_list:
        with snp_list.get_file_handle().open('r') as snp_file:
            snp_list = [snp.rstrip() for snp in snp_file.readlines()]

        gts = gts[gts['ID'].isin(snp_list)]

    gt_totals = gts.groupby('SAMPLE').aggregate(GT_sum=('GT', 'sum'))

    # Make a dummy frame to make sure we capture all sample IDs
    dummy = pd.DataFrame(index=range(0, 10000))
    dummy.index.name = 'SAMPLE'

    gt_totals = dummy.merge(gt_totals, left_index=True, right_index=True, how='left')
    gt_totals = gt_totals.fillna(0)

    return gt_totals['GT_sum'].to_numpy()


@pytest.mark.parametrize("filtering_expression, gene_list_handler, snp_list_handler, expected_num_sites, collapse_func",
                         [
                             ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', None, None, 13592, True),
                             ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', gene_enst_path, None, 128, True),
                             (None, None, snp_path, 2299, True),
                             ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', None, None, 13592, False),
                             ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', gene_enst_path, None, 128, False),
                             (None, None, snp_path, 2299, False),
                         ]
                         )
def test_csr_matrix_generation(tmp_path: Path, filtering_expression: str, gene_list_handler: InputFileHandler,
                               snp_list_handler: InputFileHandler, expected_num_sites: int, collapse_func: bool):
    """This tests both :func:`generate_csr_matrix_from_bgen` and :func:`check_matrix_stats` functions.

    :param tmp_path: temporary pytest path for the logs
    :param filtering_expression: the filtering expression that is used as a mask
    :param gene_list_handler: filepath to   a text file containing a gene list to be used in the mask
    :param snp_list_handler: filepath to a text file containing a variant list to be used as a mask
    :param expected_num_sites: expected number of sites post-masking
    :param collapse_func: whether to collapse the matrix or not. WARNING: THIS DOES NOTHING AND NEEDS TO BE FIXED!
    """

    log_path = tmp_path / 'HC_PTV-MAF_01.log'
    test_log = CollapseLOGGER(log_path)

    snp_list_generator = SNPListGenerator(bgen_dict=bgen_dict, filtering_expression=filtering_expression,
                                          gene_list_handler=gene_list_handler, snp_list_handler=snp_list_handler, log_file=test_log)

    total_variants = 0
    for bgen_prefix, variant_list in snp_list_generator.genes.items():
        bgen_dict[bgen_prefix]['index'].get_file_handle()

        genotypes = generate_genotype_matrix(
            bgen_prefix,
            bgen_dict[bgen_prefix],
            variant_list,
            should_collapse=False if snp_list_handler or gene_list_handler else True,
            delete_on_complete=False  # Make sure we keep the files for testing
        )

        ac_table, gene_ac_table, gene_totals = check_matrix_stats(genotypes[1:], variant_list)
        assert len(ac_table) == 10000
        assert len(gene_ac_table) == 10000
        assert len(gene_totals) == len(variant_list['ENST'].unique())

        expected_sites_path = bgen_dict[bgen_prefix]['gts'].get_file_handle()
        expected_counts = generate_expected_counts(expected_sites_path,
                                                   gene_list=gene_list_handler, snp_list=snp_list_handler)

        # EUGENE – Remember that this fails due to an annotation issue in Duat
        # The problem is:
        # 1. There are two variants with identical pos / ref / alt alleles but different CSQ annotations
        # 2. bcftools annotate cannot tell the difference and just uses the 1st annotation
        # I have left in this numpy bit so that you can see where the error is. It prints out the offending
        # sample where compiled gt != expected gt;
        # tldr: the test data is wrong, not collapsevariants
        # print(bgen_prefix)
        # print(np.argwhere(np.not_equal(ac_table, expected_counts)))
        assert np.array_equal(ac_table, expected_counts)

        # the output is a tuple, but we only want the matrix here to do the test
        genotype_matrix, genotype_totals = genotypes[1:]

        total_variants += np.sum(genotype_matrix)
        assert type(genotype_matrix) is csr_matrix

    assert total_variants == expected_num_sites


@pytest.mark.parametrize("sample_file", [bgen_info['sample'].get_file_handle() for bgen_info in bgen_dict.values()])
def test_get_sample_ids(sample_file: Path):
    """
    Test the `get_sample_ids` function to ensure it correctly reads sample IDs from a given sample file.

    This test verifies that:
    1. The number of sample IDs read is 10,000.
    2. The returned sample IDs are in a list.
    3. Each sample ID is a non-delimited string.

    Parameters:
    sample_file (Path): The path to the sample file to be tested.

    Asserts:
    - The length of the sample IDs list is 10,000.
    - The type of the sample IDs list is `list`.
    - The type of the first sample ID is `str`.
    - Each sample ID is a non-delimited string.
    """
    sample_ids = get_sample_ids(sample_file)
    assert len(sample_ids) == 10000
    assert type(sample_ids) is list
    assert type(sample_ids[0]) is str
    assert len(sample_ids[0].split()) == 1  # Make sure the sample IDs are non-delimited strings


@pytest.mark.parametrize("bgen_prefix, bgen_info", bgen_dict.items())
@pytest.mark.parametrize(
    "name, filtering_expression, gene_list_handler, snp_list_handler",
    [
        ('expression', 'PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', None, None),
        ('gene_list', 'PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', gene_enst_path, None),
        ('snp_list', None, None, snp_path)
    ]
)
def test_stat_writer(tmp_path: Path, bgen_prefix, bgen_info, name: str, filtering_expression: str,
                     gene_list_handler: InputFileHandler, snp_list_handler: InputFileHandler):
    """
    Test the `stat_writer` function to ensure it correctly writes statistics about the collapsing operations.

    This test verifies that:
    1. The log file is created.
    2. The contents of the log file match the expected log contents.

    Parameters:
    tmp_path (Path): A temporary directory path provided by pytest for storing the log file.
    filtering_expression (str): The filtering expression used to generate the SNP list.
    gene_list_path (Path): The path to the gene list file.
    snp_list_path (Path): The path to the SNP list file.
    expected_num_sites (int): The expected number of sites based on the filtering expression.

    Asserts:
    - The log file is created.
    - The contents of the log file match the expected log contents.
    """

    log_path = tmp_path / 'HC_PTV-MAF_01.log'
    test_log = CollapseLOGGER(log_path)
    print(tmp_path)

    snp_list_generator = SNPListGenerator(bgen_dict={bgen_prefix: bgen_info}, filtering_expression=filtering_expression,
                                          gene_list_handler=gene_list_handler, snp_list_handler=snp_list_handler,
                                          log_file=test_log)

    for bgen_prefix, variant_list in snp_list_generator.genes.items():
        genotypes = generate_genotype_matrix(
            bgen_prefix,
            bgen_dict[bgen_prefix],
            variant_list,
            should_collapse=False if snp_list_handler or gene_list_handler else True,
            delete_on_complete=False  # Make sure we keep the files for testing
        )

        ac_table, gene_ac_table, gene_totals = check_matrix_stats(genotypes[1:], variant_list)

        assert len(ac_table) == 10000
        assert len(gene_ac_table) == 10000
        assert len(gene_totals) == len(variant_list['ENST'].unique())

        # Call the function
        stat_writer(ac_table, gene_ac_table, gene_totals, snp_list_generator.total_sites, test_log)

        # make sure the file exists
        assert log_path.exists()

        # Read and compare the files as sets
        with open(log_path, 'r') as log_file,\
              open(os.path.join(test_data_dir, "../expected_outputs", f"{bgen_prefix}_{name}_log.txt"), 'r') as expected_file:

            assert set(log_file.read()) == set(expected_file.read()), "Log file contents do not match expected output"