import os
from pathlib import Path

import pytest

from collapsevariants.collapse_logger import CollapseLOGGER
from collapsevariants.collapse_utils import generate_csr_matrix_from_bgen, get_sample_ids, check_matrix_stats
from collapsevariants.parallelization_wrappers import stat_writer
from collapsevariants.snp_list_generator import SNPListGenerator

# Validated test data:
test_dir = Path(__file__).parent
test_data_dir = test_dir / 'test_data/'

# Filtering lists
snp_path = test_data_dir / 'snp_list.v2.txt'
gene_enst_path = test_data_dir / 'gene_list.ENST.txt'

# Variant information
bgen_dict = {'chr1_chunk1': {'index': test_data_dir / 'chr1_chunk1.bgen.bgi',
                             'bgen': test_data_dir / 'chr1_chunk1.bgen',
                             'sample': test_data_dir / 'chr1_chunk1.sample',
                             'vep': test_data_dir / 'chr1_chunk1.vep.tsv.gz',
                             'gts': test_data_dir / 'chr1_chunk1.gts'},
             'chr1_chunk2': {'index': test_data_dir / 'chr1_chunk2.bgen.bgi',
                             'bgen': test_data_dir / 'chr1_chunk2.bgen',
                             'sample': test_data_dir / 'chr1_chunk2.sample',
                             'vep': test_data_dir / 'chr1_chunk2.vep.tsv.gz',
                             'gts': test_data_dir / 'chr1_chunk2.gts'},
             'chr1_chunk3': {'index': test_data_dir / 'chr1_chunk3.bgen.bgi',
                             'bgen': test_data_dir / 'chr1_chunk3.bgen',
                             'sample': test_data_dir / 'chr1_chunk3.sample',
                             'vep': test_data_dir / 'chr1_chunk3.vep.tsv.gz',
                             'gts': test_data_dir / 'chr1_chunk3.gts'}}


@pytest.mark.parametrize("sample_file", [bgen_info['sample'] for bgen_info in bgen_dict.values()])
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


@pytest.mark.parametrize(
    "filtering_expression, gene_list_path, snp_list_path",
    [
        ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', None, None),
        ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', gene_enst_path, None),
        (None, None, snp_path)
    ]
)
def test_stat_writer(tmp_path: Path, filtering_expression: str, gene_list_path: Path,
                     snp_list_path: Path):
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

    # Has to be inside the test function since I have to re-open every time the test is run
    # Also make sure this is rb, as we wrap gzip around this in the function.
    vep_dict = {bgen_prefix: bgen_info['vep'].open('rb') for bgen_prefix, bgen_info in bgen_dict.items()}

    snp_list_generator = SNPListGenerator(vep_dict=vep_dict, filtering_expression=filtering_expression,
                                          gene_list_path=gene_list_path, snp_list_path=snp_list_path, log_file=test_log)

    for bgen_prefix, variant_list in snp_list_generator.genes.items():
        genotypes = generate_csr_matrix_from_bgen(variant_list,
                                                  bgen_dict[bgen_prefix]['bgen'],
                                                  bgen_dict[bgen_prefix]['sample'])

        ac_table, gene_ac_table, gene_totals = check_matrix_stats(genotypes, variant_list)

        # Call the function
        stat_writer(ac_table, gene_ac_table, gene_totals, snp_list_generator.total_sites, test_log)

        # make sure the file exists
        assert os.path.exists(log_path)

        # read in the log file
        with open(log_path, 'r') as log_file:
            new_log = log_file.read()

        # also read in the correct log file & make sure they match
        with open(os.path.join(test_data_dir, "../expected_outputs", f"{bgen_prefix}_log.txt"), 'r') as log_file:
            log_contents = log_file.read()
            assert set(log_contents) == set(new_log)
