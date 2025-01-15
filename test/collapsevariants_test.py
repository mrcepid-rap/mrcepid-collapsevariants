import numpy as np
import pytest
import filecmp
import pandas as pd

from pathlib import Path
from scipy.sparse import csr_matrix

from collapsevariants.collapse_logger import CollapseLOGGER
from collapsevariants.collapse_utils import generate_csr_matrix_from_bgen, get_sample_ids, check_matrix_stats
from collapsevariants.snp_list_generator import SNPListGenerator


def generate_expected_counts(test_data: Path, gene_list: Path = None, snp_list: Path = None) -> np.ndarray:
    gts = pd.read_csv(test_data,
        sep='\t', names=['CHROM', 'POS', 'ID', 'ENST', 'CSQ', 'LOFTEE', 'MAF', 'SAMPLE', 'GT'])

    gts['GT'] = gts['GT'].apply(lambda x: 1 if x == '0/1' else 2)
    gts['SAMPLE'] = gts['SAMPLE'].str.split('_', expand=True)[1].astype(int)

    # Filter by ENST if required
    if gene_list:
        with gene_list.open('r') as gene_file:
            gene_list = [gene.rstrip() for gene in gene_file.readlines()]

        gts = gts[gts['ENST'].isin(gene_list)]

    # Filter by SNP if required
    if snp_list:
        with snp_list.open('r') as snp_file:
            snp_list = [snp.rstrip() for snp in snp_file.readlines()]

        gts = gts[gts['ID'].isin(snp_list)]

    gt_totals = gts.groupby('SAMPLE').aggregate(GT_sum=('GT', 'sum'))

    # Make a dummy frame to make sure we capture all sample IDs
    dummy = pd.DataFrame(index=range(0, 10000))
    dummy.index.name = 'SAMPLE'

    gt_totals = dummy.merge(gt_totals, left_index=True, right_index=True, how='left')
    gt_totals = gt_totals.fillna(0)

    return gt_totals['GT_sum'].to_numpy()


# Validated test data:
test_dir = Path(__file__).parent
test_data_dir = test_dir / 'test_data/'
correct_log = test_data_dir / 'correct_log.txt'

# Filtering lists
snp_path = test_data_dir / 'snp_list.v2.txt'
gene_enst_path = test_data_dir / 'gene_list.ENST.txt'
gene_symbol_path = test_data_dir / 'gene_list.SYMBOL.txt'

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

def test_logger(tmp_path):
    """Test the logging functionality coded for tracking variants

    This function takes no inputs but uses a correctly formatted log from a previous run of this code
    (correct_log.txt) to determine accuracy of these tests. This function also tests upload and download from the
    DNANexus file system.
    """

    log_path = tmp_path / 'test_log.log'
    assert log_path.exists() is False
    LOG_FILE = CollapseLOGGER(log_path)
    assert log_path.exists()

    LOG_FILE.write_header('THIS IS A TEST')
    LOG_FILE.write_int('test int no vep', 1, False)
    LOG_FILE.write_int('test int vep', 1, True)
    LOG_FILE.write_float('test float', 0.0001)  # should be 0.000
    LOG_FILE.write_float('test float', 0.12369)  # Should be 0.124
    LOG_FILE.write_scientific('test sci', 0.0001)  # should be 1.000e-04
    LOG_FILE.write_scientific('test sci', 0.12369)  # should be 1.237e-01
    LOG_FILE.write_string('test str', 'foo')
    LOG_FILE.write_generic('test generic')
    LOG_FILE.write_spacer()
    LOG_FILE.write_histogram(0, 1)
    LOG_FILE.write_histogram(1, 5)
    LOG_FILE.write_histogram(2, 10)
    LOG_FILE.write_histogram(3, 5)
    LOG_FILE.write_histogram(4, 1)
    LOG_FILE.write_spacer()
    LOG_FILE.write_int('This text is too long to be included and should be truncated', 1)
    LOG_FILE.write_header('This text is too long to be included and should be truncated')
    LOG_FILE.close()

    # And make sure the contents match the known data exactly
    assert filecmp.cmp(correct_log, log_path)


@pytest.mark.parametrize("filtering_expression, gene_list_path, snp_list_path, expected_num_sites, check_error", [
    ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', None, None, 13592, False),
    ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', gene_enst_path, None, 34, False),
    ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', gene_symbol_path, None, 34, False),
    (None, None, snp_path, 826, False),
    (None, None, None, 0, True),
    (None, gene_enst_path, None, 0, True),
    ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', gene_enst_path, snp_path, 0, True),
    ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', None, snp_path, 0, True)])
def test_snp_list_generator(tmp_path: Path, filtering_expression: str, gene_list_path: Path, snp_list_path: Path,
                            expected_num_sites: int, check_error: bool):

    log_path = tmp_path / 'HC_PTV-MAF_01.log'
    test_log = CollapseLOGGER(log_path)

    # Has to be inside the test function since I have to re-open every time the test is run
    # Also make sure this is rb, as we wrap gzip around this in the function.
    vep_dict = {bgen_prefix: bgen_info['vep'].open('rb') for bgen_prefix, bgen_info in bgen_dict.items()}

    # Testing to make sure filtering mode parsing works correctly
    if check_error:
        with pytest.raises(ValueError):
            SNPListGenerator(vep_dict=vep_dict, filtering_expression=filtering_expression,
                             gene_list_path=gene_list_path, snp_list_path=snp_list_path, log_file=test_log)
    else:
        snp_list_generator = SNPListGenerator(vep_dict=vep_dict, filtering_expression=filtering_expression,
                                              gene_list_path=gene_list_path, snp_list_path=snp_list_path, log_file=test_log)

        assert snp_list_generator.total_sites == expected_num_sites

        table_sites = 0
        table_vars = set()
        table_genes = set()
        for bgen_prefix, variant_index in snp_list_generator.genes.items():
            table_sites += variant_index.shape[0]
            table_vars.update(variant_index['varID'].tolist())
            table_genes.update(variant_index['ENST'].tolist())

            assert bgen_prefix in vep_dict.keys()
            assert variant_index.columns.tolist() == ['varID', 'CHROM', 'POS', 'ENST']

        # Make sure we find the correct number of snps / genes if using that mode
        if snp_list_path:
            with snp_list_path.open('r') as snp_file:
                expected_vars = set([var.rstrip() for var in snp_file.readlines()])

            assert len(expected_vars.union(table_vars)) == expected_num_sites

        elif gene_list_path:
            with gene_list_path.open('r') as gene_file:
                expected_genes = set([gene.rstrip() for gene in gene_file.readlines()])

            assert len(table_genes) == len(expected_genes)

        assert table_sites == expected_num_sites

    test_log.close()


@pytest.mark.parametrize("sample_file", [bgen_info['sample'] for bgen_info in bgen_dict.values()])
def test_get_sample_ids(sample_file: Path):
    sample_ids = get_sample_ids(sample_file)
    assert len(sample_ids) == 10000
    assert type(sample_ids) is list
    assert type(sample_ids[0]) is str
    assert len(sample_ids[0].split()) == 1   # Make sure the sample IDs are non-delimited strings


@pytest.mark.parametrize("filtering_expression, gene_list_path, snp_list_path, expected_num_sites",
                         [('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', None, None, 13592),
                          ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', gene_enst_path, None, 34),
                          (None, None, snp_path, 826)])
def test_csr_matrix_generation(tmp_path: Path, filtering_expression: str, gene_list_path: Path,
                               snp_list_path: Path, expected_num_sites: int):
    """This tests both :func:`generate_csr_matrix_from_bgen` and :func:`check_matrix_stats` functions.

    :param tmp_path:
    :param filtering_expression:
    :param gene_list_path:
    :param snp_list_path:
    :param expected_num_sites:
    :return:
    """

    log_path = tmp_path / 'HC_PTV-MAF_01.log'
    test_log = CollapseLOGGER(log_path)

    # Has to be inside the test function since I have to re-open every time the test is run
    # Also make sure this is rb, as we wrap gzip around this in the function.
    vep_dict = {bgen_prefix: bgen_info['vep'].open('rb') for bgen_prefix, bgen_info in bgen_dict.items()}

    snp_list_generator = SNPListGenerator(vep_dict=vep_dict, filtering_expression=filtering_expression,
                                          gene_list_path=gene_list_path, snp_list_path=snp_list_path, log_file=test_log)

    total_variants = 0
    for bgen_prefix, variant_list in snp_list_generator.genes.items():

        genotypes = generate_csr_matrix_from_bgen(variant_list,
                                                  bgen_dict[bgen_prefix]['bgen'],
                                                  bgen_dict[bgen_prefix]['sample'])

        ac_table, gene_ac_table, gene_totals = check_matrix_stats(genotypes, variant_list)
        assert len(ac_table) == 10000
        assert len(gene_ac_table) == 10000
        assert len(gene_totals) == len(variant_list['ENST'].unique())

        expected_sites_path = bgen_dict[bgen_prefix]['gts']
        expected_counts = generate_expected_counts(expected_sites_path,
                                                   gene_list=gene_list_path, snp_list=snp_list_path)

        # EUGENE – Remember that this fails due to an annotation issue in Duat
        # The problem is:
        # 1. There are two variants with identical pos / ref / alt alleles but different CSQ annotations
        # 2. bcftools annotate cannot tell the difference and just uses the 1st annotation
        # I have left in this numpy bit so that you can see where the error is. It prints out the offending
        # sample where compiled gt != expected gt;
        # tldr: the test data is wrong, not collapsevariants
        print(bgen_prefix)
        print(np.argwhere(np.not_equal(ac_table, expected_counts)))
        assert np.array_equal(ac_table, expected_counts)

        total_variants += genotypes.shape[1]

        assert type(genotypes) is csr_matrix

    assert total_variants == expected_num_sites
