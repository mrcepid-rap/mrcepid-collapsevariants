# test_stat_dictionary.py
import pytest
import filecmp
from pathlib import Path
import pandas as pd
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler

from utilities.collapse_logger import CollapseLOGGER
from collapsevariants.snp_list_generator.snp_list_generator import SNPListGenerator
from collapsevariants.snp_list_generator.stat_dictionary import StatDictionary

# Test data set-up
test_dir = Path(__file__).parent
test_data_dir = test_dir / 'test_data/'
correct_log = test_data_dir / 'correct_log.txt'

# Filtering lists
snp_path = InputFileHandler(test_data_dir / 'snp_list.v2.txt')
gene_enst_path = InputFileHandler(test_data_dir / 'gene_list.ENST.txt')
gene_symbol_path = InputFileHandler(test_data_dir / 'gene_list.SYMBOL.txt')


# Variant information
bgen_dict = {'chr1_chunk1': {'index': InputFileHandler(test_data_dir / 'chr1_chunk1.bgen.bgi'),
                             'bgen': InputFileHandler(test_data_dir / 'chr1_chunk1.bgen'),
                             'sample': InputFileHandler(test_data_dir / 'chr1_chunk1.sample'),
                             'vep': InputFileHandler(test_data_dir / 'chr1_chunk1.vep.tsv.gz'),
                             'gts': InputFileHandler(test_data_dir / 'chr1_chunk1.gts')},
             'chr1_chunk2': {'index': InputFileHandler(test_data_dir / 'chr1_chunk2.bgen.bgi'),
                             'bgen': InputFileHandler(test_data_dir / 'chr1_chunk2.bgen'),
                             'sample': InputFileHandler(test_data_dir / 'chr1_chunk2.sample'),
                             'vep': InputFileHandler(test_data_dir / 'chr1_chunk2.vep.tsv.gz'),
                             'gts': InputFileHandler(test_data_dir / 'chr1_chunk2.gts')},
             'chr1_chunk3': {'index': InputFileHandler(test_data_dir / 'chr1_chunk3.bgen.bgi'),
                             'bgen': InputFileHandler(test_data_dir / 'chr1_chunk3.bgen'),
                             'sample': InputFileHandler(test_data_dir / 'chr1_chunk3.sample'),
                             'vep': InputFileHandler(test_data_dir / 'chr1_chunk3.vep.tsv.gz'),
                             'gts': InputFileHandler(test_data_dir / 'chr1_chunk3.gts')}}


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


def parse_log_file_to_df(log_file_path):
    """
    Parse the log file into a DataFrame
    This is easier for testing purposes
    """
    with log_file_path.open('r') as log_file:
        lines = log_file.readlines()

    data = {
        "Statistic": [],
        "Value": []
    }

    for line in lines:
        if ':' in line:
            parts = line.split(':')
            statistic = parts[0].strip()
            value = parts[1].strip()
            data["Statistic"].append(statistic)
            data["Value"].append(value)

    df = pd.DataFrame(data)
    return df


def test_stat_dictionary(tmp_path,
                         data1=pd.read_csv(test_data_dir / 'chr1_chunk1.vep.tsv.gz', sep="\t", ),
                         data2=pd.read_csv(test_data_dir / 'chr1_chunk2.vep.tsv.gz', sep="\t", )
                         ):
    """
    Test the StatDictionary class that collects statistics from the data and writes them to a log file.
    """

    # Create a CollapseLOGGER instance
    log_path = tmp_path / 'HC_PTV-MAF_01.log'
    test_log = CollapseLOGGER(log_path)

    # Create a StatDictionary instance
    stat_dict = StatDictionary(test_log)
    stat_dict.collect_stats(data1)

    # collect_stats function
    # Assert the statistics of the collect_stats function
    assert stat_dict._total_variants == 20897
    assert stat_dict._max_missingness == 12.3
    assert stat_dict._max_af == 0.001
    assert stat_dict._num_pass == 20897
    assert stat_dict._num_not_pass == 0
    assert stat_dict._loftee_counts == {'HC': 6930}
    assert stat_dict._parsed_consequence_counts == {'SYN': 7133, 'PTV': 6930, 'MISSENSE': 6834}
    assert stat_dict._vep_consequence_counts == {'SYN&testcsq1&testcsq2': 3612, 'SYN&testcsq1': 3521,
                                                 'PTV&testcsq1&testcsq2': 3514, 'MISS&testcsq1': 3429,
                                                 'PTV&testcsq1': 3416, 'MISS&testcsq1&testcsq2': 3405}

    # write_bgen_stats function
    stat_dict.write_bgen_stats()
    # make sure the output log file exists
    assert log_path.exists()
    # Parse the log file content into a DataFrame
    df = parse_log_file_to_df(log_path)
    # Expected data
    expected_data = {
        "Statistic": [
            "Total number of variants",
            "Maximum missingness",
            "Maximum Allele Frequency",
            "Total number of PASS variants",
            "Total number of non-PASS variants",
            "Number of LOFTEE HC",
            "Number of Parsed Consequence – SYN",
            "Number of Parsed Consequence – PTV",
            "Number of Parsed Consequence – MISSENSE",
            "Number of VEP Consequence - SYN&testcsq1&testcsq2",
            "Number of VEP Consequence - SYN&testcsq1",
            "Number of VEP Consequence - PTV&testcsq1&testcsq2",
            "Number of VEP Consequence - MISS&testcsq1",
            "Number of VEP Consequence - PTV&testcsq1",
            "Number of VEP Consequence - MISS&testcsq1&testcsq2"
        ],
        "Value": [
            "20897",
            "12.300",
            "1.000e-03",
            "20897",
            "0",
            "6930",
            "7133",
            "6930",
            "6834",
            "3612",
            "3521",
            "3514",
            "3429",
            "3416",
            "3405"
        ]
    }
    # Create a DataFrame from the expected data
    expected_df = pd.DataFrame(expected_data)
    # Assert the DataFrame content
    pd.testing.assert_frame_equal(df, expected_df)

    # get_total_sites function
    # Assert the total sites
    assert stat_dict.get_total_sites() == 20897
    # Increment the variant count by 5
    stat_dict._increment_variant_count(5)
    # Assert the total sites
    assert stat_dict.get_total_sites() == 20902

    # test the functions which update the statistics
    stat_dict._increment_variant_count(len(data2))
    stat_dict._update_max_missingness(data2['F_MISSING'])
    stat_dict._update_max_af(data2['AF'])
    stat_dict._update_pass_stats(data2['FILTER'])
    stat_dict._update_loftee_counts(data2['LOFTEE'])
    stat_dict._update_parsed_consequence_counts(data2['PARSED_CSQ'])
    stat_dict._update_vep_consequence_counts(data2['CSQ'])

    # check the results are as expected
    # Assert the statistics of the collect_stats function
    assert stat_dict._total_variants == 36403
    assert stat_dict._max_missingness == 12.3
    assert stat_dict._max_af == 0.001
    assert stat_dict._num_pass == 36398
    assert stat_dict._num_not_pass == 0
    assert stat_dict._loftee_counts == {'HC': 11939}
    assert stat_dict._parsed_consequence_counts == {'SYN': 12427, 'PTV': 11939, 'MISSENSE': 12032}
    assert stat_dict._vep_consequence_counts == {'SYN&testcsq1&testcsq2': 6208, 'SYN&testcsq1': 6219,
                                                 'PTV&testcsq1&testcsq2': 6027, 'MISS&testcsq1': 6058,
                                                 'PTV&testcsq1': 5912, 'MISS&testcsq1&testcsq2': 5974}


@pytest.mark.parametrize("filtering_expression, gene_list_handler, snp_list_handler, expected_num_sites, check_error", [
    ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', None, None, 13592, False),
    ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', gene_enst_path, None, 34, False),
    ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', gene_symbol_path, None, 34, False),
    (None, None, snp_path, 826, False),
    (None, None, None, 0, True),
    (None, gene_enst_path, None, 0, True),
    ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', gene_enst_path, snp_path, 0, True),
    ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', None, snp_path, 0, True)])
def test_snp_list_generator(tmp_path: Path, filtering_expression: str,
                            gene_list_handler: InputFileHandler, snp_list_handler,
                            expected_num_sites: int, check_error: bool):
    log_path = tmp_path / 'HC_PTV-MAF_01.log'
    test_log = CollapseLOGGER(log_path)

    # Testing to make sure filtering mode parsing works correctly
    if check_error:
        with pytest.raises(ValueError):
            SNPListGenerator(bgen_dict=bgen_dict, filtering_expression=filtering_expression,
                             gene_list_handler=gene_list_handler, snp_list_handler=snp_list_handler, log_file=test_log)

    else:
        snp_list_generator = SNPListGenerator(bgen_dict=bgen_dict, filtering_expression=filtering_expression,
                                              gene_list_handler=gene_list_handler, snp_list_handler=snp_list_handler,
                                              log_file=test_log)

        assert snp_list_generator.total_sites == expected_num_sites

        table_sites = 0
        table_vars = set()
        table_genes = set()
        for bgen_prefix, variant_index in snp_list_generator.genes.items():
            table_sites += variant_index.shape[0]
            table_vars.update(variant_index['varID'].tolist())
            table_genes.update(variant_index['ENST'].tolist())

            assert variant_index.columns.tolist() == ['varID', 'CHROM', 'POS', 'ENST']

        # Make sure we find the correct number of snps / genes if using that mode
        if snp_list_handler:
            with snp_list_handler.get_file_handle().open('r') as snp_file:
                expected_vars = set([var.rstrip() for var in snp_file.readlines()])

            assert len(expected_vars.union(table_vars)) == expected_num_sites

        elif gene_list_handler:
            with gene_list_handler.get_file_handle().open('r') as gene_file:
                expected_genes = set([gene.rstrip() for gene in gene_file.readlines()])

            assert len(table_genes) == len(expected_genes)

        assert table_sites == expected_num_sites

    test_log.close()
