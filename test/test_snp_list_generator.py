# test_stat_dictionary.py
import pytest
import filecmp
from pathlib import Path
import pandas as pd
from collapsevariants.collapse_logger import CollapseLOGGER
from collapsevariants.snp_list_generator import StatDictionary

# Set the test data parameters:
test_dir = Path(__file__).parent
# test_dir=Path("/Users/alish.palmos/PycharmProjects/mrcepid-collapsevariants/test")
test_data_dir = test_dir / 'test_data/'
correct_log = test_data_dir / 'correct_log.txt'


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


# data = pd.read_csv(gzip.open("/Users/alish.palmos/PycharmProjects/mrcepid-collapsevariants/test/test_data/chr1_chunk1.vep.tsv.gz"
#                              , mode='rt'), sep="\t",
#                                   index_col='varID',
#                                   dtype={'SIFT': str, 'POLYPHEN': str, 'LOFTEE': str,
#                                          'AA': str, 'AApos': str})


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
    assert stat_dict._max_cadd == 0.456
    assert stat_dict._max_revel == 0.789
    assert stat_dict._num_na_revel == 14063
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
            "Minimum CADD Score",
            "Minimum REVEL Score",
            "Number of NA REVEL Scores",
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
            "0.456",
            "0.789",
            "14063",
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
    stat_dict._update_max_cadd(data2['CADD'])
    stat_dict._update_revel(data2['REVEL'])
    stat_dict._update_pass_stats(data2['FILTER'])
    stat_dict._update_loftee_counts(data2['LOFTEE'])
    stat_dict._update_parsed_consequence_counts(data2['PARSED_CSQ'])
    stat_dict._update_vep_consequence_counts(data2['CSQ'])

    # check the results are as expected
    # Assert the statistics of the collect_stats function
    assert stat_dict._total_variants == 36403
    assert stat_dict._max_missingness == 12.3
    assert stat_dict._max_af == 0.001
    assert stat_dict._max_cadd == 0.456
    assert stat_dict._max_revel == 0.789
    assert stat_dict._num_na_revel == 24366
    assert stat_dict._num_pass == 36398
    assert stat_dict._num_not_pass == 0
    assert stat_dict._loftee_counts == {'HC': 11939}
    assert stat_dict._parsed_consequence_counts == {'SYN': 12427, 'PTV': 11939, 'MISSENSE': 12032}
    assert stat_dict._vep_consequence_counts == {'SYN&testcsq1&testcsq2': 6208, 'SYN&testcsq1': 6219,
                                                 'PTV&testcsq1&testcsq2': 6027, 'MISS&testcsq1': 6058,
                                                 'PTV&testcsq1': 5912, 'MISS&testcsq1&testcsq2': 5974}
