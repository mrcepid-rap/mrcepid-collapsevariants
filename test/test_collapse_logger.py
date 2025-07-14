import filecmp
from pathlib import Path

from utilities.collapse_logger import CollapseLOGGER

# Validated test data:
test_dir = Path(__file__).parent
test_data_dir = test_dir / 'test_data/'
correct_log = test_data_dir / 'correct_log.txt'


def test_logger(tmp_path):
    """
    Test the logging functionality coded for tracking variants

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
