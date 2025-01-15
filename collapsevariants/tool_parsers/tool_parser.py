from abc import abstractmethod, ABC
from pathlib import Path
from typing import Dict, List

import pandas as pd

from scipy.sparse import csr_matrix

from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.mrc_logger import MRCLogger


class ToolParser(ABC):
    """An Interface for creating files for burden testing.

    This Interface implements the general functionality for creating files for burden testing for all possible
    tools that are implemented downstream in this pipeline. Generally speaking it:

    1. Requests input of variant lists (:param genes:) and genotype indices (:param genotype_index:) for each BGEN file.
    Note that the genotype indices are not necessarily used by every implementing class.

    2. It sets an output list for holding various outputs of implementing classes.

    3. It uses ThreadUtility to parallelize the creation of files for burden testing using the abstract method
    :func:`_make_output_files`. It is required by implementing classes to define this method to create files that
    are specific to the tool they are implementing.

    This means that the implementing classes DO NOT need to worry about parallelizing the creation of files. This is
    automatically handled by the ToolParser class automatically parallelizing the :func:`_make_output_files` abstract
    method.

    :param genes: A dictionary containing the genes to collapse with keys of the BGEN file prefixes and values of
    :param genotype_index: A dictionary containing values of csr_matrix and keys of each BGEN file prefix.
    :param sample_ids: A list of sample IDs for processing this file.
    :param output_prefix: A string representing the prefix of the output files.
    :param tool_name: A string representing the name of the tool of the implementing class.
    """

    def __init__(self, genes: Dict[str, pd.DataFrame], genotype_index: Dict[str, csr_matrix], sample_ids: List[str],
                 output_prefix: str, tool_name: str):

        # Initialize the logger
        self._logger = MRCLogger(__name__).get_logger()

        # Dereference the input variables
        self._genes = genes
        self._genotype_index = genotype_index
        self._sample_ids = sample_ids
        self._output_prefix = output_prefix
        self._tool_name = tool_name

        # Initialize the output files list
        self._output_files = []

        # Call the multithreaded file creation method
        self._multithread_file_creation()

    def get_output_files(self) -> List[Path]:
        """Getter for the file-based output of this interface / implementing classes

        :return: A list of Path objects representing the output files created by the implementing class
        """

        return self._output_files

    def _multithread_file_creation(self):
        """Parallelize the creation of files for burden testing.

        This method uses the ThreadUtility class to parallelize the abstract method :func:`_make_output_files`. This
        method is automatically called by the constructor of this class and should not be called directly by
        implementing classes.

        :return: None
        """

        self._logger.info(f'Generating {self._tool_name} files for burden testing')

        # Use the make_output_files helper class to parallelize the creating of tool-specific outputs.
        # This is mostly for readability and to avoid having to write a lot of boilerplate code to manage threads
        thread_utility = ThreadUtility(error_message=f'Error in {self._tool_name} processing step')

        # Note that iteration MUST be keyed on self._genes!!! This is because genes tracks which input bgen files
        # actually had found variants based on requested filtering.
        for bgen_prefix in self._genes.keys():
            thread_utility.launch_job(self._make_output_files,
                                      bgen_prefix=bgen_prefix)

        for result in thread_utility:
            self._output_files.extend(result)

    @abstractmethod
    def _make_output_files(self, bgen_prefix: str) -> List[Path]:
        """Abstract method to create files for burden testing.

        This method is required to be implemented by all classes that inherit from this interface. It is the run
        automatically by the constructor of the ToolParser interface.

        :param bgen_prefix: A string representing the prefix of a BGEN file to convert to the tool-specific format.
        :return: A List of Path objects representing the output files created by the implementing class
        """
        pass