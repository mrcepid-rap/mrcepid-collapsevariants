from pathlib import Path
from typing import IO

from dxpy import DXFile
from general_utilities.association_resources import generate_linked_dx_file


class CollapseLOGGER:
    """A class to log information about filtering that was performed as part of this applet.

    This class consists of a File I/O and a number of methods to format various element types into text within this
    file. These generally come in the format (using pythong f-strings) of:

    1. text
    2. A spacer with length of <space 50 - len(text)>
    3. A colon ':'
    4. Some piece of information (int, str, etc.)

    :param file_prefix: A name to append to beginning of the output.
    """

    def __init__(self, log_path: Path):
        self._log_path = log_path
        self._LOG_FILE = self._open_writer()
        self._header_width = 30
        self._spacer = '-'
        self._line_width = 50
        self._vep_width = 130

    def _open_writer(self) -> IO:
        """Open the logfile filehandle"""

        log_file = self._log_path.open('w')
        return log_file

    def close(self):
        """Close the logfile filehandle"""
        self._LOG_FILE.close()

    def close_and_upload(self) -> DXFile:
        """Assert that the logfile has been closed AND upload it to the DNANexus platform"""

        self.close()
        linked_log_file = generate_linked_dx_file(self._log_path)
        return linked_log_file

    def write_header(self, text: str) -> None:
        """Write a header line with a title surrounded by a '-' spacer.

        This can only support header lines shorter than 50 characters in length. Otherwise it will truncate the output.

        :param text: The header text to write
        """
        write_string = f'{text:{self._spacer}^{self._header_width}.{self._header_width}}\n'
        self._write(write_string)

    def write_int(self, text: str, number: int, is_vep: bool = False) -> None:
        """Write an integer either with some whitespace or a large amount of whitespace if a VEP count

        :param text: The text label for this output
        :param number: The value to write
        :param is_vep: Do we need an extra bit of space to print long vep strings?
        """

        line_width = self._vep_width if is_vep else self._line_width
        write_string = f'{text:{line_width}.{line_width}}: {number}\n'
        self._write(write_string)

    def write_float(self, text: str, number: float) -> None:
        """Write a sprintf formatted float

        :param text: The text label for this output
        :param number: The value to write, formatted to three decimal places (`0.3f`)
        """

        write_string = f'{text:{self._line_width}.{self._line_width}}: {number:0.3f}\n'
        self._write(write_string)

    def write_scientific(self, text: str, number: float) -> None:
        """Write a sprintf formatted scientific float

        :param text: The text label for this output
        :param number: The value to write, formatted to three decimal places in scientific notation (`0.3e`)
        """

        write_string = f'{text:{self._line_width}.{self._line_width}}: {number:0.3e}\n'
        self._write(write_string)

    def write_histogram(self, bin: str, count: int) -> None:
        """Write a histogram bin

        :param bin: The bin label (x)
        :param count: The bin count (y)
        """

        write_string = f'{bin}\t{count}\n'
        self._write(write_string)

    def write_string(self, text: str, value: str) -> None:
        """Write a string with some formatted text

        :param text: The text label for this output
        :param value: The string to write
        """

        write_string = f'{text:{self._line_width}.{self._line_width}}: {value}\n'
        self._write(write_string)

    def write_generic(self, text: str) -> None:
        """Write a generic string to the log without formatting

        :param text: The string to write
        """

        write_string = f'{text}\n'
        self._write(write_string)

    def write_spacer(self):
        """Write a carriage return in the logfile"""

        write_string = '\n'
        self._write(write_string)

    def _write(self, write_string: str):
        """Actually write the formatted string to the log

        :param write_string: The formatted f-string in str format
        """

        self._LOG_FILE.write(write_string)
        self._LOG_FILE.flush()
