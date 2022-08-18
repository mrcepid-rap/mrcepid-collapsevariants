import dxpy

from typing import IO
from dxpy import DXFile


class CollapseLOGGER:

    def __init__(self, file_prefix: str):

        self._file_prefix = file_prefix
        self._LOG_FILE = self._open_writer()
        self._header_width = 30
        self._spacer = '-'
        self._line_width = 50
        self._vep_width = 130

    def _open_writer(self) -> IO:
        log_file = open(f'{self._file_prefix}.log', 'w')
        return log_file

    def close_writer(self) -> DXFile:

        self._LOG_FILE.close()
        linked_log_file = dxpy.upload_local_file(f'{self._file_prefix}.log')
        return linked_log_file

    def write_header(self, text: str) -> None:

        write_string = f'{text:{self._spacer}^{self._header_width}}\n'
        self._write(write_string)

    def write_int(self, text: str, number: int, is_vep: bool) -> None:

        write_string = f'{text:{self._vep_width if is_vep else self._line_width}}: {number}\n'
        self._write(write_string)

    def write_float(self, text: str, number: float) -> None:

        write_string = f'{text:{self._line_width}}: {number:0.2f}\n'
        self._write(write_string)

    def write_scientific(self, text: str, number: float) -> None:

        write_string = f'{text:{self._line_width}}: {number:0.3e}\n'
        self._write(write_string)

    def write_histogram(self, bin: int, count: int) -> None:

        write_string = f'{bin}\t{count}\n'
        self._write(write_string)

    def write_generic(self, text: str) -> None:

        write_string = f'{text}\n'
        self._write(write_string)

    def write_spacer(self):

        write_string = '\n'
        self._write(write_string)

    def _write(self, write_string: str):

        self._LOG_FILE.write(write_string)
