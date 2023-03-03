import csv
from pathlib import Path
from typing import Tuple, Dict

from general_utilities.association_resources import run_cmd


class STAARMergingException(Exception):
    """An exception to raise for any possible error in this class"""

    def __init__(self, message: str = "STAAR processing failed"):
        super().__init__(message)


class STAARParser:

    def __init__(self, file_prefix: str, chromosome: str):
        """A class to process files for input into STAAR

        STAAR requires me to make a "matrix" with rows of sample IDs and Columns of individual variants: rows can be
        pulled from the sample file; columns can be derived from the original SNP list file generated in
        SNPListGenerator. The actual "insides" (i.e. matrix[i][j]) are printed in parse_filters_BOLT (file name like:
        `{file_prefix}.{chromosome}.parsed.txt`). Doing it in sparse matrix (marketMatrix) format where I just give
        the sample ID (row) x var ID (col) only for sample x variant combinations with a genotype (e.g., each row of
        the file generated in parse_filters_BOLT).

        :param file_prefix: A name to append to beginning of output files.
        :param chromosome: The chromosome currently being processed. This must be the short form of the chromosome name
            (e.g., '1' not 'chr1').
        """

        self._file_prefix = file_prefix
        self._chromosome = chromosome

    def parse_filters_STAAR(self) -> None:
        """Generate input format files that can be provided to STAAR

        This proceeds in 5 steps:

        1. Make the sample dictionary: columns of sample rowNum
        2. Make the variant dictionary: columns of variantID chrom pos ENST colNum (also retained in the final output)
        3. Get the overall number of sample x variants (i.e, the total number of records in the sparseMatrix)
        4. Write the matrix in tsv. matrixMarket format
        5. Create the sparseMatrix in R (using buildSTAARmatrix.R) and save it in compressed RDS format.

        :return: None
        :raises STAARMergingException: STAAR fails sometimes for unknown reasons, so we catch all exceptions and raise
            an exception specific to STAAR so the parent process knows to retry
        """

        try:
            samples_dict_file, samples = self._make_samples_dict()
            variants_dict_file, variants = self._make_variants_dict()
            file_length = self._get_num_records()
            matrix_file = self._make_matrix_file(samples, variants, file_length)

            # And then run the R script `buildSTAARmatrix.R` to properly attach row and column names and save in .RDS
            # format for quick read/write in later parts of the pipeline. See the script itself for how it works,
            # but it takes 4 inputs:
            #
            # 1. Samples dict (samples_dict_file above) – the row names of our final matrix
            # 2. Variants dict (variants_dict_file above) – the column names of our final matrix
            # 3. Matrix file (matrix_file above) – the cells in our final matrix, formatted in %MatrixMarket format
            # 4. Out file – The name of the .rds file for final output
            #
            # And generates one output:

            # 1. A .rds file (named by the last argument) that can be read back into STAAR during
            # mrcepid-runassociationtesting

            cmd = f'Rscript /prog/buildSTAARmatrix.R ' \
                  f'/test/{samples_dict_file} ' \
                  f'/test/{variants_dict_file} ' \
                  f'/test/{matrix_file} ' \
                  f'/test/{self._file_prefix}.{self._chromosome}.STAAR.matrix.rds'
            run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting:latest',
                    docker_mounts=['/usr/bin/:/prog'])

            # Remove intermediate files so they aren't stored in the final tar:
            samples_dict_file.unlink()
            matrix_file.unlink()

        except Exception:
            raise STAARMergingException

    def _make_samples_dict(self) -> Tuple[Path, Dict[str, int]]:
        """Make a sample 'dictionary' in text format

        Make a text file of sample ID / rowNum for writing the final matrix. tsv

        :return: A Tuple containing a PathLike for the sample dictionary and a python dict() representation of this data
        """

        # First create a dictionary of sample row numbers
        # This is the raw .sample file that accompanies the bgen file
        samples = dict()
        samples_dict_file = Path(f'{self._file_prefix}.{self._chromosome}.samples.tsv')
        col_num = 1
        with samples_dict_file.open('w') as samples_dict_writer,\
                Path(f'{self._file_prefix}.{self._chromosome}.sample').open('r') as sample_file:
            samples_dict_writer.write('sampID\trow\n')
            sample_csv = csv.DictReader(sample_file, delimiter=' ', quoting=csv.QUOTE_NONE)
            for sample in sample_csv:
                if sample['ID'] != "0":  # This gets rid of the wierd header row in bgen sample files...
                    samples[sample['ID']] = col_num
                    samples_dict_writer.write(f'{sample["ID"]}\t{col_num}\n')
                    col_num += 1
            samples_dict_writer.close()

        return samples_dict_file, samples

    def _make_variants_dict(self) -> Tuple[Path, Dict[str, int]]:
        """Make a variant 'dictionary' in text format

        Makes a text file of variant ID / chrom / pos / ENST / column for use in making the final matrix .tsv but also
        as a dictionary of possible variants for downstream processing both by STAAR and other tools.

        :return: A Tuple containing a PathLike for the variant dictionary and a python dict() representation of this
            data
        """

        # ... Then a dictionary of variant column IDs
        # This file is also saved for runassociationtesting
        variants = dict()
        variants_dict_file = Path(f'{self._file_prefix}.{self._chromosome}.variants_table.STAAR.tsv')
        row_num = 1
        with variants_dict_file.open('w') as variants_dict_writer,\
                Path('snp_ENST.txt').open('r') as variants_list_file:

            variants_list_csv = csv.DictReader(variants_list_file, delimiter='\t', quoting=csv.QUOTE_NONE)
            variants_dict_csv = csv.DictWriter(variants_dict_writer, delimiter='\t', extrasaction='ignore',
                                               fieldnames=['varID', 'chrom', 'pos', 'ENST', 'column'])
            variants_dict_csv.writeheader()
            for var in variants_list_csv:
                if var['chrom'] == self._chromosome or self._chromosome == 'SNP' or self._chromosome == 'GENE':
                    variants[var['varID']] = row_num

                    var['column'] = row_num
                    variants_dict_csv.writerow(var)
                    row_num += 1

        return variants_dict_file, variants

    def _get_num_records(self) -> int:
        """ Get the number of non-homozygous reference sample x variant pairs

        Uses a simple system wc -l call to get the number of possible sample x variant pairs as required by the
        header of the matrixMarket format. wc is faster than reading the file in for very large numbers of records.

        :return: The total number of sample x variant pairs
        """

        # Need to get the file length of this file...
        wc_file = Path(f'{self._file_prefix}.{self._chromosome}_wc.txt')
        run_cmd(f'wc -l {self._file_prefix}.{self._chromosome}.parsed.txt', is_docker=False, stdout_file=wc_file)
        with wc_file.open('r') as wc_reader:
            for line in wc_reader:
                line = line.rstrip()
                file_length = int(line.split()[0])
        wc_file.unlink()

        return file_length

    def _make_matrix_file(self, samples: Dict[str, int], variants: Dict[str, int], file_length: int) -> Path:
        """Write a matrixMarket format file

        This method writes a .tsv representation of a sparse matrix in matrixMarket format containing only non-reference
        sample x genotype combinations. This matrix is the primary genotype representation used by STAAR for burden
        testing.

        :param samples: The sample dictionary from :func:`_make_samples_dict()`
        :param variants: The variants dictionary from :func:`_make_variants_dict()`
        :param file_length: The total number of sample x variant pairs as calculated by :func:`_get_num_records`
        :return: A PathLike pointing to the matrixMarket format file created by this method
        """
        row_num = max(samples.values())
        col_num = max(variants.values())

        # And format it to R Matrix (matrixMarket) spec. I can't use a DictWriter here because of the double header
        # lines which makes it not fit csv/tsv standard
        matrix_file = Path(f'{self._file_prefix}.{self._chromosome}.STAAR.matrix.R.txt')
        with Path(f'{self._file_prefix}.{self._chromosome}.parsed.txt').open('r') as sparse_matrix,\
                matrix_file.open('w') as matrix_file_writer:

            sparse_matrix = csv.DictReader(sparse_matrix, delimiter='\t', quoting=csv.QUOTE_NONE,
                                           fieldnames=['sample', 'varID', 'genotype'])

            # Write the header in %MatrixMarket format
            matrix_file_writer.write('%%MatrixMarket matrix coordinate integer general\n')
            matrix_file_writer.write(f'{row_num} {col_num} {file_length}\n')
            # Write the individual cells in the matrix
            for row in sparse_matrix:
                gt_val = 1 if row['genotype'] == '0/1' else 2
                matrix_file_writer.write(f'{samples[row["sample"]]} {variants[row["varID"]]} {gt_val}\n')
            matrix_file_writer.close()

        return matrix_file
