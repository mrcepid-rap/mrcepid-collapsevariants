import csv
import subprocess
from pathlib import Path

from general_utilities.association_resources import run_cmd


class STAARMergingException(Exception):

    def __init__(self, message: str = "Merging failed"):
        super().__init__(message)


class STAARParser:

    def __init__(self, file_prefix: str, chromosome: str):

        # STAAR merging will occasionally fail for no obvious reason, this try-catch should allow for the calling method
        # to know this and then properly re-run the code.
        try:
            self._parse_filters_STAAR(file_prefix, chromosome)
        except Exception:
            raise STAARMergingException(f'STAAR chr {chromosome} failed to merge.')

    # Generate input format files that can be provided to STAAR
    @staticmethod
    def _parse_filters_STAAR(file_prefix: str, chromosome: str) -> None:

        # STAAR requires me to make a "matrix" with rows of sample IDs and Columns of individual variants:
        # rows can be pulled from the sample file from the raw .bgen
        # columns can be pulled from the raw .bgen file
        # The actual "insides" (i.e. matrix[i][j]) are printed in parse_filters_BOLT (file name like:
        # file_prefix + "." + chromosome + ".parsed.txt). Doing it in sparse matrix format where I just give the
        # sample ID (row) x var ID (col). So each row of the file generated in parse_filters_BOLT

        # First create a dictionary of sample row numbers
        # This is the raw .sample file that accompanies the bgen file
        samples = {}
        samples_dict_file = Path(f'{file_prefix}.{chromosome}.samples.tsv')
        col_num = 1
        with samples_dict_file.open('w') as samples_dict_writer,\
                Path(f'{file_prefix}.{chromosome}.sample').open('r') as sample_file:
            samples_dict_writer.write('sampID\trow\n')
            sample_csv = csv.DictReader(sample_file, delimiter=' ', quoting=csv.QUOTE_NONE)
            for sample in sample_csv:
                if sample['ID'] != "0":  # This gets rid of the wierd header row in bgen sample files...
                    samples[sample['ID']] = col_num
                    samples_dict_writer.write(f'{sample["ID"]}\t{col_num}\n')
                    col_num += 1
            samples_dict_writer.close()

        # ... Then a dictionary of variant column IDs
        # This file is also saved for runassociationtesting
        variants = {}
        variants_dict_file = Path(f'{file_prefix}.{chromosome}.variants_table.STAAR.tsv')
        row_num = 1
        with variants_dict_file.open('w') as variants_dict_writer,\
                Path('snp_ENST.txt').open('r') as variants_list_file:

            variants_list_csv = csv.DictReader(variants_list_file, delimiter="\t", quoting=csv.QUOTE_NONE)
            variants_dict_csv = csv.DictWriter(variants_dict_writer, delimiter='\t', extrasaction='ignore',
                                               fieldnames=['varID', 'chrom', 'pos', 'ENST', 'column'])
            variants_dict_csv.writeheader()
            for var in variants_list_csv:
                if var['CHROM'] == chromosome or chromosome == 'SNP' or chromosome == 'GENE':
                    variants[var['varID']] = row_num
                    var['column'] = row_num
                    variants_dict_csv.writerow(var)
                    row_num += 1

        # Need to get the file length of this file...
        # This code is janky AF
        wc_file = Path(f'{file_prefix}.{chromosome}_wc.txt')
        run_cmd(f'wc -l {file_prefix}.{chromosome}.parsed.txt', is_docker=False, stdout_file=wc_file)
        with wc_file.open('r') as wc_reader:
            for line in wc_reader:
                line = line.rstrip()
                file_length = int(line.split()[0])
        wc_file.unlink()

        # And format it to R Matrix (matrixMarket) spec. I can't use a DictWriter here because of the double header
        # lines which makes it not fit csv/tsv standard
        matrix_file = Path(f'{file_prefix}.{chromosome}.STAAR.matrix.R.txt')
        with Path(f'{file_prefix}.{chromosome}.parsed.txt').open('r') as sparse_matrix,\
                matrix_file.open('w') as matrix_file_writer:

            sparse_matrix = csv.DictReader(sparse_matrix, delimiter='\t', quoting=csv.QUOTE_NONE,
                                           fieldnames=['sample', 'varID', 'genotype'])

            # Write the header in %MatrixMarket format
            matrix_file_writer.write('%%MatrixMarket matrix coordinate integer general\n')
            matrix_file_writer.write(f'{col_num - 1} {row_num - 1} {file_length}\n')  # -1 because of 0-based header
            # Write the individual cells in the matrix
            for row in sparse_matrix:
                gt_val = 1 if row['genotype'] == '0/1' else 2
                matrix_file_writer.write(f'{samples[row["sample"]]} {variants[row["varID"]]} {gt_val}\n')
            matrix_file_writer.close()

        # And then run the R script `buildSTAARmatrix.R` to properly attach row and column names and save in .RDS format
        # for quick read/write in later parts of the pipeline.
        # See the script itself for how it works, but it takes 4 inputs:
        # 1. Samples dict (samples_dict_file above) – the row names of our final matrix
        # 2. Variants dict (variants_dict_file above) – the column names of our final matrix
        # 3. Matrix file (matrix_file above) – the cells in our final matrix, formatted in %MatrixMarket format for ease
        # 4. Out file – The name of the .rds file for final output
        #
        # And generates one output:
        # 1. A .rds file (named by the last argument) that can be read back into STAAR during
        #   mrcepid-runassociationtesting
        cmd = f'Rscript /prog/buildSTAARmatrix.R /test/{samples_dict_file} ' \
              f'/test/{variants_dict_file} ' \
              f'/test/{matrix_file} ' \
              f'/test/{file_prefix}.{chromosome}.STAAR.matrix.rds'
        run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting:latest',
                docker_mounts=['/usr/bin/:/prog'])

        # Remove intermediate files so they aren't stored in the final tar:
        samples_dict_file.unlink()
        matrix_file.unlink()
