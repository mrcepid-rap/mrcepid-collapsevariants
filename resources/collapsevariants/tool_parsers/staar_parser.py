import os
import csv

from collapsevariants.collapse_resources import *


class STAARParser:

    def __init__(self, file_prefix: str, chromosome: str):

        self._parse_filters_STAAR(file_prefix, chromosome)

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
        samples_dict_file = file_prefix + "." + chromosome + '.samples.tsv'
        col_num = 1
        with open(samples_dict_file, 'w') as samples_dict_writer:
            samples_dict_writer.write('sampID\trow\n')
            sample_file = csv.DictReader(open(file_prefix + "." + chromosome + '.sample', 'r', newline='\n'),
                                         delimiter=" ",
                                         quoting=csv.QUOTE_NONE)
            for sample in sample_file:
                if sample['ID'] != "0": # This gets rid of the wierd header row in bgen sample files...
                    samples[sample['ID']] = col_num
                    samples_dict_writer.write('%s\t%i\n' % (sample['ID'], col_num))
                    col_num+=1
            samples_dict_writer.close()

        # ... Then a dictionary of variant column IDs
        # This file is also saved for runassociationtesting
        variants = {}
        variants_dict_file = file_prefix + "." + chromosome + '.variants_table.STAAR.tsv'
        row_num = 1
        with open(variants_dict_file, 'w') as variants_dict_writer:
            variants_file = open('snp_ENST.txt', 'r', newline='\n')
            variants_dict_writer.write('varID\tchrom\tpos\tENST\tcolumn\n')
            variants_csv = csv.DictReader(variants_file,
                                          delimiter="\t",
                                          quoting=csv.QUOTE_NONE)
            for var in variants_csv:
                if var['CHROM'] == chromosome or chromosome == 'SNP' or chromosome == 'GENE':
                    variants[var['varID']] = row_num
                    variants_dict_writer.write('%s\t%s\t%s\t%s\t%i\n' % (var['varID'], var['CHROM'], var['POS'], var['ENST'], row_num))
                    row_num += 1
            variants_dict_writer.close()

        # Need to get the file length of this file...
        # This code is janky AF
        proc = subprocess.Popen("wc -l " + file_prefix + "." + chromosome + ".parsed.txt", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        file_length = stdout.decode('utf-8')
        file_length = file_length.strip(' ').split(' ')[0]

        # And format it to R spec
        sparse_matrix = csv.DictReader(open(file_prefix + "." + chromosome + '.parsed.txt', 'r'), delimiter="\t",quoting=csv.QUOTE_NONE,fieldnames=['sample', 'varID', 'genotype'])
        matrix_file = file_prefix + "." + chromosome + '.STAAR.matrix.R.txt'
        with open(matrix_file, 'w') as matrix_file_writer:
            # Write the header in %MatrixMarket format
            matrix_file_writer.write('%%MatrixMarket matrix coordinate integer general\n')
            matrix_file_writer.write('%i %i %s\n' % ((col_num-1), (row_num-1), file_length))  # -1 because of how the iterator above works...
            # Write the individual cells in the matrix
            for row in sparse_matrix:
                gt_val = 1 if row['genotype'] == '0/1' else 2
                matrix_file_writer.write('%s %s %s\n' % (samples[row['sample']], variants[row['varID']], gt_val))
            matrix_file_writer.close()

        # And then run the R script `buildSTAARmatrix.R` to properly attach row and column names and save in .RDS format
        # for quick read/write in later parts of the pipeline.
        # See the script itself for how it works, but it takes 4 inputs:
        # 1. Samples dict (samples_dict_file above) – the row names of our final matrix
        # 2. Variants dict (variants_dict_file above) – the column names of our final matrix
        # 3. Matrix file (matrix_file above) – the cells in our final matrix, formatted in %MatrixMarket format for easy read
        # 4. Out file – The name of the .rds file for final output
        # And generates one output:
        # 1. A .rds file (named by out file) that can be read back into STAAR during mrcepid-runassociationtesting
        cmd = "Rscript /prog/buildSTAARmatrix.R /test/%s /test/%s /test/%s /test/%s" % \
              (samples_dict_file, variants_dict_file, matrix_file, file_prefix + '.' + chromosome + '.STAAR.matrix.rds')
        run_cmd(cmd, True)

        # Remove intermediate files so they aren't stored in the final tar:
        os.remove(samples_dict_file)
        os.remove(matrix_file)
