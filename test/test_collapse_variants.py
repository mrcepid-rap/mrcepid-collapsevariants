import filecmp
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.io import mmread
from scipy.sparse import csr_matrix

from collapsevariants.collapse_logger import CollapseLOGGER
from collapsevariants.collapse_utils import generate_csr_matrix_from_bgen, get_sample_ids, check_matrix_stats
from collapsevariants.parallelization_wrappers import generate_snp_or_gene_masks, generate_generic_masks
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


@pytest.mark.parametrize("filtering_expression, gene_list_path, snp_list_path, expected_num_sites",
                         [('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', None, None, 13592),
                          ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', gene_enst_path, None, 34),
                          (None, None, snp_path, 826)])
def test_csr_matrix_generation(tmp_path: Path, filtering_expression: str, gene_list_path: Path,
                               snp_list_path: Path, expected_num_sites: int):
    """This tests both :func:`generate_csr_matrix_from_bgen` and :func:`check_matrix_stats` functions.

    :param tmp_path: temporary pytest path for the logs
    :param filtering_expression: the filtering expression that is used as a mask
    :param gene_list_path: filepath to a text file containing a gene list to be used in the mask
    :param snp_list_path: filepath to a text file containing a variant list to be used as a mask
    :param expected_num_sites: expected number of sites post-masking
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


@pytest.mark.parametrize(
    "filtering_expression, gene_list_path, snp_list_path, expected_matrix, expected_samples, expected_variants",
    [
        ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', None, None, (10000, 13592), (10000, 2), (13592, 5)),
        ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', gene_enst_path, None, (10000, 34), (10000, 2), (34, 5)),
        (None, None, snp_path, (10000, 826), (10000, 2), (826, 5))
    ]
)
def test_snp_and_gene_masks(tmp_path: Path, filtering_expression: str, gene_list_path: Path,
                            snp_list_path: Path, expected_matrix: int, expected_samples: int, expected_variants: int):
    """
    Test the generation of SNP and gene masks.

    This function tests the generation of SNP and gene masks by verifying the output files and their contents.
    It checks the shape of the output matrices, sample files, and variant files to ensure they match the expected values.

    Parameters:
    tmp_path (Path): Temporary pytest path for the logs.
    filtering_expression (str): The filtering expression used as a mask.
    gene_list_path (Path): Filepath to a text file containing a gene list to be used in the mask.
    snp_list_path (Path): Filepath to a text file containing a variant list to be used in the mask.
    expected_matrix (int): Expected shape of the output matrix.
    expected_samples (int): Expected shape of the output sample file.
    expected_variants (int): Expected shape of the output variant file.

    Asserts:
    - The output files exist.
    - The output files have the correct shape.
    - The genes in the output files are filtered as required.
    """

    log_path = tmp_path / 'HC_PTV-MAF_01.log'
    test_log = CollapseLOGGER(log_path)

    # Has to be inside the test function since I have to re-open every time the test is run
    # Also make sure this is rb, as we wrap gzip around this in the function.
    vep_dict = {bgen_prefix: bgen_info['vep'].open('rb') for bgen_prefix, bgen_info in bgen_dict.items()}

    snp_list_generator = SNPListGenerator(vep_dict=vep_dict, filtering_expression=filtering_expression,
                                          gene_list_path=gene_list_path, snp_list_path=snp_list_path, log_file=test_log)

    # Generate a sparse matrix of genotypes from a BGEN file
    genotype_index = {}
    for bgen_prefix, variant_list in snp_list_generator.genes.items():
        geno_matrix = generate_csr_matrix_from_bgen(variant_list,
                                                    bgen_dict[bgen_prefix]['bgen'],
                                                    bgen_dict[bgen_prefix]['sample'])
        genotype_index[bgen_prefix] = geno_matrix

    # Get sample IDs for a single BGEN file (they are all the same)
    sample_ids = get_sample_ids(list(bgen_dict.values())[0]['sample'], )
    n_samples = len(sample_ids)
    assert n_samples == 10000

    # Now we need to filter the bgen files to only include the variants we want to keep
    if snp_list_path:
        output_files = generate_snp_or_gene_masks(genes=snp_list_generator.genes, genotype_index=genotype_index,
                                                  sample_ids=sample_ids, output_prefix='testing_output',
                                                  bgen_type='SNP')

        # Check the output files
        assert output_files[0].exists()
        assert output_files[1].exists()
        assert output_files[2].exists()

        # Check the output files are the correct shape
        matrix_outfile = mmread(output_files[0])
        assert matrix_outfile.shape == expected_matrix

        sample_outfile = pd.read_csv(output_files[1], sep='\t')
        assert sample_outfile.shape == expected_samples

        variant_outfile = pd.read_csv(output_files[2], sep='\t')
        assert variant_outfile.shape == expected_variants

    elif gene_list_path:
        output_files = generate_snp_or_gene_masks(snp_list_generator.genes, genotype_index,
                                                  sample_ids, 'testing_output', 'GENE')

        # Check the output files
        assert output_files[0].exists()
        assert output_files[1].exists()
        assert output_files[2].exists()

        # Check the output files are the correct shape
        matrix_outfile = mmread(output_files[0])
        assert matrix_outfile.shape == expected_matrix

        sample_outfile = pd.read_csv(output_files[1], sep='\t')
        assert sample_outfile.shape == expected_samples

        variant_outfile = pd.read_csv(output_files[2], sep='\t')
        assert variant_outfile.shape == expected_variants

    else:
        output_files = generate_generic_masks(snp_list_generator.genes, genotype_index, sample_ids, 'testing_output')

        file_types = [
            ('.bgen', 3),
            ('.bgen.bgi', 3),
            ('.sample', 3),
            ('.annotationFile.txt', 3),
            ('.maskfile.txt', 3),
            ('.setListFile.txt', 3),
            ('.groupFile.txt', 3),
            ('.samples_table.tsv', 3),
            ('.variants_table.tsv', 3)
        ]

        # Iterate over each file type and its expected count
        for file_type, expected_count in file_types:
            files = [path for path in output_files if str(path).endswith(file_type)]
            assert len(files) == expected_count


@pytest.mark.parametrize(
    "filtering_expression, gene_list_path, snp_list_path",
    [
        ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', None, None),
        ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', gene_enst_path, None),
        (None, None, snp_path)
    ]
)
def test_check_file_consistency(filtering_expression, gene_list_path, snp_list_path):
    """
    This is to check the files outputted by the tests above and to make sure that we are happy with them
    NOTE: I created this as a separate function because it seemed like there was a lag with writing the files in the
    function above.

    TODO: Maybe we can discuss what checks we want to include, this is a very rough go
    Would be good to maybe read in the original VCF files, apply the filtering, and make sure we get
    the correct number of genes
    """
    current_directory = Path(os.getcwd())

    # check to make sure the sample IDs are consistent
    # we're not checking all the files here - this is more of a sanity check
    snp_samples = pd.read_csv(Path(current_directory / 'testing_output.SNP.STAAR.samples_table.tsv'), sep='\t')
    gene_samples = pd.read_csv(Path(current_directory / 'testing_output.GENE.STAAR.samples_table.tsv'), sep='\t')
    assert snp_samples['sampID'].any() == gene_samples['sampID'].any()

    bolt_sample = pd.read_csv(Path(current_directory / 'testing_output.chr1_chunk1.BOLT.sample'), sep=' ')
    starr_sample = pd.read_csv(Path(current_directory / 'testing_output.chr1_chunk1.STAAR.samples_table.tsv'), sep='\t')
    assert bolt_sample['ID_1'].any() == starr_sample['sampID'].any()

    staar_sample_1 = pd.read_csv(Path(current_directory / 'testing_output.chr1_chunk1.STAAR.samples_table.tsv'),
                                 sep='\t')
    staar_sample_2 = pd.read_csv(Path(current_directory / 'testing_output.chr1_chunk2.STAAR.samples_table.tsv'),
                                 sep='\t')
    assert staar_sample_1['sampID'].any() == staar_sample_2['sampID'].any()

    # check that the variants tables are correct
    # we previously split the genetic data into chunks
    # this here makes sure that the chunking is still in effect (i.e. no funny business has occurred)
    snp_variants = pd.read_csv(Path(current_directory / 'testing_output.SNP.STAAR.variants_table.tsv'), sep='\t')
    assert snp_variants['POS'].all() < 90000000

    staar_variants_1 = pd.read_csv(Path(current_directory / 'testing_output.chr1_chunk1.STAAR.variants_table.tsv'),
                                   sep='\t')
    assert (staar_variants_1['POS'] > 0).all() and (
            staar_variants_1['POS'] < 30000000).all(), "POS column contains values outside the valid range."

    staar_variants_2 = pd.read_csv(Path(current_directory / 'testing_output.chr1_chunk2.STAAR.variants_table.tsv'),
                                   sep='\t')
    assert (staar_variants_2['POS'] > 30000000).all() and (
            staar_variants_1['POS'] < 60000000).all(), "POS column contains values outside the valid range."

    # let's also make sure that the genes in our output files are being
    # filtered as required using a manual test
    for bgen, file in bgen_dict.items():
        if filtering_expression is not None:  # if using the mask
            # read in the original VEP file
            vep = pd.read_csv(file['vep'], sep='\t')

            # filter is manually
            # 'PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001'
            filtered_df = vep[(vep['PARSED_CSQ'] == "PTV") & (vep['LOFTEE'] == 'HC') & (vep['MAF'] < 0.001)]

            # read in the output from collapsevariants
            annot = pd.read_csv(test_dir / f'testing_output.{bgen}.REGENIE.setListFile.txt', sep='\t', header=None)

            # see if they are matching
            assert len(filtered_df['ENST'].unique()) == len(annot.iloc[:, 0].unique())

        else:  # if using a SNP list
            # read in the original VEP file
            vep = pd.read_csv(file['vep'], sep='\t')

            # read in the SNPs that we are filtering on
            snps_to_filter = pd.read_csv(snp_path, header=None)
            # get the SNPs into the same format as in the output file
            snps_to_filter[['CHROM', 'POS', 'REF', 'ALT']] = snps_to_filter.iloc[:, 0].str.split(':', expand=True)

            # need to make sure the data-type is matching
            # (NOTE: should this go into the main codebase, or it doesn't really matter?)
            vep['CHROM'] = vep['CHROM'].astype(int)
            vep['POS'] = vep['POS'].astype(int)
            vep['REF'] = vep['REF'].astype(str)
            vep['ALT'] = vep['ALT'].astype(str)

            snps_to_filter['CHROM'] = snps_to_filter['CHROM'].astype(int)
            snps_to_filter['POS'] = snps_to_filter['POS'].astype(int)
            snps_to_filter['REF'] = snps_to_filter['REF'].astype(str)
            snps_to_filter['ALT'] = snps_to_filter['ALT'].astype(str)

            # merge the SNP list with our VEP file
            filtered_df2 = vep.merge(
                snps_to_filter[['CHROM', 'POS', 'REF', 'ALT']], on=['CHROM', 'POS', 'REF', 'ALT'], how='inner')

            # read in the output from collapsevariants
            annot = pd.read_csv(test_dir / f'testing_output.{bgen}.STAAR.variants_table.tsv', sep='\t')

            # see if they are matching
            assert len(filtered_df2) == len(annot['ENST'].unique())

    # delete_test_files(test_dir)


def delete_test_files(directory):
    """
    Delete all the files after we are done testing them
    """
    # Use glob to find files starting with "testing_output"
    files_to_delete = glob.glob(os.path.join(directory, 'testing_output*'))

    # Iterate and delete each file
    for file_path in files_to_delete:
        try:
            os.remove(file_path)
        except Exception as e:
            print(f"Error deleting {file_path}: {e}")
