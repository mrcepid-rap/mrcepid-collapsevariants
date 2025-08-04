import os
import shutil
from pathlib import Path

import pandas as pd
import pytest
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
from scipy.io import mmread

from collapsevariants.utilities.collapse_logger import CollapseLOGGER
from collapsevariants.utilities.collapse_utils import get_sample_ids
from collapsevariants.utilities.parallelization_wrappers import generate_genotype_matrix
from collapsevariants.tool_parsers.mask_generators import generate_generic_masks, generate_snp_or_gene_masks
from collapsevariants.snp_list_generator.snp_list_generator import SNPListGenerator

# Set this flag to True if you want to keep (copy) the temporary output files
KEEP_TEMP = False


@pytest.fixture(scope="module")
def pipeline_data():
    """
    Provides a shared dictionary for storing intermediate file names
    produced by each stage of the pipeline.
    """
    return {}


@pytest.fixture
def temporary_path(tmp_path, monkeypatch):
    """
    Prepare a temporary working directory that contains a copy of the test_data
    directory, then change the working directory to it.

    If KEEP_TEMP is True, after the test the entire temporary directory will be copied
    to a folder 'temp_test_outputs' in the project root.
    """
    # Determine where the original test_data directory is located.
    # (Assumes it is at <project_root>/test_data)
    test_data_source = Path(__file__).parent / "test_data"

    # Create the destination folder inside the tmp_path.
    destination = tmp_path / "test_data"
    destination.parent.mkdir(parents=True, exist_ok=True)

    # Copy the entire test_data directory into the temporary directory.
    shutil.copytree(test_data_source, destination)

    # Change the current working directory to the temporary directory.
    monkeypatch.chdir(tmp_path)

    # Yield the temporary directory to the test.
    yield tmp_path

    # After the test, if KEEP_TEMP is True, copy the temporary directory to a persistent location.
    if KEEP_TEMP:
        persistent_dir = Path(__file__).parent / "temp_test_outputs" / tmp_path.name
        persistent_dir.parent.mkdir(exist_ok=True)
        shutil.copytree(tmp_path, persistent_dir, dirs_exist_ok=True)


# Validated test data:
test_dir = Path(__file__).parent
test_data_dir = test_dir / 'test_data/'
correct_log = test_data_dir / 'correct_log.txt'

# Filtering lists
snp_handler = InputFileHandler(test_data_dir / 'snp_list.v2.txt')
gene_enst_handler = InputFileHandler(test_data_dir / 'gene_list.ENST.txt')
gene_symbol_handler = InputFileHandler(test_data_dir / 'gene_list.SYMBOL.txt')

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


# ======================================================================
# First test: generate the masks and store file names in pipeline_data
# ======================================================================
@pytest.mark.parametrize(
    "filtering_expression, gene_list_handler, snp_list_handler, expected_matrix, expected_samples, expected_variants,"
    "output_prefix",
    [
        ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', None, None,
         (10000, 13592), (10000, 2), (13592, 5), 'HC_PTV-MAF_001'),
        ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', gene_enst_handler, None,
         (10000, 4), (10000, 2), (34, 5), 'HC_PTV-MAF_001'),
        (None, None, snp_handler,
         (10000, 826), (10000, 2), (826, 5), "HC_PTV-MAF_001")
    ]
)
def test_snp_and_gene_masks(tmp_path, pipeline_data: pytest.fixture, filtering_expression: str,
                            gene_list_handler: InputFileHandler, snp_list_handler: InputFileHandler,
                            expected_matrix, expected_samples, expected_variants,
                            output_prefix: str):
    """
    Test the generation of SNP and gene masks by checking that the output files
    exist and have the correct shape.
    """
    log_path = tmp_path / 'HC_PTV-MAF_01.log'
    test_log = CollapseLOGGER(log_path)

    snp_list_generator = SNPListGenerator(
        bgen_dict=bgen_dict,
        filtering_expression=filtering_expression,
        gene_list_handler=gene_list_handler,
        snp_list_handler=snp_list_handler,
        log_file=test_log
    )

    # Generate a sparse matrix for each BGEN file and store the results.
    genotype_index = {}
    for bgen_prefix, variant_list in snp_list_generator.genes.items():
        bgen_dict[bgen_prefix]['index'].get_file_handle()
        geno_matrix = generate_genotype_matrix(
            bgen_prefix,
            bgen_dict[bgen_prefix],
            variant_list,
            delete_on_complete=False  # Make sure we keep the files for testing
        )
        genotype_index[bgen_prefix] = geno_matrix

    # Get sample IDs (assume they are identical across BGEN files).
    sample_ids = get_sample_ids(list(bgen_dict.values())[0]['sample'].get_file_handle())
    n_samples = len(sample_ids)
    assert n_samples == 10000

    # Generate output files based on which filtering list is provided.
    if snp_list_handler:
        output_files = generate_snp_or_gene_masks(
            genes=snp_list_generator.genes,
            genotype_index=genotype_index,
            sample_ids=sample_ids,
            output_prefix=output_prefix,
            bgen_type='SNP'
        )
        # Check existence of the 3 expected output files.
        for i in range(3):
            assert output_files[i].exists(), f"File {output_files[i]} does not exist."

        # Check that the output files are the correct shape.
        matrix_outfile = mmread(output_files[0])
        assert matrix_outfile.shape == expected_matrix, f"Expected {expected_matrix} got {matrix_outfile.shape}"

        sample_outfile = pd.read_csv(output_files[1], sep='\t')
        assert sample_outfile.shape == expected_samples, f"Expected {expected_samples} got {sample_outfile.shape}"

        variant_outfile = pd.read_csv(output_files[2], sep='\t')
        assert variant_outfile.shape == expected_variants, f"Expected {expected_variants} got {variant_outfile.shape}"

    elif gene_list_handler:
        output_files = generate_snp_or_gene_masks(
            snp_list_generator.genes, genotype_index,
            sample_ids, output_prefix, 'GENE'
        )
        for i in range(3):
            assert output_files[i].exists(), f"File {output_files[i]} does not exist."

        matrix_outfile = mmread(output_files[0])
        assert matrix_outfile.shape == expected_matrix, f"Expected {expected_matrix} got {matrix_outfile.shape}"

        sample_outfile = pd.read_csv(output_files[1], sep='\t')
        assert sample_outfile.shape == expected_samples, f"Expected {expected_samples} got {sample_outfile.shape}"

        variant_outfile = pd.read_csv(output_files[2], sep='\t')
        assert variant_outfile.shape == expected_variants, f"Expected {expected_variants} got {variant_outfile.shape}"

    else:
        output_files = generate_generic_masks(
            snp_list_generator.genes, genotype_index, sample_ids, output_prefix
        )
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
        # Verify that each file type has the expected number of files.
        for file_type, expected_count in file_types:
            files = [path for path in output_files if str(path).endswith(file_type)]
            assert len(files) == expected_count, f"Expected {expected_count} {file_type} files, got {len(files)}"

    # Ensure the key exists in pipeline_data
    if 'output_snp_and_gene_masks' not in pipeline_data:
        pipeline_data['output_snp_and_gene_masks'] = []

    # Append new files from this run
    pipeline_data['output_snp_and_gene_masks'].extend(
        [tmp_path / output_file.name for output_file in output_files])

    for file in pipeline_data['output_snp_and_gene_masks']:
        assert file.exists(), f"File not found when storing: {file}"


# ======================================================================
# Second test: verify file consistency using the stored pipeline_data.
# ======================================================================
@pytest.mark.parametrize(
    "filtering_expression, gene_list_handler, snp_list_handler, output_prefix",
    [
        ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', None, None, 'HC_PTV-MAF_001'),
        ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', gene_enst_handler, None, 'HC_PTV-MAF_001'),
        (None, None, snp_handler, 'HC_PTV-MAF_001')
    ]
)
def test_check_file_consistency(temporary_path, pipeline_data, filtering_expression,
                                gene_list_handler: InputFileHandler, snp_list_handler: InputFileHandler,
                                output_prefix: str):
    """
    Check that the files outputted by the previous test are consistent.
    This test rebases the file names stored in pipeline_data to the current temporary_path,
    then performs several sanity checks on sample and variant files.

    Note: If you don't run the tests from test_snp_and_gene_masks, this test will fail because the data does not exist.
    This includes if the previous tests passed, and only this method was run!
    """
    stored_files = pipeline_data.get('output_snp_and_gene_masks')
    # Verify that files exist before storing them
    for file in stored_files:
        assert file.exists(), f"File does not exist: {file}"
    assert stored_files is not None, "Missing output files in pipeline_data"

    # let's copy these into the temporary directory, otherwise the function won't work
    for file in stored_files:
        target = temporary_path / file.name
        if not target.exists():
            shutil.copy(file, target)

    # -----------------------------
    # Sanity check on sample tables:
    # -----------------------------
    snp_samples = pd.read_csv(f'{output_prefix}.SNP.STAAR.samples_table.tsv', sep='\t')
    gene_samples = pd.read_csv(f'{output_prefix}.GENE.STAAR.samples_table.tsv', sep='\t')
    # Using equals() to compare the entire column (instead of any())
    assert snp_samples['sampID'].equals(gene_samples['sampID']), "SNP and GENE sample IDs do not match."

    bolt_sample = pd.read_csv(f'{output_prefix}.chr1_chunk1.BOLT.sample', sep=' ')
    # bolt has an extra row at the top, so let's delete it to make sure the dataframes are matching
    bolt_sample = bolt_sample.iloc[1:].reset_index(drop=True)
    staar_sample = pd.read_csv(f'{output_prefix}.chr1_chunk1.STAAR.samples_table.tsv', sep='\t')
    assert bolt_sample['ID_1'].equals(staar_sample['sampID']), "BOLT and STAAR sample IDs do not match."

    staar_sample_1 = pd.read_csv(f'{output_prefix}.chr1_chunk1.STAAR.samples_table.tsv', sep='\t')
    staar_sample_2 = pd.read_csv(f'{output_prefix}.chr1_chunk2.STAAR.samples_table.tsv', sep='\t')
    assert staar_sample_1['sampID'].equals(staar_sample_2['sampID']), "Chunk1 and Chunk2 sample IDs differ."

    # -----------------------------
    # Sanity check on variant tables:
    # -----------------------------
    snp_variants = pd.read_csv(f'{output_prefix}.SNP.STAAR.variants_table.tsv', sep='\t')
    assert (snp_variants['POS'] < 90000000).all(), "Some SNP variant positions are >= 90000000."

    staar_variants_1 = pd.read_csv(f'{output_prefix}.chr1_chunk1.STAAR.variants_table.tsv', sep='\t')
    assert ((staar_variants_1['POS'] > 0) & (staar_variants_1['POS'] < 30000000)).all(), \
        "Chunk1 variant positions are out of expected range (0, 30000000)."

    staar_variants_2 = pd.read_csv(f'{output_prefix}.chr1_chunk2.STAAR.variants_table.tsv', sep='\t')
    assert ((staar_variants_2['POS'] > 30000000) & (staar_variants_2['POS'] < 60000000)).all(), \
        "Chunk2 variant positions are out of expected range (30000000, 60000000)."

    # -----------------------------
    # Sanity check on gene/SNP filtering:
    # -----------------------------
    for bgen, file in bgen_dict.items():
        if filtering_expression is not None:
            # For a mask filter, read in the original VEP file.
            vep = pd.read_csv(file['vep'].get_file_handle(), sep='\t')
            # Apply the filter manually.
            filtered_df = vep[(vep['PARSED_CSQ'] == "PTV") &
                              (vep['LOFTEE'] == 'HC') &
                              (vep['MAF'] < 0.001)]
            # Read the output file from the pipeline.
            annot = pd.read_csv(f'{output_prefix}.{bgen}.REGENIE.setListFile.txt',
                                sep='\t', header=None)
            # Compare the number of unique genes.
            assert len(filtered_df['ENST'].unique()) == len(annot.iloc[:, 0].unique()), \
                f"Filtered gene counts differ for BGEN prefix {bgen}."
        else:
            # If using a SNP list, merge the original VEP file with the SNP list.
            vep = pd.read_csv(file['vep'].get_file_handle(), sep='\t')
            snps_to_filter = pd.read_csv(snp_handler.get_file_handle(), header=None)
            snps_to_filter[['CHROM', 'POS', 'REF', 'ALT']] = snps_to_filter.iloc[:, 0].str.split(':', expand=True)
            # Ensure proper dtypes for merging.
            for col in ['CHROM', 'POS']:
                vep[col] = vep[col].astype(int)
                snps_to_filter[col] = snps_to_filter[col].astype(int)
            for col in ['REF', 'ALT']:
                vep[col] = vep[col].astype(str)
                snps_to_filter[col] = snps_to_filter[col].astype(str)
            filtered_df2 = vep.merge(
                snps_to_filter[['CHROM', 'POS', 'REF', 'ALT']],
                on=['CHROM', 'POS', 'REF', 'ALT'],
                how='inner'
            )
            annot = pd.read_csv(f'{output_prefix}.{bgen}.STAAR.variants_table.tsv', sep='\t')
            assert len(filtered_df2) == len(annot['ENST'].unique()), \
                f"Filtered SNP counts differ for BGEN prefix {bgen}."
