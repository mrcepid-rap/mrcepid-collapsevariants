import os
import shutil
from pathlib import Path

import pandas as pd
import pytest
from scipy.io import mmread

from collapsevariants.collapse_logger import CollapseLOGGER
from collapsevariants.collapse_utils import generate_csr_matrix_from_bgen, get_sample_ids
from collapsevariants.parallelization_wrappers import generate_snp_or_gene_masks, generate_generic_masks
from collapsevariants.snp_list_generator import SNPListGenerator

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
        print(f"Temporary output files have been copied to: {persistent_dir}")


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


# ======================================================================
# First test: generate the masks and store file names in pipeline_data
# ======================================================================
@pytest.mark.parametrize(
    "filtering_expression, gene_list_path, snp_list_path, expected_matrix, expected_samples, expected_variants",
    [
        ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', None, None,
         (10000, 13592), (10000, 2), (13592, 5)),
        ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', gene_enst_path, None,
         (10000, 34), (10000, 2), (34, 5)),
        (None, None, snp_path,
         (10000, 826), (10000, 2), (826, 5))
    ]
)
def test_snp_and_gene_masks(temporary_path, pipeline_data, filtering_expression: str,
                            gene_list_path: Path, snp_list_path: Path,
                            expected_matrix, expected_samples, expected_variants):
    """
    Test the generation of SNP and gene masks by checking that the output files
    exist and have the correct shape.
    """

    log_path = temporary_path / 'HC_PTV-MAF_01.log'
    test_log = CollapseLOGGER(log_path)

    # Open the VEP files from bgen_dict; we assume that bgen_dict is a dict of info per BGEN file.
    vep_dict = {bgen_prefix: bgen_info['vep'].open('rb')
                for bgen_prefix, bgen_info in bgen_dict.items()}

    snp_list_generator = SNPListGenerator(
        vep_dict=vep_dict,
        filtering_expression=filtering_expression,
        gene_list_path=gene_list_path,
        snp_list_path=snp_list_path,
        log_file=test_log
    )

    # Generate a sparse matrix for each BGEN file and store the results.
    genotype_index = {}
    for bgen_prefix, variant_list in snp_list_generator.genes.items():
        geno_matrix = generate_csr_matrix_from_bgen(
            variant_list,
            bgen_dict[bgen_prefix]['bgen'],
            bgen_dict[bgen_prefix]['sample']
        )
        genotype_index[bgen_prefix] = geno_matrix

    # Get sample IDs (assume they are identical across BGEN files).
    sample_ids = get_sample_ids(list(bgen_dict.values())[0]['sample'])
    n_samples = len(sample_ids)
    assert n_samples == 10000

    # Generate output files based on which filtering list is provided.
    if snp_list_path:
        output_files = generate_snp_or_gene_masks(
            genes=snp_list_generator.genes,
            genotype_index=genotype_index,
            sample_ids=sample_ids,
            output_prefix='testing_output',
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

    elif gene_list_path:
        output_files = generate_snp_or_gene_masks(
            snp_list_generator.genes, genotype_index,
            sample_ids, 'testing_output', 'GENE'
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
            snp_list_generator.genes, genotype_index, sample_ids, 'testing_output'
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
        [temporary_path / output_file.name for output_file in output_files])

    # Debugging print
    print("Updated stored files in pipeline_data:", pipeline_data['output_snp_and_gene_masks'])

    for file in pipeline_data['output_snp_and_gene_masks']:
        assert file.exists(), f"File not found when storing: {file}"


# ======================================================================
# Second test: verify file consistency using the stored pipeline_data.
# ======================================================================
@pytest.mark.parametrize(
    "filtering_expression, gene_list_path, snp_list_path",
    [
        ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', None, None),
        ('PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001', gene_enst_path, None),
        (None, None, snp_path)
    ]
)
def test_check_file_consistency(temporary_path, pipeline_data, filtering_expression, gene_list_path, snp_list_path):
    """
    Check that the files outputted by the previous test are consistent.
    This test rebases the file names stored in pipeline_data to the current temporary_path,
    then performs several sanity checks on sample and variant files.
    """
    stored_files = pipeline_data.get('output_snp_and_gene_masks')
    # Verify that files exist before storing them
    for file in stored_files:
        assert file.exists(), f"File does not exist: {file}"
    assert stored_files is not None, "Missing output files in pipeline_data"

    # let's copy these into the current directory, otherwise the function won't work
    current_dir = Path(os.getcwd())
    for file in stored_files:
        target = current_dir / file.name
        if not target.exists():
            shutil.copy(file, target)
            print(f"Copied {file} to {target}")

    # -----------------------------
    # Sanity check on sample tables:
    # -----------------------------
    snp_samples = pd.read_csv('testing_output.SNP.STAAR.samples_table.tsv', sep='\t')
    gene_samples = pd.read_csv('testing_output.GENE.STAAR.samples_table.tsv', sep='\t')
    # Using equals() to compare the entire column (instead of any())
    assert snp_samples['sampID'].equals(gene_samples['sampID']), "SNP and GENE sample IDs do not match."

    bolt_sample = pd.read_csv('testing_output.chr1_chunk1.BOLT.sample', sep=' ')
    # bolt has an extra row at the top, so let's delete it to make sure the dataframes are matching
    bolt_sample = bolt_sample.iloc[1:].reset_index(drop=True)
    staar_sample = pd.read_csv('testing_output.chr1_chunk1.STAAR.samples_table.tsv', sep='\t')
    assert bolt_sample['ID_1'].equals(staar_sample['sampID']), "BOLT and STAAR sample IDs do not match."

    staar_sample_1 = pd.read_csv('testing_output.chr1_chunk1.STAAR.samples_table.tsv', sep='\t')
    staar_sample_2 = pd.read_csv('testing_output.chr1_chunk2.STAAR.samples_table.tsv', sep='\t')
    assert staar_sample_1['sampID'].equals(staar_sample_2['sampID']), "Chunk1 and Chunk2 sample IDs differ."

    # -----------------------------
    # Sanity check on variant tables:
    # -----------------------------
    snp_variants = pd.read_csv('testing_output.SNP.STAAR.variants_table.tsv', sep='\t')
    assert (snp_variants['POS'] < 90000000).all(), "Some SNP variant positions are >= 90000000."

    staar_variants_1 = pd.read_csv('testing_output.chr1_chunk1.STAAR.variants_table.tsv', sep='\t')
    assert ((staar_variants_1['POS'] > 0) & (staar_variants_1['POS'] < 30000000)).all(), \
        "Chunk1 variant positions are out of expected range (0, 30000000)."

    staar_variants_2 = pd.read_csv('testing_output.chr1_chunk2.STAAR.variants_table.tsv', sep='\t')
    assert ((staar_variants_2['POS'] > 30000000) & (staar_variants_2['POS'] < 60000000)).all(), \
        "Chunk2 variant positions are out of expected range (30000000, 60000000)."

    # -----------------------------
    # Sanity check on gene/SNP filtering:
    # -----------------------------
    for bgen, file in bgen_dict.items():
        if filtering_expression is not None:
            # For a mask filter, read in the original VEP file.
            vep = pd.read_csv(file['vep'], sep='\t')
            # Apply the filter manually.
            filtered_df = vep[(vep['PARSED_CSQ'] == "PTV") &
                              (vep['LOFTEE'] == 'HC') &
                              (vep['MAF'] < 0.001)]
            # Read the output file from the pipeline.
            annot = pd.read_csv(f'testing_output.{bgen}.REGENIE.setListFile.txt',
                                sep='\t', header=None)
            # Compare the number of unique genes.
            assert len(filtered_df['ENST'].unique()) == len(annot.iloc[:, 0].unique()), \
                f"Filtered gene counts differ for BGEN prefix {bgen}."
        else:
            # If using a SNP list, merge the original VEP file with the SNP list.
            vep = pd.read_csv(file['vep'], sep='\t')
            snps_to_filter = pd.read_csv(snp_path, header=None)
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
            annot = pd.read_csv(f'testing_output.{bgen}.STAAR.variants_table.tsv', sep='\t')
            assert len(filtered_df2) == len(annot['ENST'].unique()), \
                f"Filtered SNP counts differ for BGEN prefix {bgen}."
