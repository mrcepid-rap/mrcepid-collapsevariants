import pickle
from pathlib import Path
from typing import Dict, List, Tuple

import dxpy
import pandas as pd
from general_utilities.import_utils.file_handlers.export_file_handler import ExportFileHandler
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
from general_utilities.job_management.joblauncher_factory import joblauncher_factory
from scipy.io import mmwrite, mmread
from scipy.sparse import csr_matrix, hstack

from collapsevariants.tool_parsers.bolt_parser import BOLTParser
from collapsevariants.tool_parsers.regenie_parser import REGENIEParser
from collapsevariants.tool_parsers.saige_parser import SAIGEParser
from collapsevariants.tool_parsers.staar_parser import STAARParser
from collapsevariants.utilities.collapse_utils import GenotypeInfo


# def generate_generic_masks(genes: Dict[str, pd.DataFrame], genotype_index: Dict[str, Tuple[csr_matrix, Dict[str, GenotypeInfo]]],
#                            sample_ids: List[str], output_prefix: str) -> List[Path]:
#     """Wrapper to help generate output files for each tool.
#
#     :param genes: A dictionary containing the genes to collapse with keys of the BGEN file prefixes and values of
#         a Pandas DataFrame containing per-variant information.
#     :param genotype_index: A dictionary containing values of csr_matrix and keys of each BGEN file prefix.
#     :param sample_ids: A list of sample IDs for processing this file.
#     :param output_prefix: A string representing the prefix of the output files.
#     :return: A list of Path objects representing the output files created by the implementing classes.
#     """
#
#     # Generate output files for each tool
#     output_files = []
#     tool_methods = [BOLTParser, SAIGEParser, REGENIEParser, STAARParser]
#     for tool in tool_methods:
#
#         tool_instance = tool(genes, genotype_index, sample_ids, output_prefix)
#         output_files.extend(tool_instance.get_output_files())
#
#     return output_files

def generate_generic_masks(genes: Dict[str, pd.DataFrame], genotype_index: Dict[str, Tuple[csr_matrix, Dict[str, GenotypeInfo]]],
                           sample_ids: List[str], output_prefix: str) -> List[Path]:
    """Wrapper to help generate output files for each tool.

    :param genes: A dictionary containing the genes to collapse with keys of the BGEN file prefixes and values of
        a Pandas DataFrame containing per-variant information.
    :param genotype_index: A dictionary containing values of csr_matrix and keys of each BGEN file prefix.
    :param sample_ids: A list of sample IDs for processing this file.
    :param output_prefix: A string representing the prefix of the output files.
    :return: A list of Path objects representing the output files created by the implementing classes.
    """

    # set the launcher
    launcher = joblauncher_factory()

    # set the exporter
    exporter = ExportFileHandler(delete_on_upload=False)

    for chunk, df in genes.items():

        # get the gene data for this key
        gene_path = f"{chunk}.csv"
        df.to_csv(gene_path, index=False)
        gene_path = exporter.export_files(gene_path)
        # get the matrix for this key
        mmwrite(f"{chunk}.mtx", genotype_index[chunk][0])
        matrix_path = exporter.export_files(f"{chunk}.mtx")
        # expor the summary dict
        with open(f"{chunk}_summary_dict.pkl", "wb") as f:
            pickle.dump(genotype_index[chunk][1], f)
        summary_path = exporter.export_files(f"{chunk}_summary_dict.pkl")
        # get the samples for this key
        with open(f'sample_list_{chunk}.txt', 'w') as f:
            for item in sample_ids:
                f.write(f"{item}\n")
        sample_path = exporter.export_files(f'sample_list_{chunk}.txt')

        launcher.launch_job(function=multithread_generic_mask_generation,
                            inputs={
                                'chunk': chunk,
                                'gene_path': gene_path,
                                'matrix_path': matrix_path,
                                'summary_path': summary_path,
                                'sample_path': sample_path,
                                'output_prefix': output_prefix
                            },
                            outputs=['output_files']
                            )
    launcher.submit_and_monitor()
    output_files = []
    for result in launcher:
        output_files.extend(result['output_files'])

    return output_files


@dxpy.entry_point('multithread_generic_mask_generation')
def multithread_generic_mask_generation(chunk: str, gene_path, matrix_path, summary_path,
                                        sample_path, output_prefix: str):
    """Placeholder for future multithreading implementation of generic mask generation."""

    # read out genes data in
    genes = {}
    gene_path = InputFileHandler(gene_path).get_file_handle()
    df = pd.read_csv(gene_path)
    genes[chunk] = df

    # read the matrix data in
    matrix_path = InputFileHandler(matrix_path).get_file_handle()
    genotype_index = {}
    matrix = mmread(matrix_path)
    matrix = csr_matrix(matrix)  # Convert to subscriptable format
    # Load summary_dict from file
    summary_path = InputFileHandler(summary_path).get_file_handle()
    with open(summary_path, "rb") as f:
        summary_dict = pickle.load(f)
    genotype_index[chunk] = (matrix, summary_dict)

    # read the sample data in
    sample_path = InputFileHandler(sample_path).get_file_handle()
    with open(sample_path, 'r') as f:
        loaded_lst = [line.strip() for line in f]

    output_files = []
    tool_methods = [BOLTParser, SAIGEParser, REGENIEParser, STAARParser]
    for tool in tool_methods:

        tool_instance = tool(genes, genotype_index, loaded_lst, output_prefix)
        output_files.extend(tool_instance.get_output_files())

    return output_files


def generate_snp_or_gene_masks(genes: Dict[str, pd.DataFrame], genotype_index: Dict[str, Tuple[csr_matrix, Dict[str, GenotypeInfo]]],
                               sample_ids: List[str], output_prefix: str, bgen_type: str) -> List[Path]:
    """
    Wrapper similar to generate_generic_masks, but for SNP and GENE masks, creating output inputs for various tools.

    This method takes the output of the collapsing process and 'stacks' the resulting matrices. Since only one matrix
    is required for each GENE / SNP mask, we do not need to create multiple outputs for each .bgen as when running
    a filtering expression.

    :param genes: A dictionary containing the genes to collapse with keys of the BGEN file prefixes and values of
        a Pandas DataFrame containing per-variant information.
    :param genotype_index: A dictionary containing values of csr_matrix and keys of each BGEN file prefix.
    :param sample_ids: A list of sample IDs for processing this file.
    :param output_prefix: A string representing the prefix of the output files.
    :param bgen_type: A string representing the type of BGEN file (e.g., 'SNP' or 'GENE').
    :return: A list of Path objects representing the output files created by the implementing classes.
    """
    # Need to concatenate all matrices into a single matrix while ensuring concatenation is done in the same order for
    # variant indices AND genotype matrices

    # 1. Collect submatrices and variant indices
    final_variant_index_list = []
    matrix_list = []

    for bgen_prefix, variant_index in genes.items():
        current_matrix, _ = genotype_index[bgen_prefix]
        # Collect the variant index DataFrame
        final_variant_index_list.append(variant_index)
        # Collect the sparse matrix for later concatenation
        matrix_list.append(current_matrix)

    # 2. Concatenate all variant indices
    final_variant_index = pd.concat(final_variant_index_list)

    # 3. Perform one hstack on the list of sparse matrices
    final_genotype_matrix = hstack(matrix_list)

    # 4. Write the final data to disk (Matrix Market format)
    matrix_output_path = Path(f'{output_prefix}.{bgen_type}.STAAR.mtx')

    ## make sure the output is in the correct format
    mmwrite(matrix_output_path, final_genotype_matrix)

    sample_output_path = STAARParser.make_samples_dict(output_prefix, bgen_type, sample_ids)
    variant_output_path = STAARParser.make_variants_dict(output_prefix, bgen_type, final_variant_index)

    return [matrix_output_path, sample_output_path, variant_output_path]
