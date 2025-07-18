import csv

import pytest
from pathlib import Path

from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler

from collapsevariants.utilities.ingest_data import IngestData, download_bgen

test_dir = Path(__file__).parent
test_data_dir = test_dir / 'test_data/'

filtering_expression = 'PARSED_CSQ=="PTV" & LOFTEE=="HC" & MAF < 0.001'
snp_list = Path(__file__).parent / 'test_data' / 'snp_list.txt'
gene_list = Path(__file__).parent / 'test_data' / 'gene_list.SYMBOL.txt'


@pytest.fixture(scope="function")
def bgen_file(tmpdir) -> Path:

    bgen_tsv = tmpdir / 'bgen_index.tsv'

    with bgen_tsv.open('w') as bgen_writer:
        bgen_csv = csv.DictWriter(bgen_writer, fieldnames=['prefix', 'bgen_id', 'bgen_index_id', 'sample_id', 'vep_id'], delimiter="\t")
        bgen_csv.writeheader()

        for i in range(1,4):
            prefix = f'chr1_chunk{i}'
            bgen_csv.writerow({
                'prefix': prefix,
                'bgen_id': f'{test_data_dir}/{prefix}.bgen',
                'bgen_index_id': f'{test_data_dir}/{prefix}.bgen.bgi',
                'sample_id': f'{test_data_dir}/{prefix}.sample',
                'vep_id': f'{test_data_dir}/{prefix}.vep.tsv.gz',
            })

    return Path(bgen_tsv)


@pytest.mark.parametrize(
    "expression, snp_file, gene_file",
    [
        (filtering_expression, None, None),
        (filtering_expression, snp_list, None),
        (filtering_expression, None, gene_list),
    ]
)
def test_ingest_data(bgen_file, expression, snp_file, gene_file):

    ingested_data = IngestData(bgen_file, expression, snp_file, gene_file)
    expected_file_types = {'index', 'sample', 'bgen', 'vep'}

    assert len(ingested_data.sample_ids) == 10_000
    assert len(ingested_data.bgen_dict.keys()) == 3

    for bgen_prefix, bgen_info in ingested_data.bgen_dict.items():

        assert len(expected_file_types.difference(bgen_info.keys())) == 0

        for file_type in expected_file_types:
            assert isinstance(bgen_info[file_type], InputFileHandler)
            assert bgen_info[file_type].get_file_handle().exists(), f"{file_type} file {bgen_info[file_type].get_file_handle()} does not exist"

    if snp_file:
        assert ingested_data.snp_list_handler.get_file_handle().exists()

    if gene_file:
        assert ingested_data.gene_list_handler.get_file_handle().exists()

# def test_download_bgen(bgen_index):
#
#     download_bgen()