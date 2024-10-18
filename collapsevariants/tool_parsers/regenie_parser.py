import csv
from pathlib import Path


def parse_filters_REGENIE(tarball_prefix: str, chromosome: str) -> None:

    # This is used to print the set list file (2) below
    gene_dict = {}
    prefix = f'{tarball_prefix}.{chromosome}'

    # 1. Annotation File
    with Path(f'{prefix}.REGENIE.annotationFile.tsv').open('w') as annotation_file, \
            Path(f'{prefix}.variants_table.STAAR.tsv').open('r') as table_file:

        table_reader = csv.DictReader(table_file, delimiter='\t')
        annotation_writer = csv.DictWriter(annotation_file,
                                           delimiter='\t',
                                           fieldnames=['ID', 'ENST', 'annotation'],
                                           extrasaction='ignore',
                                           lineterminator='\n')  # REGENIE is very fussy about line terminators.
        last_var = None  # Need to check for small number of duplicate variants...
        for variant in table_reader:
            if last_var != variant['varID']:
                variant['annotation'] = tarball_prefix
                annotation_writer.writerow(variant)
                last_var = variant['varID']
                # And build gene_dict while we iterate...
                if variant['ENST'] in gene_dict:
                    gene_dict[variant['ENST']]['varIDs'].append(variant['varID'])
                else:
                    gene_dict[variant['ENST']] = {'chrom': variant['chrom'],
                                                  'pos': variant['pos'],
                                                  'varIDs': [variant['varID']],
                                                  'ENST': variant['ENST']}

    # 2. Set list file
    with Path(f'{prefix}.REGENIE.setListFile.tsv').open('w') as set_list_file:
        set_list_writer = csv.DictWriter(set_list_file,
                                         delimiter="\t",
                                         fieldnames=['ENST', 'chrom', 'pos', 'varIDs'],
                                         lineterminator='\n')
        for gene in sorted(list(gene_dict.values()), key=lambda item: item['pos']):
            gene['varIDs'] = ','.join(gene['varIDs'])
            set_list_writer.writerow(gene)

    # 3. This makes the mask name file. Just needs to be the name of the mask (tarball prefix) used in file #1
    with Path(f'{prefix}.REGENIE.maskfile.tsv').open('w') as mask_file:
        mask_file.write(f'{tarball_prefix}\t{tarball_prefix}\n')
