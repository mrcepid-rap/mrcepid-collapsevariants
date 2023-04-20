import csv
from pathlib import Path
from typing import Tuple, TypedDict, List, Dict

from general_utilities.association_resources import run_cmd


class GeneDict(TypedDict):
    """A TypedDict containing information about all variants collapsed into a given gene

    :cvar CHROM: The chromosome this gene is on
    :cvar poss: A list of coordinates for all variants within this gene
    :cvar varIDs: A list of all variant IDs within this gene
    :cvar min_poss: The minimum coordinate for a variant within this gene. Is not a perfect proxy for gene start, but
        is good enough for the purposes of the code that needs that information.
    """

    CHROM: str
    poss: List[int]
    varIDs: List[str]
    min_poss: int


def parse_filters_SAIGE(file_prefix: str, chromosome: str) -> Tuple[Dict[str, GeneDict], Dict[str, str]]:
    """Generate input format files that can be provided to SAIGE

    For the way SAIGE is implemented downstream of this applet, SAIGE requires a standard .bcf file and a 'groupFile'.
    This groupFile is a simple tab-delimited file with individual genes (here defined as ENSTs) as the first column,
    followed by all the variant IDs that should be included in that gene per our mask definition.

    Variant IDs are slightly different from that included in other files produced by this applet, in that they must
    follow the format of CHR:POS_REF/ALT rather than CHR:POS:REF:ALT as defined in the original .bgen files produced
    prior to running this applet.

    This method returns two dictionaries:

    1. keys equal to ENST IDs and values equal to the GeneDict class
    2. keys equal to variant IDs and values equal to ENST IDs

    :param file_prefix: A name to append to beginning of output files.
    :param chromosome: The chromosome currently being processed. This must be the short form of the chromosome name
        (e.g., '1' not 'chr1').
    :return: A tuple containing two dictionaries of ENSTs mapped to variants and variants mapped to ENSTs, respectively
    """

    # Easier to run SAIGE with a BCF file as I already have that pipeline set up
    cmd = f'plink2 --memory 9000 --threads 1 --bgen /test/{file_prefix}.{chromosome}.bgen \'ref-last\' ' \
          f'--sample /test/{file_prefix}.{chromosome}.sample ' \
          f'--export bcf ' \
          f'--out /test/{file_prefix}.{chromosome}.SAIGE'
    run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting:latest')

    # and index...
    cmd = f'bcftools index /test/{file_prefix}.{chromosome}.SAIGE.bcf'
    run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting:latest')

    # Need to make the SAIGE groupFile. I can use the file 'snp_ENST.txt' created above to generate this...
    with Path('snp_ENST.txt').open('r') as snp_reader,\
            Path(f'{file_prefix}.{chromosome}.SAIGE.groupFile.txt').open('w') as output_setfile_SAIGE:

        genes: Dict[str, GeneDict] = dict()
        snp_gene_map: Dict[str, str] = dict()
        snp_csv = csv.DictReader(snp_reader, delimiter='\t')
        for snp in snp_csv:
            if snp['chrom'] == chromosome:
                snp_gene_map[snp['varID']] = snp['ENST']
                if snp['ENST'] in genes:
                    genes[snp['ENST']]['poss'].append(int(snp['pos']))
                    genes[snp['ENST']]['varIDs'].append(snp['varID'])
                else:
                    genes[snp['ENST']] = {'CHROM': snp['chrom'],
                                          'poss': [int(snp['pos'])],
                                          'varIDs': [snp['varID']]}

        for gene in genes:
            min_pos = min(genes[gene]['poss'])
            genes[gene]['min_poss'] = min_pos

        for gene in genes:
            # This is just using *args to place the four values that will always be here as {0} .. {3} automatically
            # into string formatting. Could probably write it a more functional way, but don't want to risk
            # disturbing this code that I know works properly
            id_string = "\t".join(["{0}:{1}_{2}/{3}".format(*item2) for item2 in
                                   [item.split(":") for item in genes[gene]['varIDs']]])
            output_setfile_SAIGE.write(f'{gene}\t{id_string}\n')

        Path(f'{file_prefix}.{chromosome}.SAIGE.log').unlink()

    return genes, snp_gene_map
