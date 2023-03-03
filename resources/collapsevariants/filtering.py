from pathlib import Path

from general_utilities.association_resources import run_cmd


def run_filtering(bgenprefix: str, chromosome: str, file_prefix: str) -> int:
    """A wrapper method around bgenix to filter WES .bgen files to variants as requested by SNPListGenerator

    This class will take a set of SNPs, hardcoded in this applet into the file 'include_snps.txt', and filter the
    provided bgen files down to only this set of SNPs in bgenv1.2, 'ref-last', 8-bit format. This class also creates
    a companion .sample file, a new .bgi index, and calculates the total number of variants parsed by the class.

    :param bgenprefix: The prefix name of the bgen file. Will often be identical to `chromosome`, but is included as a
        separate parameter to ensure compatibility with potential other data sources.
    :param chromosome: The chromosome currently being processed. This must be the short form of the chromosome name
        (e.g., '1' not 'chr1').
    :param file_prefix: The prefix provided at runtime for the final output name.

    :return: Total number of variants passing provided filters for this chromosome
    """

    # Simple bgenix command that includes variants from the filtering expression and just outputs a new "filtered"
    # bgen file
    cmd = f'bgenix -g /test/{bgenprefix}.bgen -incl-rsids /test/include_snps.txt > {file_prefix}.{chromosome}.bgen'
    run_cmd(cmd, is_docker=True, docker_image="egardner413/mrcepid-burdentesting:latest")

    cmd = f'cp {bgenprefix}.sample {file_prefix}.{chromosome}.sample'
    run_cmd(cmd, is_docker=False)

    cmd = f'bgenix -index -g /test/{file_prefix}.{chromosome}.bgen'
    run_cmd(cmd, is_docker=True, docker_image="egardner413/mrcepid-burdentesting:latest")

    # This just helps to get the total number of variants:
    cmd = f'bgenix -list -g /test/{file_prefix}.{chromosome}.bgen > {file_prefix}.{chromosome}.snps'
    run_cmd(cmd, is_docker=True, docker_image="egardner413/mrcepid-burdentesting:latest")

    total_vars = 0
    with Path(f'{file_prefix}.{chromosome}.snps').open('r') as snp_file:
        for line in snp_file:
            line = line.rstrip()
            if '#' not in line and 'alternate_ids' not in line:
                total_vars += 1
        snp_file.close()

    return total_vars
