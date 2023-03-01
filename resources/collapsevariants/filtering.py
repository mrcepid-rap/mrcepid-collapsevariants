from pathlib import Path

from general_utilities.association_resources import run_cmd


class Filter:

    def __init__(self, bgenprefix: str, chromosome: str, file_prefix: str):

        self.num_variants = self._run_filtering(bgenprefix, chromosome, file_prefix)

    # This runs the per-chromosome side of filtering
    @staticmethod
    def _run_filtering(bgenprefix: str, chromosome: str, file_prefix: str) -> int:

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
