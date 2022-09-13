from collapse_resources import *

class Filter:

    def __init__(self, bgenprefix: str, chromosome: str, file_prefix: str):

        self.num_variants = self._run_filtering(bgenprefix, chromosome, file_prefix)

    # This runs the per-chromosome side of filtering
    @staticmethod
    def _run_filtering(bgenprefix: str, chromosome: str, file_prefix: str) -> int:

        # Simple bgenix command that includes variants from the filtering expression and just outputs a new "filtered"
        # bgen file
        cmd = "bgenix -g /test/" + bgenprefix + ".bgen -incl-rsids /test/include_snps.txt > " + file_prefix + "." + chromosome + ".bgen"
        run_cmd(cmd, True)

        cmd = "cp " + bgenprefix + ".sample " + file_prefix + "." + chromosome + ".sample"
        run_cmd(cmd)

        cmd = "bgenix -index -g /test/" + file_prefix + "." + chromosome + ".bgen"
        run_cmd(cmd, True)

        # This just helps to get the total number of variants:
        cmd = "bgenix -list -g /test/" + file_prefix + "." + chromosome + ".bgen > " + file_prefix + "." + chromosome + ".snps"
        run_cmd(cmd, True)

        total_vars = 0
        with open(file_prefix + "." + chromosome + '.snps') as snp_file:
            for line in snp_file:
                line = line.rstrip()
                if '#' not in line and 'alternate_ids' not in line:
                    total_vars += 1
            snp_file.close()

        return total_vars
