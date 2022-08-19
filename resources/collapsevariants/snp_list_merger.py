import shutil

from tool_parsers.bolt_parser import *
from tool_parsers.staar_parser import *

class SNPMerger:

    def __init__(self, valid_chromosomes: list, file_prefix: str, found_genes: bool):

        self._merge_across_snplist(valid_chromosomes, file_prefix, found_genes)


    @staticmethod
    # This is used regardless of if a SNP of Gene list (apologies for naming...)
    # Note that the default is to merge a SNP list. An extra flag for found genes can be added to do a gene list
    def _merge_across_snplist(valid_chromosomes: list, file_prefix: str, found_genes: bool = False):

        # file prefix to distinguish output from gene list or snp list
        if found_genes:
            tar_type = "GENE"
        else:
            tar_type = "SNP"

        # First generate a list of vcfs that we are going to mash together and get variant names:
        cmd = "bcftools concat -Ob -o /test/" + file_prefix + "." + tar_type + ".SAIGE.pre.bcf "
        variant_IDs = ['ENST99999999999' if found_genes else 'ENST00000000000']  #1st ENST replacement
        snp_gene_map = {}
        for chrom in valid_chromosomes:
            cmd += "/test/" + file_prefix + "." + chrom + ".SAIGE.bcf "
            with open(file_prefix + "." + chrom + ".SAIGE.groupFile.txt", 'r') as group_file:
                for line in group_file:
                    line = line.rstrip()
                    variants = line.split("\t")
                    variant_IDs.extend(variants[1:len(variants)])
                    for variant in variants[1:len(variants)]:
                        bolt_format_ID = variant.replace('_',':').replace('/',':')
                        snp_gene_map[bolt_format_ID] = 'ENST00000000000'
                        if found_genes:                                     # 2nd ENST replacement
                            snp_gene_map[bolt_format_ID] = 'ENST99999999999'

        # Combine with bcftools concat
        run_cmd(cmd, True)
        # Make sure sorted properly...
        cmd = "bcftools sort -Ob -o /test/" + file_prefix + "." + tar_type + ".SAIGE.bcf /test/" + file_prefix + "." + tar_type + ".SAIGE.pre.bcf"
        run_cmd(cmd, True)
        os.remove(file_prefix + "." + tar_type + ".SAIGE.pre.bcf")
        # And index:
        cmd = "bcftools index /test/" + file_prefix + "." + tar_type + ".SAIGE.bcf"
        run_cmd(cmd, True)

        # Write new groupFile:
        with open(file_prefix + "." + tar_type + ".SAIGE.groupFile.txt", "w") as snp_groupfile:
            snp_groupfile.write("\t".join(variant_IDs))
            snp_groupfile.close()

        # Trick the already made BOLT code above to build a new merged BOLT file:
        genes = {}
        if found_genes:
            genes['ENST99999999999'] = {'CHROM': 1, 'min_poss': 1}
        else:
            genes['ENST00000000000'] = {'CHROM': 1, 'min_poss': 1}

        # This copy is slightly dodgy, as it assumes at least one chromosome has come through in the variable
        # 'chrom' from the loop above
        shutil.copy(file_prefix + "." + chrom + ".sample", file_prefix + "." + tar_type + ".sample")
        BOLTParser(file_prefix, tar_type, genes, snp_gene_map)

        # Trick the already made STAAR code above to build a new merged set of STAAR files
        STAARParser(file_prefix, tar_type)

        # Delete old files to avoid confusion:
        for chrom in valid_chromosomes:
            os.remove(file_prefix + "." + chrom + ".SAIGE.bcf")
            os.remove(file_prefix + "." + chrom + ".SAIGE.bcf.csi")
            os.remove(file_prefix + "." + chrom + ".SAIGE.groupFile.txt")
            os.remove(file_prefix + "." + chrom + ".BOLT.bgen")
            os.remove(file_prefix + "." + chrom + ".BOLT.sample")
            os.remove(file_prefix + "." + chrom + ".STAAR.matrix.rds")
            os.remove(file_prefix + "." + chrom + ".variants_table.STAAR.tsv")
