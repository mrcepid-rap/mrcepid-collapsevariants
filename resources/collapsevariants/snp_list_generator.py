import dxpy
import pandas as pd

from pathlib import Path
from typing import List, Dict

from ingest_data import BGENIndex
from collapse_logger import CollapseLOGGER


# The whole purpose of this class is to process the self.variant_index pandas DataFrame.
# Runs a filtering expression provided by the user to filter to a subset of variants from the associated filtered VCF
# file. The expression MUST be in a format parsable by pandas.query
class SNPListGenerator:

    def __init__(self, bgen_index: Dict[str, BGENIndex], filtering_expression: str, found_snps: bool, found_genes: bool,
                 LOG_FILE: CollapseLOGGER):

        self.bgen_index = bgen_index
        self.filtering_expression = filtering_expression
        self.found_snps = found_snps
        self.found_genes = found_genes
        self.LOG_FILE = LOG_FILE

        self.variant_index = self._load_variant_index()
        self._query_variant_index()
        self.chromosomes = self._write_variant_lists()
        self.total_sites = self._check_bgen_stats()

    # Load the entire variant index stored in the .vep tsv files
    def _load_variant_index(self) -> pd.DataFrame:

        variant_index = []
        # Open all chromosome indicies and load them into a list
        # Specify dtype for SIFT/POLYPHEN as pandas returns warnings when loading these due to weird formating
        for chromosome in self.bgen_index.keys():
            variant_index.append(
                pd.read_csv(f'filtered_bgen/chr{chromosome}.filtered.vep.tsv.gz', sep="\t",
                            dtype={'SIFT': str, 'POLYPHEN': str}))

        variant_index = pd.concat(variant_index)
        variant_index = variant_index.set_index('varID')

        return variant_index

    # Actually run the query and modify the variant_index loaded above based on this/these queries
    def _query_variant_index(self):

        # Three possibilities, each with an if/ifelse:
        # 1. Filtering Expression + Gene List = Select specific genes to collapse on
        # 2. Filtering Expression = Any variant which qualifies under the given filtering expression
        # 3. SNPList = Select specific SNPs to collapse on

        # If using a gene list with a filtering expression...
        if self.filtering_expression is not None and self.found_genes:
            my_genelist_file = open("gene_list.genes", "r")
            genelist = my_genelist_file.readlines()
            genelist = [gene.rstrip() for gene in genelist]

            # First filter to the gene Symbols we care about
            self.variant_index = self.variant_index.query("SYMBOL == @genelist")
            all_genelist = list(self.variant_index["SYMBOL"].unique())
            found_genes_symbol_present = list(set.intersection(set(genelist), set(all_genelist)))

            # Now further filter down to the filtering_expression we are interested in
            self.variant_index = self.variant_index.query(self.filtering_expression)

            # Set all gene ENSTs to ENST99999999999
            # This is a dummy value so that association_testing knows we are running a gene list
            self.variant_index['ENST'] = 'ENST99999999999'

            # Following code is just to report stats about genes that were/were not found
            queried_genes = []
            self.LOG_FILE.write_header('Per-gene stats')
            for key, value in self.variant_index['SYMBOL'].value_counts().iteritems():
                queried_genes.append(key)
                self.LOG_FILE.write_int(key, value)
            self.LOG_FILE.write_spacer()

            for gene in genelist:
                if gene not in found_genes_symbol_present:
                    print(f'Asked for gene – {gene} – gene symbol was not found in the variant index.')
                elif gene not in queried_genes:
                    print(f'Asked for gene – {gene} – no variants in this gene after applying filtering expression.')

        # This one is simple – just need to identify ALL variants that fit a given expression
        elif self.filtering_expression is not None:
            self.variant_index = self.variant_index.query(self.filtering_expression)

        # This is if we want to extract specific SNPs
        elif self.found_snps:
            self.LOG_FILE.write_generic('Variants not Found:')
            with Path('snp_list.snps').open('r') as snplist:
                snps = []
                not_found = []
                for snp in snplist:
                    snp = snp.rstrip()
                    if snp in self.variant_index.index:
                        snps.append(snp)
                    else:
                        not_found.append(snp)
                        self.LOG_FILE.write_generic(snp)
                snplist.close()
            self.LOG_FILE.write_spacer()

            if len(snps) == 0:
                raise dxpy.AppError('No SNPs remain after using SNPlist!')
            else:
                # Print to logfile variants that we found...
                self.LOG_FILE.write_generic('Variants found')
                for snp in snps:
                    self.LOG_FILE.write_generic(snp)
                self.LOG_FILE.write_spacer()

            self.LOG_FILE.write_int('Total variants found', len(snps))
            self.LOG_FILE.write_spacer()

            # And finally extract variants here
            self.variant_index = self.variant_index.loc[snps]

            # Here we create a fake 'GENE' to collapse on later to make the process of generating a mask easier.
            # We also use this gene ID to notify runassociationtesting that we are doing a SNP-based mask
            self.variant_index['ENST'] = 'ENST00000000000'

    # Write the lists that will be used for filtering and generate a list of chromosomes that are valid for actual
    # processing.
    def _write_variant_lists(self) -> List[str]:

        # Now prep a list of variants to include in any analysis. We need two files:
        # 1. Assigns a variant to a gene (gene_id_file)
        # 2. Just includes variant IDs (snp_id_file)
        snp_id_file = open('include_snps.txt', 'w')
        gene_id_file = open('snp_ENST.txt', 'w')
        gene_id_file.write("varID\tCHROM\tPOS\tENST\n")
        for row in self.variant_index.iterrows():
            # row[0] in this context is the varID since it is the 'index' in the pandas DataFrame
            # All other information is stored in a dictionary that is list element [1]
            snp_id_file.write(f'{row[0]}\n')
            gene_id_file.write(f'{row[0]}\t{row[1]["CHROM"].replace("chr", "")}\t{row[1]["POS"]}\t{row[1]["ENST"]}')
        snp_id_file.close()
        gene_id_file.close()

        # Get chromosomes we want to actually run to save time downstream:
        chromosomes = self.variant_index['CHROM'].value_counts()
        chromosomes = chromosomes.index.to_list()
        chromosomes = [str.replace(chrom, "chr", "") for chrom in chromosomes]

        return chromosomes

    # This gets information relating to included variants in bgen format files (per-gene)
    def _check_bgen_stats(self) -> int:

        self.LOG_FILE.write_header('Overall Statistics')
        total_sites = self.variant_index['CHROM'].count()  # Assign as a variable so I can return it below for further checking
        self.LOG_FILE.write_int('Total number of variants', total_sites)
        self.LOG_FILE.write_float('Maximum missingness', self.variant_index['F_MISSING'].max() * 100)
        self.LOG_FILE.write_scientific('Maximum Allele Frequency', self.variant_index['AF'].max(skipna=True))
        self.LOG_FILE.write_float('Minimum CADD Score', self.variant_index['CADD'].min(skipna=True))
        self.LOG_FILE.write_float('Minimum REVEL Score', self.variant_index['REVEL'].min(skipna=True))
        self.LOG_FILE.write_int('Number of NA REVEL Scores',
                                self.variant_index[self.variant_index['REVEL'].isna()]['CHROM'].count())
        self.LOG_FILE.write_int('Total number of PASS variants',
                                self.variant_index[self.variant_index['FILTER'] == 'PASS']['CHROM'].count())
        self.LOG_FILE.write_int('Total number of non-PASS variants',
                                self.variant_index[self.variant_index['FILTER'] != 'PASS']['CHROM'].count())

        # LOFTEE:
        for key, value in self.variant_index['LOFTEE'].value_counts().iteritems():
            key = f'Number of LOFTEE {key}'
            self.LOG_FILE.write_int(key, value)
        self.LOG_FILE.write_spacer()

        # Parsed Consequences:
        self.LOG_FILE.write_header('Consequence statistics')
        for key, value in self.variant_index['PARSED_CSQ'].value_counts().iteritems():
            key = f'Number of Parsed Consequence – {key}'
            self.LOG_FILE.write_int(key, value)

        # VEP Consequences:
        for key, value in self.variant_index['CSQ'].value_counts().iteritems():
            key = f'Number of VEP Consequence - {key}'
            self.LOG_FILE.write_int(key, value, is_vep=True)
        self.LOG_FILE.write_spacer()

        return total_sites
