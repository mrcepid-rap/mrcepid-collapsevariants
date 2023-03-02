import csv
import dxpy
import pandas as pd

from pathlib import Path
from typing import List, Dict

from ingest_data import BGENIndex
from collapse_logger import CollapseLOGGER
from general_utilities.mrc_logger import MRCLogger


class SNPListGenerator:
    """Generate a SNP list to create a variant Mask

    The purpose of this class is to process the self.variant_index pandas DataFrame into a final set of filtered SNPs.
    This can be done in one of three ways:

    1. Via a filtering expression – Runs a filtering expression provided by the user to filter to a subset of variants from the associated filtered VCF
        file. The expression MUST be in a format parsable by pandas.query
    2. Via a gene list – Uses a filtering expression (as above in [1]), but also takes a list of genes and collapses
        them into a single fake 'GENE' with the ID 'ENST99999999999'
    3. Via a SNP list – Takes a list of SNPs and collapses them into a single fake 'GENE' with the ID 'ENST00000000000'.
        Unlike the previous two methods, this approach DOES NOT use a filtering expression and relies on the user to
        filter SNPs.

    :param bgen_index: A dict that contains chromosomes as keys and a BGENIndex TypedDict containing information on
        filepaths containing filtered & annotated variants.
    :param filtering_expression: A pandas.query() compatible expression to filter variants on, 'None' if not provided
    :param found_snps: Were SNPs found to filter against in collapse_variants.ingest_data (IngestData class)
    :param found_genes: Were SNPs found to filter against in collapse_variants.ingest_data (IngestData class)
    :param LOG_FILE: A class of CollapseLOGGER to store information on variant filtering
    """

    def __init__(self, bgen_index: Dict[str, BGENIndex], filtering_expression: str, found_snps: bool, found_genes: bool,
                 LOG_FILE: CollapseLOGGER):

        self._logger = MRCLogger(__name__).get_logger()
        self.bgen_index = bgen_index
        self.filtering_expression = filtering_expression
        self.found_snps = found_snps
        self.found_genes = found_genes
        self.LOG_FILE = LOG_FILE

        self.variant_index = self._query_variant_index()
        self.chromosomes = self._write_variant_lists()
        self.total_sites = self._check_bgen_stats()

    # Load the entire variant index stored in the .vep tsv files
    def _load_variant_index(self) -> pd.DataFrame:
        """Load vep annotated tsv.gz files into a pandas DataFrame

        Will iterate over all chromosomes that were found by collapse_variants.ingest_data, NOT the actual chromosomes
        in the human genome (e.g., 1..23, X, Y).

        :return: A pandas.DataFrame containing variants loaded from all provided chromosomes
        """

        variant_index = []
        # Open all chromosome indicies and load them into a list
        # Specify dtype for SIFT/POLYPHEN as pandas returns warnings when loading these due to weird formating
        for chromosome in self.bgen_index.keys():
            variant_index.append(
                pd.read_csv(f'filtered_bgen/chr{chromosome}.filtered.vep.tsv.gz', sep='\t',
                            dtype={'SIFT': str, 'POLYPHEN': str}))

        variant_index = pd.concat(variant_index)
        variant_index = variant_index.set_index('varID')

        return variant_index

    # Actually run the query and modify the variant_index loaded above based on this/these queries
    def _query_variant_index(self) -> pd.DataFrame:
        """Query self.variant_index based on the logic outlined in the __init__ method for this class

        <<INPUT>>

        As described for the class, there are three possibilities, each with an if/ifelse, and in this logical order:
        1. Filtering Expression + Gene List = Select specific genes to collapse on
        2. Filtering Expression = Any variant which qualifies under the given filtering expression
        3. SNPList = Select specific SNPs to collapse on

        :return: A pandas.DataFrame containing variants filtered based on provided input parameters
        """

        # First load
        variant_index = self._load_variant_index()

        # 1. Filtering expression + Gene List
        if self.filtering_expression is not None and self.found_genes:
            variant_index = self._query_gene_list(variant_index)

        # 2. Filtering expression
        elif self.filtering_expression is not None:
            variant_index = self._query_filtering_expression(variant_index)

        # 3. SNP List
        elif self.found_snps:
            variant_index = self._query_snp_list(variant_index)

        return variant_index

    def _query_gene_list(self, variant_index: pd.DataFrame) -> pd.DataFrame:
        """Filter variants by a provided gene list

        :param variant_index: A pandas.DataFrame containing variants loaded from all provided chromosomes
        :return: A modified version of the variant_index pandas.DataFrame AFTER filtering on provided Gene list
        """

        with Path('gene_list.genes').open('r') as my_genelist_file:
            genelist = my_genelist_file.readlines()
            genelist = [gene.rstrip() for gene in genelist]

        # First filter to the gene Symbols we care about
        variant_index = variant_index.query('SYMBOL == @genelist')
        all_genelist = list(variant_index['SYMBOL'].unique())
        found_genes_symbol_present = list(set.intersection(set(genelist), set(all_genelist)))

        # Now further filter down to the filtering_expression we are interested in
        variant_index = variant_index.query(self.filtering_expression)

        # Set all gene ENSTs to ENST99999999999
        # This is a dummy value so that association_testing knows we are running a gene list
        variant_index['ENST'] = 'ENST99999999999'

        # Following code is just to report stats about genes that were/were not found
        queried_genes = []
        self.LOG_FILE.write_header('Per-gene stats')
        for key, value in variant_index['SYMBOL'].value_counts().items():
            queried_genes.append(key)
            self.LOG_FILE.write_int(key, value)
        self.LOG_FILE.write_spacer()

        for gene in genelist:
            if gene not in found_genes_symbol_present:
                self._logger.warning(f'Asked for gene – {gene} – gene symbol was not found in the variant index.')
            elif gene not in queried_genes:
                self._logger.warning(f'Asked for gene – {gene} – no variants in this gene after applying '
                                     f'filtering expression.')

        return variant_index

    def _query_filtering_expression(self, variant_index: pd.DataFrame) -> pd.DataFrame:
        """Filter variants by a pandas.DataFrame.query compatible expression

        This one is simple – just need to identify ALL variants that fit a given expression. All we do is
        DataFrame.query('expression') and return the resulting pared-down pd.DataFrame.

        :param variant_index: A pandas.DataFrame containing variants loaded from all provided chromosomes
        :return: A modified version of the variant_index pandas.DataFrame AFTER filtering on
            provided filtering_expression
        """

        variant_index = variant_index.query(self.filtering_expression)
        return variant_index

    def _query_snp_list(self, variant_index: pd.DataFrame) -> pd.DataFrame:
        """Filter variants by a provided SNP list

        :param variant_index: A pandas.DataFrame containing variants loaded from all provided chromosomes
        :return: A modified version of the variant_index pandas.DataFrame AFTER filtering on provided SNP list
        """
        self.LOG_FILE.write_generic('Variants not Found:')
        with Path('snp_list.snps').open('r') as snplist:
            snps = []
            not_found = []
            for snp in snplist:
                snp = snp.rstrip()
                if snp in variant_index.index:
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
        variant_index = variant_index.loc[snps]

        # Here we create a fake 'GENE' to collapse on later to make the process of generating a mask easier.
        # We also use this gene ID to notify runassociationtesting that we are doing a SNP-based mask
        variant_index['ENST'] = 'ENST00000000000'

        return variant_index

    # Write the lists that will be used for filtering and generate a list of chromosomes that are valid for actual
    # processing.
    def _write_variant_lists(self) -> List[str]:

        # Now prep a list of variants to include in any analysis. We need two files:
        # 1. Assigns a variant to a gene (gene_id_file)
        # 2. Just includes variant IDs (snp_id_file)
        with Path('include_snps.txt').open('w') as snp_id_file,\
                Path('snp_ENST.txt').open('w') as gene_id_file:

            snp_id_csv = csv.DictWriter(snp_id_file, delimiter="\t", fieldnames=['varID'], extrasaction='ignore')
            gene_id_csv = csv.DictWriter(gene_id_file, delimiter='\t', fieldnames=['varID', 'CHROM', 'POS', 'ENST'])
            gene_id_csv.writeheader()

            for row in self.variant_index.iterrows():
                # row[0] in this context is the varID since it is the 'index' in the pandas DataFrame
                # All other information is stored in a dictionary that is list element [1]
                line_dict = {'varID': row[0],
                             'CHROM': row[1]['CHROM'].replace('chr', ''),
                             'POS': row[1]['POS'],
                             'ENST': row[1]['ENST']}
                snp_id_csv.writerow(line_dict)
                gene_id_csv.writerow(line_dict)

        # Get chromosomes we want to actually run to save time downstream:
        chromosomes = self.variant_index['CHROM'].value_counts()
        chromosomes = chromosomes.index.to_list()
        chromosomes = [str.replace(chrom, 'chr', '') for chrom in chromosomes]

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
        for key, value in self.variant_index['LOFTEE'].value_counts().items():
            key = f'Number of LOFTEE {key}'
            self.LOG_FILE.write_int(key, value)
        self.LOG_FILE.write_spacer()

        # Parsed Consequences:
        self.LOG_FILE.write_header('Consequence statistics')
        for key, value in self.variant_index['PARSED_CSQ'].value_counts().items():
            key = f'Number of Parsed Consequence – {key}'
            self.LOG_FILE.write_int(key, value)

        # VEP Consequences:
        for key, value in self.variant_index['CSQ'].value_counts().items():
            key = f'Number of VEP Consequence - {key}'
            self.LOG_FILE.write_int(key, value, is_vep=True)
        self.LOG_FILE.write_spacer()

        return total_sites
