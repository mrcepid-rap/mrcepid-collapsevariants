import csv
import dxpy
import numpy as np
import pandas as pd

from enum import Enum, auto
from pathlib import Path
from typing import Tuple

from collapsevariants.ingest_data import IngestData
from collapsevariants.collapse_logger import CollapseLOGGER
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.mrc_logger import MRCLogger


class FilteringMode(Enum):
    GENE_LIST = auto()
    FILTERING_EXPRESSION = auto()
    SNP_LIST = auto()


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

    :param ingested_data: An object from CollapseIngester containing: 1) A dict that contains chromosomes as keys and
        a BGENIndex TypedDict containing information on filepaths containing filtered & annotated variants. 2) A
        pandas.query() compatible expression to filter variants on, 'None' if not provided. 3) Were SNPs
        found to filter against in collapse_variants.ingest_data (IngestData class) 4) Were SNPs found
        to filter against in collapse_variants.ingest_data (IngestData class)
    :param log_file: A class of CollapseLOGGER to store information on variant filtering
    """

    def __init__(self, ingested_data: IngestData, log_file: CollapseLOGGER):

        self._bgen_index = ingested_data.bgen_index
        self._filtering_expression = ingested_data.filtering_expression
        self._found_genes = ingested_data.found_genes  # This is a boolean if we found a gene list
        self._gene_list = None  # Placeholder for a numpy array of gene symbols
        self._found_snps = ingested_data.found_snps  # This is a boolean if we found a SNP list
        self._snp_list = None  # Placeholder for a numpy array of SNP IDs

        self._logger = MRCLogger(__name__).get_logger()
        self.log_file = log_file

        # Decide on the type of filtering we are doing (SNP, Gene, or expression)
        self._filtering_mode = self._decide_filtering_mode()

        # Iterate through all possible bgens in parallel and filter them
        self._variant_index = []
        thread_utility = ThreadUtility()
        for prefix in self._bgen_index.keys():
            thread_utility.launch_job(self._query_variant_index,
                                      prefix=prefix)

        self.chromosomes = dict()
        for thread in thread_utility:
            current_table, prefix, chrom, vars_found = thread.result
            if vars_found:
                self._variant_index.append(thread.result)
                if chrom not in self.chromosomes:
                    self.chromosomes[chrom] = [prefix]
                else:
                    self.chromosomes[chrom].append(prefix)

        self._variant_index = pd.concat(self._variant_index)

        # Write the variant lists to be used for filtering
        self.snp_enst_index_path, self.include_path = self._write_variant_lists()

        # Check the stats of the bgen files
        self.total_sites = self._check_bgen_stats()

    def _decide_filtering_mode(self) -> FilteringMode:

        if self._filtering_expression is not None and self._found_genes:

            with Path('gene_list.genes').open('r') as my_genelist_file:
                gene_list = my_genelist_file.readlines()
                self._gene_list = np.array([gene.rstrip() for gene in gene_list])

            return FilteringMode.FILTERING_EXPRESSION

        # 2. Filtering expression
        elif self._filtering_expression is not None:
            return FilteringMode.FILTERING_EXPRESSION

        # 3. SNP List
        elif self._found_snps:

            with Path('snp_list.snps').open('r') as my_snplist_file:
                snp_list = my_snplist_file.readlines()
                self._snp_list = np.array([snp.rstrip() for snp in snp_list])

            return FilteringMode.SNP_LIST

        else:
            raise ValueError('No filtering expression, SNP list, or Gene list provided!')

    @staticmethod
    def _load_variant_index(prefix: str) -> pd.DataFrame:
        """Load vep annotated *.tsv.gz files into a pandas DataFrame

        :return: A pandas.DataFrame containing variants loaded from all provided chromosomes
        """

        current_vep = pd.read_csv(f'{prefix}.filtered.vep.tsv.gz', sep='\t',
                                  dtype={'SIFT': str, 'POLYPHEN': str})
        variant_index = current_vep.set_index('varID')

        return variant_index

    # Actually run the query and modify the variant_index loaded above based on this/these queries
    def _query_variant_index(self, prefix: str) -> Tuple[pd.DataFrame, str, str, bool]:
        """Query self.variant_index based on the logic outlined in the __init__ method for this class

        As described in the __init__ for the class, there are three possibilities, each with an if/ifelse, and in this
        logical order:

        1. Filtering Expression + Gene List = Select specific genes to collapse on
        2. Filtering Expression = Any variant which qualifies under the given filtering expression
        3. SNPList = Select specific SNPs to collapse on

        :param prefix: Prefix of the file to load. Handled  by :func:`_load_variant_index`
        :return: A Tuple of four parts: 1) pandas.DataFrame containing variants from :param variant_index: filtered based on provided
            input parameters 2) The prefix of the file loaded 3) The chromosome of all variants found 4) A boolean
            indicating if any variants were found after filtering.
        """

        variant_index = self._load_variant_index(prefix)

        # 1. Filtering expression + Gene List
        if self._filtering_mode == FilteringMode.GENE_LIST:
            variant_index = self._query_gene_list(variant_index)

        # 2. Filtering expression
        elif self._filtering_mode == FilteringMode.FILTERING_EXPRESSION:
            variant_index = self._query_filtering_expression(variant_index)

        # 3. SNP List
        elif self._filtering_mode == FilteringMode.SNP_LIST:
            variant_index = self._query_snp_list(variant_index)

        poss_chrom = variant_index['CHROM'].unique()
        if len(poss_chrom) > 1:
            raise ValueError(f'More than one chromosome found in bgen {prefix}!')
        else:
            poss_chrom = poss_chrom[0]

        return variant_index, prefix, poss_chrom, len(variant_index) > 0

    def _query_gene_list(self, variant_index: pd.DataFrame) -> pd.DataFrame:
        """Filter variants by a provided gene list

        :param variant_index: A pandas.DataFrame containing variants loaded from all provided chromosomes
        :return: A modified version of the variant_index pandas.DataFrame AFTER filtering on provided Gene list
        """

        # First filter to the gene Symbols we care about
        variant_index = variant_index[variant_index['SYMBOL'].isin(self._gene_list)]

        # Now further filter down to the filtering_expression we are interested in
        variant_index = variant_index.query(self._filtering_expression)

        # Set all gene ENSTs to ENST99999999999
        # This is a dummy value so that association_testing knows we are running a gene list
        variant_index['ENST'] = 'ENST99999999999'

        return variant_index

    def _query_filtering_expression(self, variant_index: pd.DataFrame) -> pd.DataFrame:
        """Filter variants by a pandas.DataFrame.query compatible expression

        This one is simple – just need to identify ALL variants that fit a given expression. All we do is
        DataFrame.query('expression') and return the resulting pared-down pd.DataFrame.

        :param variant_index: A pandas.DataFrame containing variants loaded from all provided chromosomes
        :return: A modified version of the variant_index pandas.DataFrame AFTER filtering on
            provided filtering_expression
        """

        variant_index = variant_index.query(self._filtering_expression)
        return variant_index

    def _query_snp_list(self, variant_index: pd.DataFrame) -> pd.DataFrame:
        """Filter variants by a provided SNP list

        :param variant_index: A pandas.DataFrame containing variants loaded from all provided chromosomes
        :return: A modified version of the variant_index pandas.DataFrame AFTER filtering on provided SNP list
        """

        # And finally extract variants here
        variant_index = variant_index.loc[self._snp_list]

        # Here we create a fake 'GENE' to collapse on later to make the process of generating a mask easier.
        # We also use this gene ID to notify runassociationtesting that we are doing a SNP-based mask
        variant_index['ENST'] = 'ENST00000000000'

        return variant_index

    def _write_variant_lists(self) -> Tuple[Path, Path]:
        """Write the lists that will be used for filtering and generate a list of chromosomes that are valid for actual
            processing.

        Generate two files that will be used later:

        1. A list of variants -> Genes (gene_ENST.txt)
        2. Just a list of variant IDs (include_snps.txt)

        Also keep track of chromosomes that we found so we don't process chromosomes without a qualifying variant.

        :return: A Tuple consisting of two Path objects – one for the snp_enst_index (links snps to ENSTs) and one
            for the include_path (just a list of snps).
        """

        # Sort the variant index:
        self._variant_index.sort_values(by=['CHROM', 'POS'], inplace=True)

        # Now prep a list of variants to include in any analysis. We need two files:
        # 1. Assigns a variant to a gene (gene_id_file)
        # 2. Just includes variant IDs (snp_id_file)
        snp_enst_index_path = Path('gene_ENST.txt')
        self._variant_index[['CHROM', 'POS', 'ENST']].to_csv(snp_enst_index_path, sep='\t', header=True, index=True)

        include_path = Path('include_snps.txt')
        self._variant_index[[]].to_csv(include_path, sep='\t', header=False, index=True)  # Double bracket drops all cols

        return snp_enst_index_path, include_path

    def _check_bgen_stats(self) -> int:
        """Gets information relating to included variants in bgen format files (per-gene)

        This method just uses standard pandas.DataFrame methods to query the underlying table to get stats and report
        them in the applet log_file

        :return: The total number of sites after filtering
        """

        self.log_file.write_header('Overall Statistics')
        total_sites = self._variant_index[
            'CHROM'].count()  # Assign as a variable so I can return it below for further checking
        self.log_file.write_int('Total number of variants', total_sites)
        self.log_file.write_float('Maximum missingness', self._variant_index['F_MISSING'].max() * 100)
        self.log_file.write_scientific('Maximum Allele Frequency', self._variant_index['AF'].max(skipna=True))
        self.log_file.write_float('Minimum CADD Score', self._variant_index['CADD'].min(skipna=True))
        self.log_file.write_float('Minimum REVEL Score', self._variant_index['REVEL'].min(skipna=True))
        self.log_file.write_int('Number of NA REVEL Scores',
                                self._variant_index[self._variant_index['REVEL'].isna()]['CHROM'].count())
        self.log_file.write_int('Total number of PASS variants',
                                self._variant_index[self._variant_index['FILTER'] == 'PASS']['CHROM'].count())
        self.log_file.write_int('Total number of non-PASS variants',
                                self._variant_index[self._variant_index['FILTER'] != 'PASS']['CHROM'].count())

        # LOFTEE:
        for key, value in self._variant_index['LOFTEE'].value_counts().items():
            key = f'Number of LOFTEE {key}'
            self.log_file.write_int(key, value)
        self.log_file.write_spacer()

        # Parsed Consequences:
        self.log_file.write_header('Consequence statistics')
        for key, value in self._variant_index['PARSED_CSQ'].value_counts().items():
            key = f'Number of Parsed Consequence – {key}'
            self.log_file.write_int(key, value)

        # VEP Consequences:
        for key, value in self._variant_index['CSQ'].value_counts().items():
            key = f'Number of VEP Consequence - {key}'
            self.log_file.write_int(key, value, is_vep=True)
        self.log_file.write_spacer()

        # Append filtering information on SNP / GENE lists if required:
        if self._filtering_mode == FilteringMode.SNP_LIST or self._filtering_mode == FilteringMode.GENE_LIST:
            self._check_list_filtering()

        return total_sites

    def _check_list_filtering(self) -> None:
        """This method will check against either the original SNP or GENE list to determine if variants / genes were
        found in the final set of filtered variants
        """

        ## For SNP list
        if self._filtering_mode == FilteringMode.SNP_LIST:

            # Make a boolean array of SNPs that are in the final list
            snp_boolean = np.isin(self._snp_list, self._variant_index.index)

            # Then use the boolean to subset the snp array
            missing_snps = self._snp_list[snp_boolean == False]
            found_snps = self._snp_list[snp_boolean == True]

            if len(found_snps) == 0:
                raise dxpy.AppError('No SNPs remain after using SNPlist!')

            else:

                self.log_file.write_int('Total variants found', len(found_snps))
                self.log_file.write_int('Total variants not found', len(missing_snps))
                # Print to logfile variants that we found...
                self.log_file.write_generic('Variants not Found:')
                for snp in missing_snps:
                    self.log_file.write_generic(snp)
                self.log_file.write_spacer()

        elif self._filtering_mode == FilteringMode.GENE_LIST:
            ## For Gene list
            # Following code is just to report stats about genes that were/were not found
            all_genelist = self._variant_index['SYMBOL'].unique()  # np.ndarray
            found_genes_symbol_present = all_genelist[np.isin(self._gene_list, all_genelist)]

            queried_genes = []
            self.log_file.write_header('Per-gene stats')
            for key, value in self._variant_index['SYMBOL'].value_counts().items():
                queried_genes.append(key)
                self.log_file.write_int(key, value)
            self.log_file.write_spacer()

            for gene in self._gene_list:
                if gene not in found_genes_symbol_present:
                    self._logger.warning(f'Asked for gene – {gene} – gene symbol was not found in the variant index.')
                elif gene not in queried_genes:
                    self._logger.warning(f'Asked for gene – {gene} – no variants in this gene after applying '
                                         f'filtering expression.')
