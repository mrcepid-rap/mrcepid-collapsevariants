import gzip
from typing import Dict, Optional

import pandas as pd
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.mrc_logger import MRCLogger

from collapsevariants.snp_list_generator.filtering_mode import FilteringMode
from collapsevariants.snp_list_generator.stat_dictionary import StatDictionary
from collapsevariants.utilities.collapse_logger import CollapseLOGGER
from collapsevariants.utilities.ingest_data import BGENIndex


class SNPListGenerator:
    """Generate a SNP list to create a variant Mask

    The purpose of this class is to process the provided VEP annotations for all provided bgen files into a final set
    of filtered SNPs. This can be done in one of three ways:

    1. Via a filtering expression – Runs a filtering expression provided by the user to filter to a subset of variants
        from the associated filtered VCF file. The expression MUST be in a format parsable by pandas.query
    2. Via a gene list – Uses a filtering expression (as above in [1]), but also takes a list of genes and collapses
        them into a single fake 'GENE' with the ID 'ENST99999999999'
    3. Via a SNP list – Takes a list of SNPs and collapses them into a single fake 'GENE' with the ID 'ENST00000000000'.
        Unlike the previous two methods, this approach DOES NOT use a filtering expression and relies on the user to
        filter SNPs.

    :param bgen_dict: A dict that contains the bgen file prefixes as keys and a BGENIndex describing the bgen file
        as values. BGENIndex includes all information about the bgen file, including the vep file to use for annotation.
    :param filtering_expression: A pandas.query() compatible expression to filter variants on, 'None' if not provided.
    :param gene_list_handler: A InputFileHandler to a list of genes to filter against.
    :param snp_list_handler: A InputFileHandler to a list of SNPs to filter against.
    :param log_file: A class of CollapseLOGGER to store information on variant filtering
    """

    def __init__(self, bgen_dict: Dict[str, BGENIndex], filtering_expression: str,
                 gene_list_handler: Optional[InputFileHandler], snp_list_handler: Optional[InputFileHandler],
                 log_file: CollapseLOGGER):

        self._logger = MRCLogger(__name__).get_logger()

        # Set defaults for all class variables
        self._bgen_dict = bgen_dict
        self._filtering_expression = filtering_expression

        self._gene_list_path = None if gene_list_handler is None else gene_list_handler.get_file_handle()  # This is a Path to a found Gene list
        self._gene_list = set()  # Placeholder for a numpy array of gene symbols
        self._found_gene_dict = dict()

        self._snp_list_path = None if snp_list_handler is None else snp_list_handler.get_file_handle()  # This is a Path to a found SNP list
        self._snp_list = set()  # Placeholder for a numpy array of SNP IDs
        self._found_snp_set = set()  # Placeholder for a numpy array of SNP IDs

        # stat_dict is a dictionary that will be used to store the stats for each variant index parsed
        self._stat_dict = StatDictionary(log_file)

        # Decide on the type of filtering we are doing (SNP, Gene, or expression)
        self._filtering_mode = self._decide_filtering_mode()

        # Iterate through all possible bgens in parallel and filter them
        thread_utility = ThreadUtility()
        for prefix, bgen_info in self._bgen_dict.items():
            thread_utility.launch_job(self._query_variant_index,
                                      vep_handle=bgen_info['vep'],
                                      prefix=prefix)

        # Next we want to take the filtered result and process into a dictionary with keys of chromosomes and values of
        # genes. Genes will also be a dictionary containing SNPs and positions for later filtering.
        self.genes = dict()
        for result_dict in thread_utility:
            if result_dict['vars_found']:

                self.genes[result_dict['prefix']] = self._make_gene_dict(result_dict['variant_index'])

        # Check the stats of the bgen files
        self._stat_dict.write_bgen_stats()
        self._stat_dict.check_list_filtering(self._filtering_mode,
                                             self._snp_list, self._gene_list,
                                             self._found_snp_set, self._found_gene_dict)
        self.total_sites = self._stat_dict.get_total_sites()

    def _decide_filtering_mode(self) -> FilteringMode:
        """Decide on the filtering mode based on the provided input

        :return: The filtering mode selected in the form of a FilteringMode Enum
        """

        if self._filtering_expression is not None and self._snp_list_path is None and self._gene_list_path is None:
            # For logging purposes output the filtering expression provided by the user
            self._logger.info('Current Filtering Expression:')
            self._logger.info(self._filtering_expression)
            filtering_mode = FilteringMode.FILTERING_EXPRESSION

        elif self._filtering_expression is not None and self._snp_list_path is None and self._gene_list_path is not None:
            self._logger.info('Gene list provided together with filtering expression - will filter for variant '
                              'classes of interest in gene set')
            self._logger.info('Current Filtering Expression:')
            self._logger.info(self._filtering_expression)

            # Read the gene list into a set
            if self._gene_list_path.exists():
                with self._gene_list_path.open('r') as my_genelist_file:
                    gene_list = my_genelist_file.readlines()
                    self._gene_list = set([gene.rstrip() for gene in gene_list])

                filtering_mode = FilteringMode.GENE_LIST

            else:
                raise FileNotFoundError(f'Provided SNP list {self._snp_list_path} does not exist!')

        elif self._filtering_expression is None and self._snp_list_path is not None and self._gene_list_path is None:
            self._logger.info('Using provided SNPlist to generate a mask...')

            # Read the SNP list into a set
            if self._snp_list_path.exists():
                with self._snp_list_path.open('r') as my_snplist_file:
                    snp_list = my_snplist_file.readlines()
                    self._snp_list = set([snp.rstrip() for snp in snp_list])

                filtering_mode = FilteringMode.SNP_LIST
            else:
                raise FileNotFoundError(f'Provided SNP list {self._snp_list_path} does not exist!')

        else:
            raise ValueError('Incorrect input for snp/gene/filtering expression provided... exiting!')

        self._logger.info(f'Selected filtering mode is: {filtering_mode.name.lower()}')
        return filtering_mode

    @staticmethod
    def _load_variant_index(vep_handle: InputFileHandler) -> pd.DataFrame:
        """Load vep annotated *.tsv.gz files into a pandas DataFrame

        :param vep_handle: Pre-opened IO to a vep index file
        :return: A pandas.DataFrame containing variants loaded from all provided chromosomes
        """

        current_vep = pd.read_csv(gzip.open(vep_handle.get_file_handle(), mode='rt'), sep="\t",
                                  index_col='varID',
                                  dtype={'SIFT': str, 'POLYPHEN': str, 'LOFTEE': str,
                                         'AA': str, 'AApos': str})

        return current_vep

    # Actually run the query and modify the variant_index loaded above based on this/these queries
    def _query_variant_index(self, vep_handle: InputFileHandler, prefix: str) -> dict:
        """Query self.variant_index based on the logic outlined in the __init__ method for this class

        As described in the __init__ for the class, there are three possibilities, each with an if/ifelse, and in this
        logical order:

        1. Filtering Expression + Gene List = Select specific genes to collapse on
        2. Filtering Expression = Any variant which qualifies under the given filtering expression
        3. SNPList = Select specific SNPs to collapse on

        :param vep_handle: InputFileHandler for the vep index to load. Handled  by :func:`_load_variant_index`
        :param prefix: The prefix of the bgen file loaded
        :return: A dict containing: 1) pandas.DataFrame containing variants from :param variant_index: filtered based on provided
            input parameters 2) The prefix of the file loaded 3) The chromosome of all variants found 4) A boolean
            indicating if any variants were found after filtering 5) The minimum position of the entire bgen file.
        """

        variant_index = self._load_variant_index(vep_handle)

        # 1. Filtering expression + Gene List
        if self._filtering_mode == FilteringMode.GENE_LIST:
            variant_index = self._query_gene_list(variant_index)

        # 2. Filtering expression
        elif self._filtering_mode == FilteringMode.FILTERING_EXPRESSION:
            variant_index = self._query_filtering_expression(variant_index)
            if variant_index['CHROM'].nunique() > 1:
                raise ValueError(f'More than one chromosome found in bgen {prefix}!')

        # 3. SNP List
        elif self._filtering_mode == FilteringMode.SNP_LIST:
            variant_index = self._query_snp_list(variant_index)

        return {'variant_index': variant_index, 'prefix': prefix, 'vars_found': len(variant_index) > 0}

    def _query_gene_list(self, variant_index: pd.DataFrame) -> pd.DataFrame:
        """Filter variants by a provided gene list

        :param variant_index: A pandas.DataFrame containing variants loaded from all provided chromosomes
        :return: A modified version of the variant_index pandas.DataFrame AFTER filtering on provided Gene list
        """

        # First filter to the gene Symbols / ENSTs we care about
        # Check if we are using ENST or SYMBOL:
        query_mode = None
        is_mixed = [0, 0]
        for name in self._gene_list:
            if name.startswith('ENST'):
                is_mixed[0] += 1
                query_mode = 'ENST'
            else:
                is_mixed[1] += 1
                query_mode = 'SYMBOL'
        if is_mixed[0] > 0 and is_mixed[1] > 1:
            raise ValueError('Mixed ENST and gene SYMBOLs. Please run gene filtering with only ONE type!')

        variant_index = variant_index[variant_index[query_mode].isin(self._gene_list)]

        # Catalogue genes (and counts) that we did find:
        symbol_counts = variant_index[query_mode].value_counts().to_dict()
        for symbol, count in symbol_counts.items():
            if symbol in self._found_gene_dict:
                self._found_gene_dict[symbol] += count
            else:
                self._found_gene_dict[symbol] = count

        # Now further filter down to the filtering_expression we are interested in
        variant_index = variant_index.query(self._filtering_expression)

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
        # Have to do this in two steps as the snps in the list may not be in the index, which breaks .loc:
        union_snps = self._snp_list.intersection(variant_index.index)
        variant_index = variant_index.loc[union_snps]

        # Catalogue SNPs that we did find:
        self._found_snp_set.update(variant_index.index)

        return variant_index

    def _make_gene_dict(self, variant_index: pd.DataFrame) -> Dict[str, pd.DataFrame]:
        """Convert variant indices into a pandas DataFrame and check totals for this filtering mask.

        :param variant_index: A pandas.DataFrame containing variants AFTER applying any filtering expression
        :return: The updated genes dictionary taken from :param genes:
        """

        # This is NOT thread safe! Do not add it to a threaded process!
        self._stat_dict.collect_stats(variant_index)

        variant_index = variant_index[['CHROM', 'POS', 'ENST']]
        variant_index = variant_index.sort_values(by='POS')  # Make sure sorted by position
        variant_index = variant_index.reset_index()  # Make sure we have a numerical index in POS order

        return variant_index
