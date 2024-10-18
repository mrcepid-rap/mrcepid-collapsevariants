import gzip

import dxpy
import pandas as pd

from enum import Enum, auto
from pathlib import Path
from typing import Dict, TypedDict, List, Set, Optional

from collapsevariants.ingest_data import BGENIndex
from collapsevariants.collapse_logger import CollapseLOGGER
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.mrc_logger import MRCLogger


class GeneDict(TypedDict):
    """A TypedDict containing information about all variants collapsed into a given gene

    :cvar chrom: The chromosome this gene is on
    :cvar varIDs: A list of all variant IDs within this gene
    :cvar min_poss: The minimum coordinate for a variant within this gene. Is not a perfect proxy for gene start, but
        is good enough for the purposes of the code that needs that information.
    :cvar prefixes: A list of prefixes for the bgen files that contain variants for this gene
    """

    chrom: str
    varIDs: List[str]
    min_poss: int
    prefixes: List[str]


class FilteringMode(Enum):
    GENE_LIST = auto()
    FILTERING_EXPRESSION = auto()
    SNP_LIST = auto()


class StatDictionary:
    """A class to store statistics about the variants in a given run of this applet.

    This class was created to remove the 'stat tracking' functionality from the SNPListGenerator class and put it
    into its own class. This was done to make the SNPListGenerator class more readable and to separate functions.

    In brief, this class will be 'updated' as each variant index is processed by the SNPListGenerator class via the
    :func:`collect_stats` method. This method will update all attributes of this class based on the variant index.

    :param log_file: The CollapseLOGGER object to write stats to.
    """

    def __init__(self, log_file: CollapseLOGGER):

        self.log_file = log_file
        self._total_variants = 0
        self._max_missingness = 0
        self._max_af = 0
        self._max_cadd = 0
        self._max_revel = 0
        self._num_na_revel = 0
        self._num_pass = 0
        self._num_not_pass = 0
        self._loftee_counts = {}
        self._parsed_consequence_counts = {}
        self._vep_consequence_counts = {}
        self._collected_genes = []

    def collect_stats(self, variant_index: pd.DataFrame) -> None:
        """Method upates all the attributes of the class based on the variant index provided

        :param variant_index: A pandas.DataFrame containing variants from some file loaded by :func:`_load_variant_index`
        :return: None
        """

        self._increment_variant_count(len(variant_index))
        self._update_max_missingness(variant_index['F_MISSING'])
        self._update_max_af(variant_index['AF'])
        self._update_max_cadd(variant_index['CADD'])
        self._update_revel(variant_index['REVEL'])
        self._update_pass_stats(variant_index['FILTER'])
        self._update_loftee_counts(variant_index['LOFTEE'])
        self._update_parsed_consequence_counts(variant_index['parsed_csq'])
        self._update_vep_consequence_counts(variant_index['CSQ'])

    def write_bgen_stats(self) -> None:
        """Gets information relating to included variants in bgen format files (per-gene)

        This method just uses standard pandas.DataFrame methods to query the underlying table to get stats and report
        them in the applet log_file

        :return: None
        """

        self.log_file.write_header('Overall Statistics')
        self.log_file.write_int('Total number of variants', self._total_variants)
        self.log_file.write_float('Maximum missingness', self._max_missingness)
        self.log_file.write_scientific('Maximum Allele Frequency', self._max_af)
        self.log_file.write_float('Minimum CADD Score', self._max_cadd)
        self.log_file.write_float('Minimum REVEL Score', self._max_revel)
        self.log_file.write_int('Number of NA REVEL Scores', self._num_na_revel)
        self.log_file.write_int('Total number of PASS variants', self._num_pass)
        self.log_file.write_int('Total number of non-PASS variants', self._num_not_pass)

        # LOFTEE:
        for key, value in self._loftee_counts.items():
            key = f'Number of LOFTEE {key}'
            self.log_file.write_int(key, value)
        self.log_file.write_spacer()

        # Parsed Consequences:
        self.log_file.write_header('Consequence statistics')
        for key, value in self._parsed_consequence_counts.items():
            key = f'Number of Parsed Consequence – {key}'
            self.log_file.write_int(key, value)

        # VEP Consequences:
        for key, value in self._vep_consequence_counts.items():
            key = f'Number of VEP Consequence - {key}'
            self.log_file.write_int(key, value, is_vep=True)
        self.log_file.write_spacer()

    def check_list_filtering(self, filtering_mode: FilteringMode,
                             snp_list: Set[str], gene_list: Set[str],
                             found_snp_list: Set[str], found_gene_dict: Dict[str, int]) -> None:
        """This method will check against either the original SNP or GENE list to determine if variants / genes were
        found in the final set of filtered variants.

        This function is a bit strange in that it is taking final information from the SNPListGenerator class and
        using it to determine if genes / snps were found in the final set of filtered variants. This is done to
        enable parallel processing of the variant indexes and to keep the SNPListGenerator class clean.

        :param filtering_mode: The mode of filtering used in this run
        :param snp_list: A set of SNPs to filter against
        :param gene_list: A set of genes to filter against
        :param found_snp_list: A set of SNPs found in the final set of filtered variants
        :param found_gene_dict: A dictionary of genes and counts of those genes found in the final set of filtered
            variants
        """

        ## For SNP list
        if filtering_mode == FilteringMode.SNP_LIST:

            # Then use the boolean to subset the snp array
            missing_snps = snp_list.difference(found_snp_list)
            found_snps = snp_list.intersection(found_snp_list)

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

        ## For Gene list
        elif filtering_mode == FilteringMode.GENE_LIST:
            # Following code is just to report stats about genes that were/were not found
            missing_genes = gene_list.difference(found_gene_dict.keys())

            self.log_file.write_header('Per-gene stats')
            for gene, count in found_gene_dict.items():
                self.log_file.write_int(gene, count)
            self.log_file.write_spacer()

            self.log_file.write_generic('Variants not Found:')
            for gene in missing_genes:
                self.log_file.write_generic(gene)

        else:
            # Do nothing if we are not filtering on a list
            pass

    def get_total_sites(self) -> int:
        """Getters for the total number of sites after filtering

        :return: The total number of sites after filtering
        """
        return self._total_variants

    def _increment_variant_count(self, num_variants: int) -> None:
        """Increment the total number of variants by the number of variants in the current variant index

        :param num_variants: Number of variants in the current variant index
        :return: None
        """
        self._total_variants += num_variants

    def _update_max_missingness(self, missingness_list: pd.Series) -> None:
        """Use the max function to update the maximum missingness for this run

        :param missingness_list: A pandas.Series containing the missingness column from the variant index
        :return: None
        """
        self._max_missingness = max(self._max_missingness, missingness_list.max() * 100)

    def _update_max_af(self, af_list: pd.Series) -> None:
        """Use the max function to update the maximum allele frequency for this run

        :param af_list: A pandas.Series containing the allele frequency column from the variant index
        :return: None
        """
        self._max_af = max(self._max_af, af_list.max(skipna=True))

    def _update_max_cadd(self, cadd_list: pd.Series) -> None:
        """Use the max function to update the maximum CADD score for this run

        :param cadd_list: A pandas.Series containing the CADD column from the variant index
        :return: None
        """
        self._max_cadd = max(self._max_cadd, cadd_list.max(skipna=True))

    def _update_revel(self, revel_list: pd.Series) -> None:
        """Update REVEL-based stats for this run

        Note: This method updates MULTIPLE stats at once (REVEL max and number of NA values)

        :param revel_list: A pandas.Series containing the REVEL column from the variant index
        :return: None
        """

        self._max_revel = max(self._max_revel, revel_list.max(skipna=True))
        num_na = revel_list.isna().sum()
        self._num_na_revel += num_na

    def _update_pass_stats(self, filter_list: pd.Series) -> None:
        """Update PASS / FAIL counts for this run.

        This method updates MULTIPLE stats at once (PASS / FAIL counts).

        :param filter_list: A pandas.Series containing the FILTER column from the variant index
        :return: None
        """

        filter_counts = filter_list.value_counts()
        self._num_pass += filter_counts.get('PASS', 0)
        self._num_not_pass += filter_counts.get('FAIL', 0)

    def _update_loftee_counts(self, loftee_counts: pd.Series) -> None:
        """Count the number of variants per LOFTEE filter (e.g., HC / LC)

        :param loftee_counts: A pandas.Series containing the LOFTEE column from the variant index
        :return: None
        """

        for loftee_filter, count in loftee_counts.value_counts().items():
            if loftee_filter in self._loftee_counts:
                self._loftee_counts[loftee_filter] += count
            else:
                self._loftee_counts[loftee_filter] = count

    def _update_parsed_consequence_counts(self, parsed_consequence_counts: pd.Series) -> None:
        """Update the parsed consequence counts for this run.

        The method will iterate through the parsed_consequence_counts Series and update the counts for each consequence

        :param parsed_consequence_counts: A pandas.Series containing the parsed consequence column from the variant index
        :return: None
        """

        for parsed_csq, count in parsed_consequence_counts.value_counts().items():
            if parsed_csq in self._parsed_consequence_counts:
                self._parsed_consequence_counts[parsed_csq] += count
            else:
                self._parsed_consequence_counts[parsed_csq] = count

    def _update_vep_consequence_counts(self, vep_consequence_counts: pd.Series) -> None:
        for vep_csq, count in vep_consequence_counts.value_counts().items():
            if vep_csq in self._vep_consequence_counts:
                self._vep_consequence_counts[vep_csq] += count
            else:
                self._vep_consequence_counts[vep_csq] = count


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

    :param bgen_index: A dict that contains the bgen file prefixes as keys and a BGENIndex TypedDict containing
        information on filepaths containing filtered & annotated variants.
    :param filtering_expression: A pandas.query() compatible expression to filter variants on, 'None' if not provided.
    :param gene_list_path: A path to a list of genes to filter against.
    :param snp_list_path: A path to a list of SNPs to filter against.
    :param log_file: A class of CollapseLOGGER to store information on variant filtering
    """

    def __init__(self, bgen_index: Dict[str, BGENIndex], filtering_expression: str,
                 gene_list_path: Optional[Path], snp_list_path: Optional[Path], log_file: CollapseLOGGER):

        self._logger = MRCLogger(__name__).get_logger()

        # Set defaults for all class variables
        self._bgen_index = bgen_index
        self._filtering_expression = filtering_expression
        self._gene_list_path = gene_list_path  # This is a boolean if we found a gene list
        self._gene_list = set()  # Placeholder for a numpy array of gene symbols
        self._found_gene_dict = dict()
        self._snp_list_path = snp_list_path  # This is a boolean if we found a SNP list
        self._snp_list = set()  # Placeholder for a numpy array of SNP IDs
        self._found_snp_set = set()  # Placeholder for a numpy array of SNP IDs

        # stat_dict is a dictionary that will be used to store the stats for each variant index parsed
        self._stat_dict = StatDictionary(log_file)

        # Decide on the type of filtering we are doing (SNP, Gene, or expression)
        self._filtering_mode = self._decide_filtering_mode()

        # Iterate through all possible bgens in parallel and filter them
        thread_utility = ThreadUtility()
        for prefix, data_dict in self._bgen_index.items():
            thread_utility.launch_job(self._query_variant_index,
                                      vep_id=data_dict['vep'],
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

        if self._filtering_expression is not None and self._gene_list_path:

            with self._gene_list_path.open('r') as my_genelist_file:
                gene_list = my_genelist_file.readlines()
                self._gene_list = set([gene.rstrip() for gene in gene_list])

            return FilteringMode.GENE_LIST

        # 2. Filtering expression
        elif self._filtering_expression is not None:
            return FilteringMode.FILTERING_EXPRESSION

        # 3. SNP List
        elif self._snp_list_path:

            with self._snp_list_path.open('r') as my_snplist_file:
                snp_list = my_snplist_file.readlines()
                self._snp_list = set([snp.rstrip() for snp in snp_list])

            return FilteringMode.SNP_LIST

        else:
            raise ValueError('No filtering expression, SNP list, or Gene list provided!')

    @staticmethod
    def _load_variant_index(vep_id: str) -> pd.DataFrame:
        """Load vep annotated *.tsv.gz files into a pandas DataFrame

        Note that this method uses the dxpy.open_dxfile method, which streams the file from the DNAnexus platform
        rather than downloading it directly to the local filesystem.

        :param vep_id: DXFile ID for the vep index to load
        :return: A pandas.DataFrame containing variants loaded from all provided chromosomes
        """

        current_vep = pd.read_csv(gzip.open(dxpy.open_dxfile(vep_id, mode='rb'), mode='rt'), sep="\t",
                                  index_col='varID',
                                  dtype={'SIFT': str, 'POLYPHEN': str})

        return current_vep

    # Actually run the query and modify the variant_index loaded above based on this/these queries
    def _query_variant_index(self, vep_id: str, prefix: str) -> dict:
        """Query self.variant_index based on the logic outlined in the __init__ method for this class

        As described in the __init__ for the class, there are three possibilities, each with an if/ifelse, and in this
        logical order:

        1. Filtering Expression + Gene List = Select specific genes to collapse on
        2. Filtering Expression = Any variant which qualifies under the given filtering expression
        3. SNPList = Select specific SNPs to collapse on

        :param vep_id: DXFile ID for the vep index to load. Handled  by :func:`_load_variant_index`
        :param prefix: The prefix of the bgen file loaded
        :return: A dict containing: 1) pandas.DataFrame containing variants from :param variant_index: filtered based on provided
            input parameters 2) The prefix of the file loaded 3) The chromosome of all variants found 4) A boolean
            indicating if any variants were found after filtering 5) The minimum position of the entire bgen file.
        """

        variant_index = self._load_variant_index(vep_id)

        # 1. Filtering expression + Gene List
        if self._filtering_mode == FilteringMode.GENE_LIST:
            variant_index = self._query_gene_list(variant_index)
            poss_chrom = 'GENE'

        # 2. Filtering expression
        elif self._filtering_mode == FilteringMode.FILTERING_EXPRESSION:
            variant_index = self._query_filtering_expression(variant_index)
            poss_chrom = variant_index['CHROM'].unique()
            if len(poss_chrom) > 1:
                raise ValueError(f'More than one chromosome found in bgen {prefix}!')
            else:
                poss_chrom = poss_chrom[0]

        # 3. SNP List
        elif self._filtering_mode == FilteringMode.SNP_LIST:
            variant_index = self._query_snp_list(variant_index)
            poss_chrom = 'SNP'

        return {'variant_index': variant_index, 'prefix': prefix, 'chrom': poss_chrom,
                'vars_found': len(variant_index) > 0, 'min_pos': min_pos}

    def _query_gene_list(self, variant_index: pd.DataFrame) -> pd.DataFrame:
        """Filter variants by a provided gene list

        :param variant_index: A pandas.DataFrame containing variants loaded from all provided chromosomes
        :return: A modified version of the variant_index pandas.DataFrame AFTER filtering on provided Gene list
        """

        # First filter to the gene Symbols we care about
        variant_index = variant_index[variant_index['SYMBOL'].isin(self._gene_list)]

        # Catalogue genes (and counts) that we did find:
        symbol_counts = variant_index['SYMBOL'].value_counts().to_dict()
        for symbol, count in symbol_counts.items():
            if symbol in self._found_gene_dict:
                self._found_gene_dict[symbol] += count
            else:
                self._found_gene_dict[symbol] = count

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

        # Catalogue SNPs that we did find:
        self._found_snp_set.update(variant_index.index)

        # Here we create a fake 'GENE' to collapse on later to make the process of generating a mask easier.
        # We also use this gene ID to notify runassociationtesting that we are doing a SNP-based mask
        variant_index['ENST'] = 'ENST00000000000'

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
        variant_index.reset_index(inplace=True)

        return variant_index
