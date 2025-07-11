from typing import Set, Dict

import pandas as pd

from snp_list_generator.filtering_mode import FilteringMode
from utilities.collapse_logger import CollapseLOGGER


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
        self._update_pass_stats(variant_index['FILTER'])
        self._update_loftee_counts(variant_index['LOFTEE'])
        self._update_parsed_consequence_counts(variant_index['PARSED_CSQ'])
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
                raise ValueError('No SNPs remain after using SNPlist!')

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
