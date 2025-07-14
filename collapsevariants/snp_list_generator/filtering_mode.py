from enum import Enum, auto


class FilteringMode(Enum):
    """Enum to represent the mode of filtering used in this run of the applet."""
    GENE_LIST = auto()
    FILTERING_EXPRESSION = auto()
    SNP_LIST = auto()

