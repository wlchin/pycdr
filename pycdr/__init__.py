"""CDR-g: Condition-specific gene expression programs in single-cell data."""

__version__ = "0.0.3"

from .pycdr import run_CDR_analysis
from .experimental import calculate_enrichment, binarize_gset_on_adata
from .utils import (
    filter_genecounts_percent,
    filter_genecounts_numcells,
    output_results,
    get_top_genes,
)
