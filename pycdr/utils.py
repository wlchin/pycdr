import logging
import numpy as np
import scipy.sparse as ss
import pandas as pd


def filter_genecounts_percent(adata, cell_fraction, median_count_above_zero):
    """Filter genes based on total expression relative to cell fraction and median count.

    Implements a SCENIC-style gene count filter. Genes whose total count
    falls below ``n_cells * cell_fraction * (median_nonzero + median_count_above_zero)``
    are removed.

    Args:
        adata (anndata.AnnData): AnnData object to filter.
        cell_fraction (float): Fraction of cells used to compute the count threshold.
        median_count_above_zero (float): Value added to the median of nonzero counts
            when computing the threshold.

    Returns:
        anndata.AnnData: Filtered AnnData object (subset of genes).
    """
    if ss.issparse(adata.X):
        matdense = adata.X.toarray()
    else:
        matdense = adata.X

    abovezero = matdense[matdense > 0]
    thresh = np.median(abovezero) + median_count_above_zero
    total_gene_count_thresh = np.round(matdense.shape[0] * cell_fraction * thresh)
    adata.uns["total_gene_thresh"] = total_gene_count_thresh
    adata = adata[:,(matdense.sum(0) > total_gene_count_thresh)]

    return adata

def filter_genecounts_numcells(adata, count_threshold, min_expressed_cells):
    """Filter genes by the number of cells expressing them above a count threshold.

    Args:
        adata (anndata.AnnData): AnnData object.
        count_threshold (int): Minimum count for a gene to be considered expressed in a cell.
        min_expressed_cells (int): Minimum number of cells that must express the gene.

    Returns:
        anndata.AnnData: Filtered AnnData object (subset of genes).
    """

    num_cells_thresh = min_expressed_cells

    if ss.issparse(adata.X):
        matdense = adata.X.toarray()
    else:
        matdense = adata.X

    num_cells_filter_indices = (np.greater(matdense, count_threshold).sum(0) > num_cells_thresh)

    adata = adata[:,num_cells_filter_indices]
    adata.uns["num_cells_thresh"] = num_cells_thresh

    return adata

def get_top_genes(adata, i):
    """Return a ranked table of gene statistics for a given factor loading.

    Extracts z-scores, factor loading differences, and permutation p-values
    for all genes in the specified factor, sorted by descending z-score.

    Args:
        adata (anndata.AnnData): AnnData object after ``run_CDR_analysis``.
        i (int): Factor index (column index into ``adata.uns["zscores"]``).

    Returns:
        pandas.DataFrame: Table with columns ``z_score``, ``Fs_diff``, ``pval``,
            indexed by gene name and sorted by descending z-score.
    """

    sigs = adata.var.index.to_list()
    zscore = adata.uns["zscores"][:,i].tolist()
    floadings = adata.uns["Fs_diff"][:,i].tolist()
    pvals = adata.uns["pval_mat"][:,i].tolist()
    hum = pd.DataFrame([zscore, floadings, pvals]).T
    hum.index = sigs
    hum.columns = ["z_score", "Fs_diff", "pval"]
    return hum.sort_values("z_score", ascending = False)


def output_results(adata):
    """Extract and combine CDR-g results into a single DataFrame.

    Collects factor loading gene lists, optional enrichment terms, and
    optional enrichment statistics (from ``adata.uns["enrichment_stats"]``)
    into a combined table.

    Args:
        adata (anndata.AnnData): AnnData object after ``run_CDR_analysis``
            and optionally ``calculate_enrichment``.

    Returns:
        pandas.DataFrame: Combined results indexed by factor name, or ``None``
            if no CDR-g analysis results are found.
    """
    enrichment = False
    stats = False

    try:
        dict_variable = {key:",".join(value) for (key,value) in adata.uns["factor_loadings"].items()}
        genes = pd.DataFrame.from_dict(dict_variable, orient='index', columns=['genes'])
    except KeyError:
        logging.getLogger(__name__).warning("No CDR-g analysis identified. Have you run the pipeline?")
        return None

    try:
        dict_variable = {key:",".join(value) for (key,value) in adata.uns["enrichment_results"].items()}
        terms = pd.DataFrame.from_dict(dict_variable, orient='index', columns=['terms'])
        enrichment = True
    except KeyError:
        logging.getLogger(__name__).info("No enrichment results provided. Run enrichment_utils if required.")

    try:
        df_stats = adata.uns["enrichment_stats"]
        stats = True
    except KeyError:
        logging.getLogger(__name__).info("No enrichment stats calculated. Run if required.")

    if not enrichment and stats:
        df = pd.concat([genes, df_stats], join='inner', axis = 1)

    if enrichment and not stats:
        df = pd.concat([genes, terms], join='inner', axis = 1)

    if not enrichment and not stats:
        df = genes

    if enrichment and stats:
        df = pd.concat([genes, terms, df_stats], join='inner', axis = 1)

    return df
