import logging
import numpy as np
import scipy.sparse as ss
import pandas as pd
from scipy.stats import rankdata


def create_rank_matrix(X):
    """Create a gene-by-cell ranking matrix for ssGSEA scoring.

    Args:
        X (array): Expression matrix (cells x genes), can be sparse or dense.

    Returns:
        numpy.ndarray: Ranking matrix of shape (n_genes, n_cells).
    """
    arr = X
    if ss.issparse(arr):
        arr = arr.toarray()
    return rankdata(arr.T, axis=0, method="min")


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
    data = [zscore, floadings, pvals]
    columns = ["z_score", "Fs_diff", "pval"]
    fdr_mat = adata.uns.get("fdr_mat")
    if fdr_mat is not None:
        data.append(fdr_mat[:, i].tolist())
        columns.append("fdr")
    hum = pd.DataFrame(data).T
    hum.index = sigs
    hum.columns = columns
    return hum.sort_values("z_score", ascending = False)


def output_results(adata):
    """Extract and combine CDR-g results into a single DataFrame.

    Collects factor loading gene lists (ranked by z-score), gene counts,
    optional enrichment terms, and optional enrichment statistics into a
    combined table.

    Args:
        adata (anndata.AnnData): AnnData object after ``run_CDR_analysis``
            and optionally ``calculate_enrichment``.

    Returns:
        pandas.DataFrame: Combined results indexed by factor name with columns
            ``n_genes``, ``genes`` (list), ``top_genes`` (str), and optionally
            enrichment columns. Returns ``None`` if no CDR-g analysis found.
    """
    factor_loadings = adata.uns.get("factor_loadings")
    if factor_loadings is None:
        logging.getLogger(__name__).warning("No CDR-g analysis identified. Have you run the pipeline?")
        return None

    zscores = adata.uns.get("zscores")
    gene_names = adata.var_names.tolist()
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}

    rows = {}
    for fname, gene_list in factor_loadings.items():
        # Rank genes by z-score within this factor
        ranked_genes = list(gene_list)  # copy
        try:
            fi = int(fname.split(".")[-1])
        except (ValueError, IndexError):
            fi = None

        if zscores is not None and fi is not None and fi < zscores.shape[1] and len(ranked_genes) > 0:
            # Sort by z-score descending
            def _zscore(g):
                idx = gene_to_idx.get(g)
                if idx is not None and idx < zscores.shape[0]:
                    return zscores[idx, fi]
                return 0.0
            ranked_genes.sort(key=_zscore, reverse=True)

        top5 = ranked_genes[:5]
        rows[fname] = {
            "n_genes": len(ranked_genes),
            "genes": ranked_genes,
            "top_genes": ", ".join(top5),
        }

    df = pd.DataFrame.from_dict(rows, orient="index")

    # Add Fs-based dominant condition as fallback
    fs_dominant = adata.uns.get("dominant_condition")
    if fs_dominant is not None:
        # Build a mapping from factor name to dominant condition
        n_factors = len(fs_dominant)
        dc_map = {f"factor.{i}": fs_dominant[i] for i in range(n_factors)}
        df["dominant_condition"] = df.index.map(dc_map)

    # Join enrichment terms if available
    enrich_results = adata.uns.get("enrichment_results")
    if enrich_results is not None:
        terms_dict = {
            key: ", ".join(value)
            for key, value in enrich_results.items()
        }
        terms_df = pd.DataFrame.from_dict(
            terms_dict, orient="index", columns=["terms"]
        )
        df = pd.concat([df, terms_df], join="inner", axis=1)

    # Join enrichment stats if available
    enrich_stats = adata.uns.get("enrichment_stats")
    if enrich_stats is not None and isinstance(enrich_stats, pd.DataFrame):
        # If enrichment provides dominant_condition, drop Fs-based fallback first
        if "dominant_condition" in enrich_stats.columns and "dominant_condition" in df.columns:
            df = df.drop(columns=["dominant_condition"])
        df = pd.concat([df, enrich_stats], join="inner", axis=1)

    return df
