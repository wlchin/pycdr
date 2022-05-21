import numpy as np
import scipy.sparse as ss
import pandas as pd
import anndata as ad


def filter_genecounts_percent(adata, cell_fraction, median_count_above_zero):
    """filter function for counts

    implements a gene count filter based on percentage of cells and median count,
    as per SCENIC.
    
    Args:
        adata: anndata object to be filtered
        pheno: phenotype to filter on
        percent_cells: the percent of cells which should contain the gene for total gene filtering
        small_pheno_frac: the fraction of the smallest phenotype containing the gene
        count_above_zero: count above the median that is used for total gene filtering
    
    Returns 
        adata: filtered anndata object
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
    """filters cells based on gene content

    Args:
        adata (anndata): anndata object
        count_threshold (int): number of counts as cutoff
        min_expressed_cells (int): desired cutoff for cells

    Returns:
        anndata: inplace modification of anndata
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
    """filter function for counts

    implements a gene count filter based on percentage of cells and median count,
    as per SCENIC.
    
    Args:
        adata: anndata object to be filtered
        pheno: phenotype to filter on
        percent_cells: the percent of cells which should contain the gene for total gene filtering
        small_pheno_frac: the fraction of the smallest phenotype containing the gene
        count_above_zero: count above the median that is used for total gene filtering
    
    Returns 
        adata: filtered anndata object
    """
    
    import pandas as pd
    sigs = adata.var.index.to_list()
    zscore = adata.uns["zscores"][:,i].tolist()
    floadings = adata.uns["Fs_diff"][:,i].tolist()
    pvals = adata.uns["pval_mat"][:,i].tolist()
    hum = pd.DataFrame([zscore, floadings, pvals]).T
    hum.index = sigs
    hum.columns = ["z_score", "Fs_diff", "pval"]
    return hum.sort_values("z_score", ascending = False)
