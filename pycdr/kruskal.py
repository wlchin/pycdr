import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
import tqdm

from .utils import create_rank_matrix


def calculate_enrichment_single_geneset(geneset, arr_index, arrrank):
    """Compute ssGSEA enrichment scores for a single gene set.

    Args:
        geneset (list): Gene names belonging to the gene set.
        arr_index (pandas.Index): Gene names corresponding to rows of arrrank.
        arrrank (numpy.ndarray): Ranking matrix of shape (n_genes, n_cells).

    Returns:
        numpy.ndarray: Enrichment score vector of length n_cells.
    """
    genelength = arrrank.shape[0]
    ind = arr_index.isin(geneset)
    matreal = (arrrank[ind].mean(0)/genelength) - 0.5

    return matreal


def calculate_enrichment_all_sets(dict_gene, arr_index, arrank):
    """Compute ssGSEA enrichment scores for all gene sets.

    Args:
        dict_gene (dict): Mapping of gene set names to gene lists.
        arr_index (pandas.Index): Gene names corresponding to rows of arrank.
        arrank (numpy.ndarray): Ranking matrix of shape (n_genes, n_cells).

    Returns:
        pandas.DataFrame: Enrichment scores with cells as rows and gene sets as columns.
    """

    res_dict = {}

    for i, geneset in tqdm.tqdm(dict_gene.items()):
        res = calculate_enrichment_single_geneset(geneset,  arr_index, arrank)
        res_dict[i] = res

    enrichment_df = pd.DataFrame.from_dict(res_dict)

    return enrichment_df


def calculate_kruskal_wallis(enrichment_ser, pheno_ser):
    """Test whether enrichment scores differ across phenotype groups using Kruskal-Wallis.

    Args:
        enrichment_ser (pandas.Series): Enrichment scores for one gene set.
        pheno_ser (pandas.Series): Phenotype labels for each cell.

    Returns:
        tuple: (H-statistic, p-value) from the Kruskal-Wallis test.
    """
    # break al values into list of lists
    # then use kruskal wallis for this
    testdf = pd.concat([enrichment_ser, pheno_ser], axis=1)
    testdf.columns = ["gene_val", "pheno"]
    hum = testdf.groupby('pheno')['gene_val'].apply(list).to_list()
    stat, pval = stats.kruskal(*hum)
    return stat, pval


def calculate_kruskal_wallis_all_sets(enrichment_df, pheno_ser):
    """Run Kruskal-Wallis tests across all gene sets with FDR correction.

    Args:
        enrichment_df (pandas.DataFrame): Enrichment score matrix (cells x gene sets).
        pheno_ser (pandas.Series): Phenotype labels for each cell.

    Returns:
        pandas.DataFrame: Results with columns ``stat``, ``pval``, ``fdr``,
            sorted by descending statistic and ascending FDR.
    """
    enrichment_df.index = pheno_ser.index
    factors = enrichment_df.columns

    results = {}
    for i in factors:
        enrichment_ser = enrichment_df[i]
        stat, pval = calculate_kruskal_wallis(enrichment_ser, pheno_ser)
        medians = enrichment_ser.groupby(pheno_ser).median()
        dominant = str(medians.idxmax()) if medians.notna().any() else ""
        results[i] = [stat, pval, dominant]

    res = pd.DataFrame.from_dict(results, orient='index')
    res.columns = ["stat", "pval", "dominant_condition"]
    res["fdr"] = fdrcorrection(res["pval"])[1]
    res_sorted = res.sort_values(['stat', 'fdr'], ascending=[False, True])

    return res_sorted


def calculate_enrichment(adata, pheno):
    """Run Kruskal-Wallis enrichment analysis on CDR-g factor loadings.

    Computes ssGSEA enrichment scores for each factor loading gene set,
    then tests for differential enrichment across phenotype groups using
    a Kruskal-Wallis test with FDR correction.

    Args:
        adata (anndata.AnnData): AnnData object after ``run_CDR_analysis``.
            Must contain ``adata.uns["factor_loadings"]``.
        pheno (str): Column name in ``adata.obs`` identifying the condition of interest.

    Returns:
        pandas.DataFrame: Results with columns ``stat``, ``pval``, ``fdr``,
            indexed by factor name and sorted by descending statistic.

    Side effects:
        Stores the following keys on *adata*:

        - ``adata.layers["rank_matrix"]``: gene-by-cell ranking matrix (transposed).
        - ``adata.obsm["enrichment_score_matrix"]``: enrichment score matrix.
        - ``adata.uns["enrichment_stats"]``: the returned DataFrame (for ``output_results``).
    """
    dict_gene = adata.uns["factor_loadings"]
    arr_index = adata.var.index
    pheno_ser = adata.obs[pheno]

    arrank = create_rank_matrix(adata.X)
    enrichment_df = calculate_enrichment_all_sets(dict_gene, arr_index, arrank)
    results = calculate_kruskal_wallis_all_sets(enrichment_df, pheno_ser)

    # storage
    adata.layers['rank_matrix'] = arrank.T
    adata.obsm["enrichment_score_matrix"] = enrichment_df.to_numpy()
    adata.uns["enrichment_stats"] = results

    return results


def binarize_gset(arrrank, arr_index, geneset, nperm = 50, threshold=0.1, seed=42):
    """Determine binary gene set activation via permutation testing.

    Compares each cell's enrichment score against a null distribution
    built from *nperm* row-permutations of the ranking matrix. Cells
    whose score exceeds the null at the given threshold are considered
    active.

    Args:
        arrrank (numpy.ndarray): Ranking matrix of shape (n_genes, n_cells).
        arr_index (pandas.Index): Gene names corresponding to rows of arrrank.
        geneset (list): Gene names belonging to the gene set.
        nperm (int, optional): Number of permutations. Defaults to 50.
        threshold (float, optional): P-value cutoff for activation. Defaults to 0.1.
        seed (int, optional): Random seed for reproducibility. Defaults to 42.

    Returns:
        tuple: (pmat, matreal, active_cells) where *pmat* is the permutation
            p-value array, *matreal* is the observed enrichment score, and
            *active_cells* is a boolean array indicating activation.
    """
    rng = np.random.default_rng(seed)

    matreal = calculate_enrichment_single_geneset(geneset, arr_index, arrrank)

    genelength = arrrank.shape[0]
    mato = np.empty((nperm, arrrank.shape[1]))
    for i in range(nperm):
        perm_idx = rng.permutation(genelength)
        mato[i] = calculate_enrichment_single_geneset(geneset, arr_index, arrrank[perm_idx])
    pmat = (nperm - np.sum(matreal > mato, 0))/nperm
    active_cells = pmat < threshold

    return pmat, matreal, active_cells

def binarize_gset_on_adata(adata, factor_list, **kwargs):
    """Compute binary gene set activation for each factor on an AnnData object.

    Args:
        adata (anndata.AnnData): AnnData object with ``factor_loadings`` in ``.uns``.
        factor_list (list): Factor keys to process from ``adata.uns["factor_loadings"]``.
        **kwargs: Passed to :func:`binarize_gset` (e.g. *nperm*, *threshold*, *seed*).

    Returns:
        numpy.ndarray: Boolean array of shape (n_cells, n_factors).
    """
    arrrank = create_rank_matrix(adata.X)
    arr_index = adata.var.index
    agg = []
    for factor in factor_list:
        geneset = adata.uns["factor_loadings"][factor]
        _, _, active_cells = binarize_gset(arrrank, arr_index, geneset, **kwargs)
        agg.append(active_cells)

    return np.column_stack(agg)
