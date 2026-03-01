import logging

import numpy as np
import pandas as pd
from statsmodels.stats.proportion import proportions_chisquare
from statsmodels.stats.multitest import fdrcorrection
import tqdm

from .utils import create_rank_matrix as _create_rank_matrix

logger = logging.getLogger(__name__)

def create_rank_matrix(adata):
    """create ranking matrix

    Args:
        adata (anndata): anndata of interest

    Returns:
        arrank: ranking matrix
    """
    return _create_rank_matrix(adata.X)

def permute_matrix(adata, arrrank, factor, nperm, genecol, seed = 42):
    """Compute observed ssGSEA score and a permutation null for one factor.

    Args:
        adata (anndata.AnnData): AnnData object with ``factor_loadings`` in ``.uns``.
        arrrank (numpy.ndarray): Ranking matrix of shape (n_genes, n_cells).
        factor (str): Key into ``adata.uns["factor_loadings"]``.
        nperm (int): Number of permutations.
        genecol (str): Column in ``adata.var`` containing gene names.
        seed (int, optional): Random seed. Defaults to 42.

    Returns:
        tuple: (pmat, matreal) — permutation p-values and observed scores per cell.
    """
    ind = adata.var[genecol].isin(adata.uns["factor_loadings"][factor])

    genelength = arrrank.shape[0]

    matreal = (arrrank[ind].mean(0)/genelength) - 0.5

    rng = np.random.default_rng(seed)
    # Batch: generate all permutation indices at once
    all_perm_idx = np.array([rng.permutation(genelength) for _ in range(nperm)])
    # (nperm, genelength, n_cells) -> select rows matching ind -> mean
    all_permuted = arrrank[all_perm_idx][:, ind, :]
    mato = (all_permuted.mean(axis=1) / genelength) - 0.5

    pmat = (nperm - np.sum(matreal > mato, 0))/nperm

    return pmat, matreal

def calculate_proportions(adata, cols, pmat, thresh = 0.05):
    """Test whether the proportion of active cells differs across phenotype groups.

    Args:
        adata (anndata.AnnData): AnnData object.
        cols (str): Phenotype column in ``adata.obs``.
        pmat (numpy.ndarray): Permutation p-values per cell.
        thresh (float, optional): Activation threshold. Defaults to 0.05.

    Returns:
        tuple: (pval_obj, ct, nobs, nobs2, truthcol) — chi-square result,
            counts of active cells, total cells per group, group sizes as
            Series, and the boolean activation column.
    """
    ased = pd.DataFrame(adata.obs[cols])
    ased["phenovec"] = (pmat < thresh)
    truthcol = ased["phenovec"]
    ct = ased.groupby(cols).sum('phenovec')["phenovec"].tolist()
    nobs = ased.groupby(cols).size().tolist()
    nobs2 = ased.groupby(cols).size()
    pval_obj = proportions_chisquare(ct, nobs)
    return pval_obj, ct, nobs, nobs2, truthcol

def calculate_enrichment(adata, cols, factor_list, nperm, genecol, thresh, seed = 42):
    """Peforms enrichment test on factor loadings

        For each factor loadings, this function performs ssGSEA 
        (single sample gsea) on each assesses statistical significance
        using a test of proportions 

    Args:
        adata (AnnData): anndata object in question
        cols (str): phenotype for testing
        factor_list (list): list of factors in the uns for testing
        nperm (int): number of permutations for ssGSEA
        genecol (str): name of column for genes
        thresh (float): threshold for "active" gene set
        seed (int, optional): for reproducibility. Defaults to 42.

    Returns:
        tuple: (dict_res_prop, vals_p) — proportions dictionary and p-value
            dictionary. Results are also stored in ``adata.uns`` and ``adata.obs``.

    .. note::
        This is the legacy permutation-based enrichment function. For a
        simpler Kruskal-Wallis approach, use
        :func:`pycdr.kruskal.calculate_enrichment`.
    """
    dict_res = {}
    dict_res_prop = {}
    
    Xarr = create_rank_matrix(adata)
    
    for i,j in tqdm.tqdm(enumerate(factor_list)):
        if i % 100 == 0:
            logger.info("num factors processed:: %s", i)
        pmat, matreal = permute_matrix(adata, Xarr, j, nperm, genecol, seed = seed)
        pval_obj, ct, nobs, df, truthcol = calculate_proportions(adata, cols, pmat, thresh = thresh)
        i_ = str(j)
        ii_ = str(j) + "_score"
        dict_res[i_] = pval_obj
        dict_res_prop[i_] = [ct, nobs, df.index.to_list()] # changed index
        adata.obs[i_] = truthcol
        adata.obs[ii_] = matreal
        adata.obs[i_] = adata.obs[i_].astype(int).astype("category")
    
    
    adata.uns["dict_res_prop"] = dict_res_prop
    vals_p = {i:list(j[0:2]) for (i,j) in dict_res.items()} # this has problems
    adata.uns["pval_dict"] = vals_p
    adata.uns["enrichment_stats"] = get_df_loadings(adata).set_index("factor_loading")

    return dict_res_prop, vals_p

def get_df_loadings(adata):
    """Combine proportion and p-value results into a single DataFrame.

    Requires ``adata.uns["dict_res_prop"]`` and ``adata.uns["pval_dict"]``
    from :func:`calculate_enrichment`.

    Args:
        adata (anndata.AnnData): AnnData object after legacy ``calculate_enrichment``.

    Returns:
        pandas.DataFrame: Merged results with FDR-corrected p-values.

    .. note::
        Legacy helper — only needed when using the permutation-based enrichment.
    """
    a_ = adata.uns["dict_res_prop"]
    b_ = adata.uns["pval_dict"]
    a = tidy_up_pheno(a_)
    b = tidy_up_dict(b_)
    whole = a.merge(b)
    whole = whole.dropna()
    whole = whole.rename(columns={"max_P": "dominant_condition"})
    ps = whole["pvalue"]
    whole["fdr"] = fdrcorrection(ps)[1]

    return whole
    
def tidy_up_dict(dict_obj):
    """Convert the p-value dictionary into a DataFrame with statistic and pvalue columns.

    Args:
        dict_obj (dict): Mapping of factor names to [statistic, pvalue] lists.

    Returns:
        pandas.DataFrame: Table with columns ``statistic``, ``pvalue``, ``factor_loading``.
    """
    output = dict_obj
    df = pd.DataFrame.from_dict(output, orient = "index")
    df.columns = ["statistic", "pvalue"]
    ps = df["pvalue"]
    df["factor_loading"] = df.index
    
    return df

def tidy_up_pheno(dict2_obj):
    """Convert the proportions dictionary into a DataFrame of activation statistics.

    Args:
        dict2_obj (dict): Mapping of factor names to [counts, totals, labels] lists.

    Returns:
        pandas.DataFrame: Table with columns ``max_P``, ``a_max``, ``a_range``,
            ``a_mean``, ``factor_loading``.
    """
    dfd = {}
    
    for i,j in dict2_obj.items():
        a = np.array(j[0]).astype("int")
        b = np.array(j[1]).astype("int")
        props = np.true_divide(a, b)
        maxpheno = j[2][np.argmax(props)]
        diff = props[np.argmax(props)] - props[np.argmin(props)]
        max_a = props[np.argmax(props)]
        av_activation = np.mean(props)
        dfd[i] = (maxpheno, max_a, diff, av_activation)
    df = pd.DataFrame.from_dict(dfd, orient = "index")
    df.columns = ["max_P", "a_max", "a_range","a_mean"]
    df["factor_loading"] = df.index
        
    return df

