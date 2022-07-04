from scipy.stats import rankdata
import tqdm
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
import numpy as np


def create_rank_matrix(X):
    """create rank matrix for ssgsea

    Args:
        X (array): can be sparse or numpy array

    Returns:
        array: ranking matrix for all genes
    """
    arr = X
    import scipy
    if scipy.sparse.issparse(arr):
        arr = arr.toarray()
    arrrank = rankdata(arr.T, axis=0, method="min")
    return arrrank


def calculate_enrichment_single_geneset(geneset, arr_index, arrrank):
    """get a list of a single geneset 
    Args:
        gmtfile
    Returns:
        vector of gene set enrichment 
    """
    genelength = arrrank.shape[0]
    ind = arr_index.isin(geneset)
    matreal = (arrrank[ind].mean(0)/genelength) - 0.5

    return matreal


def calculate_enrichment_all_sets(dict_gene, arr_index, arrank):
    """perform enrichment on all the sets 
    Args:
        gmtfile or dictionary
        rank matrix
    Returns:
        mat of gene set enrichment 
    """
    
    res_dict = {}
  
    for i, geneset in tqdm.tqdm(dict_gene.items()):
        res = calculate_enrichment_single_geneset(geneset,  arr_index, arrank)
        res_dict[i] = res
    
    enrichment_df = pd.DataFrame.from_dict(res_dict)
    
    return enrichment_df


def calculate_kruskal_wallis(enrichment_ser, pheno_ser):
    """_summary_

    Args:
        enrichment_ser (_type_): _description_
        pheno_ser (_type_): _description_

    Returns:
        _type_: _description_
    """
    # break al values into list of lists
    # then use kruskal wallis for this
    testdf = pd.concat([enrichment_ser, pheno_ser], axis=1)
    testdf.columns = ["gene_val", "pheno"]
    hum = testdf.groupby('pheno')['gene_val'].apply(list).to_list()
    stat, pval = stats.kruskal(*hum)
    return stat, pval


def calculate_kruskal_wallis_all_sets(enrichment_df, pheno_ser):
    """_summary_

    Args:
        enrichment_df (_type_): _description_
        pheno_ser (_type_): _description_

    Returns:
        _type_: _description_
    """
    enrichment_df.index = pheno_ser.index
    factors = enrichment_df.columns
    
    results = {}
    for i in factors:
        enrichment_ser = enrichment_df[i]
        stat, pval = calculate_kruskal_wallis(enrichment_ser, pheno_ser)
        results[i] = [stat, pval]
    
    res = pd.DataFrame.from_dict(results, orient='index')
    res.columns = ["stat", "pval"]
    res["fdr"] = fdrcorrection(res["pval"])[1]
    res_sorted = res.sort_values(['stat', 'fdr'], ascending=[False, True])
        
    return res_sorted


def calculate_enrichment(adata, pheno):
    """_summary_

    Args:
        adata (_type_): _description_
        pheno (_type_): _description_

    Returns:
        _type_: _description_
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
    
    return results


def binarize_gset(arrrank, arr_index, geneset, nperm, threshold=0.1, seed=42):
    """calculate binary activation of gene set through permutation

    Calculates pvalues through comparing with permutation of matrices.
    Cells with active geneset have pvals below threshold

    Args:
        arrrank (_type_): _description_
        geneset (_type_): _description_
        nperm (_type_): _description_
        threshold (float, optional): _description_. Defaults to 01.
        seed (int, optional): _description_. Defaults to 42.

    Returns:
        _type_: _description_
    """
    perm_list = []
    np.random.seed(seed)

    matreal = calculate_enrichment_single_geneset(geneset, arr_index, arrrank)

    for i in range(nperm):
        arr_permuted = np.random.permutation(arrrank.T).T  # issue np permute
        matperm = calculate_enrichment_single_geneset(geneset, arr_index, arr_permuted)
        perm_list.append(matperm)

    mato = np.vstack(perm_list)  # allowing vectorisation
    pmat = (nperm - np.sum(matreal > mato, 0))/nperm
    active_cells = pmat < threshold
    
    return pmat, matreal, active_cells