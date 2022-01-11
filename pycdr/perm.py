import pandas as pd
from scipy.stats import rankdata
import numpy as np
import anndata as ad
import scipy
import tqdm
import logging

logger = logging.getLogger(__name__)

def create_rank_matrix(adata):
    arr = adata.X
    import scipy
    if scipy.sparse.issparse(arr):
        arr = arr.toarray()
    arrrank = rankdata(arr.T, axis = 0, method = "min")
    return arrrank

def permute_matrix(adata, arrrank, factor, nperm, genecol, seed = 42):
    
    ind = adata.var[genecol].isin(adata.uns["factor_loadings"][factor])

    num_perm = nperm
    
    genelength = arrrank.shape[0]

    #arrrank = rankdata(arr, axis = 0, method = "min")
    matreal = (arrrank[ind].mean(0)/genelength) - 0.5

    perm_list = []
    np.random.seed(seed)
    for i in range(num_perm):
        oops1 = arrrank.T
        oops = np.random.permutation(oops1).T
        matperm = (oops[ind].mean(0)/genelength) - 0.5
        perm_list.append(matperm)

    mato = np.vstack(perm_list)

    pmat = (num_perm - np.sum(matreal > mato, 0))/num_perm
    
    return pmat, matreal

def calculate_proportions(adata, cols, pmat, thresh = 0.05):
    ased = pd.DataFrame(adata.obs[cols])
    ased["phenovec"] = (pmat < thresh)
    truthcol = ased["phenovec"]
    ct = ased.groupby(cols).sum('phenovec')["phenovec"].tolist()
    nobs = ased.groupby(cols).size().tolist()
    nobs2 = ased.groupby(cols).size()
    from statsmodels.stats.proportion import proportions_chisquare
    pval_obj = proportions_chisquare(ct, nobs)
    return pval_obj, ct, nobs, nobs2, truthcol

def calculate_enrichment(adata, cols, factor_list, nperm, genecol, thresh, seed = 42):
    
    dict_res = {}
    dict_res_prop = {}
    
    Xarr = create_rank_matrix(adata)
    
    for i,j in tqdm.tqdm(enumerate(factor_list)):
        if i % 100 == 0:
            logger.info("num factors processed:: %s", i)
        pmat, matreal = permute_matrix(adata, Xarr, j, nperm, genecol, seed = seed)
        pval_obj, ct, nobs, df, truthcol = calculate_proportions(adata, cols, pmat, thresh = thresh)
        #i_ = "factor." + str(i)
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

    return dict_res_prop, vals_p

def get_df_loadings(adata):
    a_ = adata.uns["dict_res_prop"]
    b_ = adata.uns["pval_dict"]
    a = tidy_up_pheno(a_)
    b = tidy_up_dict(b_)
    whole = a.merge(b)
    whole = whole.dropna()
    ps = whole["pvalue"]
    from statsmodels.stats.multitest import fdrcorrection
    whole["fdr"] = fdrcorrection(ps)[1]

    return whole
    
def tidy_up_dict(dict_obj):
    oink = dict_obj
    df = pd.DataFrame.from_dict(oink, orient = "index")
    df.columns = ["statistic", "pvalue"]
    ps = df["pvalue"]
    df["factor_loading"] = df.index
    
    return df

def tidy_up_pheno(dict2_obj):
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

