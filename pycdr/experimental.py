from scipy.stats import rankdata
import tqdm
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
   
def create_rank_matrix(X):
    """create ranking matrix
    Args:
        adata (anndata): anndata of interest
    Returns:
        arrank: ranking matrix, stored in layers section
    """
    arr = X
    import scipy
    if scipy.sparse.issparse(arr):
        arr = arr.toarray()
    arrrank = rankdata(arr.T, axis = 0, method = "min")
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
  
    for i,geneset in tqdm.tqdm(dict_gene.items()):
        res = calculate_enrichment_single_geneset(geneset,  arr_index, arrank)
        res_dict[i] = res
    
    enrichment_df = pd.DataFrame.from_dict(res_dict)
    
    return enrichment_df

def calculate_kruskal_wallis(enrichment_ser, pheno_ser):
    """calculate for a single df 
    Args:
        gmtfile
    Returns:
        vector of gene set enrichment 
    """
    # break al values into list of lists
    # then use kruskal wallis for this
    testdf = pd.concat([enrichment_ser, pheno_ser], axis = 1)
    testdf.columns = ["gene_val", "pheno"]
    hum = testdf.groupby('pheno')['gene_val'].apply(list).to_list()
    stat, pval = stats.kruskal(*hum)
    return stat, pval

def calculate_kruskal_wallis_all_sets(enrichment_df, pheno_ser):
    enrichment_df.index = pheno_ser.index
    factors = enrichment_df.columns
    
    results = {}
    for i in factors:
        enrichment_ser = enrichment_df[i]
        stat, pval = calculate_kruskal_wallis(enrichment_ser, pheno_ser)
        results[i] = [stat, pval]
    
    res = pd.DataFrame.from_dict(results, orient = 'index')
    res.columns = ["stat", "pval"]
    res["fdr"] = fdrcorrection(res["pval"])[1]
    res_sorted = res.sort_values(['stat', 'fdr'], ascending=[False, True])
        
    return res_sorted

def calculate_enrichment(adata, pheno):
    """Driver function for anndata object
    
       Will store the ranking matrix in layers
    
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