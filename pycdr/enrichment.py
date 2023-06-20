import scanpy as sc
import decoupler as dc
import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, chi2
import gseapy as gp

def parse_values(dict_key, dict_value):
    if len(dict_value) > 0:
        df_vals = pd.DataFrame(dict_value)
        df_vals["factor_loading"] = dict_key
        df_vals.columns = ["gene", "factor_loading"]
    else:
        df_vals = pd.DataFrame(columns = ["gene", "factor_loading"])
    return df_vals

def flatten_factor_loading_matrix(adata):
    factor_loading_dict = adata.uns["factor_loadings"]
    df_genes = pd.concat([parse_values(x,y) for x,y in factor_loading_dict.items()], axis = 0)
    return df_genes

def transfer_dict(CD8TEM, adata, X_new):
    old_dict = CD8TEM.uns["factor_loadings"]
    your_keys = X_new
    dict_you_want = {key: old_dict[key] for key in your_keys}
    adata.uns["factor_loadings"] = dict_you_want
    return adata

def select_important_fetures(adata, target_column, num_features = 5):
    
    test = flatten_factor_loading_matrix(adata)

    dc.run_ora(
        mat=adata,
        net=test,
        source='factor_loading',
        target='gene',
        min_n=3,
        verbose=True,
        use_raw=False
    )

    y = adata.obs[target_column]
    X = adata.obsm['ora_estimate'].astype('float32')

    max_e = np.nanmax(X[np.isfinite(X)])
    X[~np.isfinite(X)] = max_e

    X_new = SelectKBest(chi2, k=num_features).fit(X, y).get_feature_names_out()
    adata_ora = sc.AnnData(X[X_new])
    adata_ora.obs["PFS_6M_Timepoint"] = y.astype("category")
    sc.pp.scale(adata_ora)
    
    adata_ora = transfer_dict(adata, adata_ora, adata_ora.var.index)
    
    return adata_ora

def analyse_single_factor(adata, factor, ontology):
    """
    str:factor = factor of interest
    str: ontology = gseapy's ontologies
    adata: adata = the output from the 
    
    """
    gene_list = adata.uns["factor_loadings"][factor]

    enr = gp.enrichr(gene_list=gene_list, # or "./tests/data/gene_list.txt",
                     gene_sets=[ontology],
                     organism='human', # don't forget to set organism to the one you desired! e.g. Yeast
                     outdir=None, # don't write to disk
                    )

    output = enr.res2d[enr.res2d["Adjusted P-value"] < 0.05]
    output["factor"] = factor
    return output

def perform_enrichment(adata_res, ontology = 'MSigDB_Hallmark_2020'):
    refdf = pd.concat([analyse_single_factor(adata_res, i, ontology) for i in adata_res.var.index])
    return refdf
