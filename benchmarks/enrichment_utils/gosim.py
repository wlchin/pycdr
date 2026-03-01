# module for calculating semantic similarity
# aim to condense variation in anndata object

# download gene ontology and go slim

from goatools import obo_parser
from goatools.semantic import semantic_similarity
import pandas as pd
from sklearn.cluster import KMeans
import numpy as np

def load_go_dag(path):
    """load go DAG from file

    Args:
        path (str): filepath for obo file

    Returns:
        goDAG: GOATOOLS-parsed DAG
    """
    go = obo_parser.GODag(path)
    return go

def extract_go_terms(adata, factor):
    """given anndata, extract terms from factor

    Args:
        adata (anndata): analysed anndata object
        factor (str): factor of interest

    Returns:
        list: list of go terms
    """
    factor_dict = adata.uns["enrichment_details"][factor]
    go_terms1 = [i["goid"] for i in factor_dict.values()]
    return go_terms1

def get_similarity_matrix(go_terms1, go_terms2, go):
    """given two lists of go terms, return similarity

    Args:
        go_terms1 (list): go terms from factor loading
        go_terms2 (list): go terms from factor loading
        go (_type_): goDAG

    Returns:
        np.array: similarity matrix of terms
    """
    term_list = []

    for i in go_terms1:
        for j in go_terms2:
            try:
                sim = semantic_similarity(i, j, go)
            except:
                sim = 0
            term_list.append(sim)

    mat_x_axis = len(go_terms1)

    sim_matrix = np.array(term_list).reshape(mat_x_axis, -1)
    return sim_matrix

def apply_similarity_metric(sim_matrix):
    """condenses term similarity to single float value

    Args:
        sim_matrix (np.array): similarity matrix of terms

    Returns:
        float: similarity value
    """
    sim = np.mean(sim_matrix)
    return sim

def get_non_empty_factors(adata):
    """inspects anndata and returns factor list

    Args:
        adata (anndata): processed anndata
    Returns:
        list: factors with enrichment
    """
    factor_list = [i for i,j in adata.uns["enrichment_results"].items() if len(j) > 0]
    return factor_list

def get_similarity_adata(mono, go):
    """similarity matrix on factors

    Args:
        mono (anndata): processed anndata
        go (GOATOOOLS dag): goDAG object

    Returns:
        np.array: similarity matrix on factors
    """

    sims = []

    non_empty_factors = get_non_empty_factors(mono)

    for i in non_empty_factors:
        for j in non_empty_factors:
            try:
                first_go_list =  extract_go_terms(mono, i)
                second_go_list = extract_go_terms(mono, j)
                sim_mat = get_similarity_matrix(first_go_list, second_go_list, go)
                comb_sim_value = apply_similarity_metric(sim_mat)
            except:
                comb_sim_value = 0
            sims.append(comb_sim_value)

    dim_similarity_matrix = len(non_empty_factors)
    arr = np.array(sims).reshape(dim_similarity_matrix, -1)
    return arr

def group_redundant_factors(arr_sim, nclusters, mono):
    """kmeans on sim matrix

    Args:
        arr_sim (np.array): similarity matrix on factors
        nclusters (int): target groups
        mono (anndata): processed anndata

    Returns:
        pd.DataFrame: pandas dataframe
    """

    kmeans = KMeans(n_clusters=nclusters, random_state=0).fit(arr_sim)
    factor_df = pd.DataFrame([get_non_empty_factors(mono), kmeans.labels_]).T
    return factor_df
