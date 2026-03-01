import numpy as np
import pandas as pd
import anndata as ad

def locate_enriched_factors(adata):
    facs = [i for i,j in adata.uns["enrichment_results"].items() if len(j) > 0]
    return facs

def retrieve_terms_for_selected_factors(x, factorlist):
    # earmarked removal move to other package
    en = x.uns["enrichment_results"]

    terms = []

    for i in factorlist:
        termie = en[i]
        terms.append(termie)

    flat_list = set([item for sublist in terms for item in sublist])

    len(set(flat_list))

    return set(flat_list)

def search_factors_for_string(x, target_string):
    present = []
    for i,j in x.uns["enrichment_results"].items():
        matching = [s for s in j if target_string in s]
        if len(matching) != 0:
            present.append(i)
    return present

def search_factors_for_genes(x, target_string):
    present = []
    for i,j in x.uns["factor_loadings"].items():
        try:
            matching = [s for s in j if target_string in s]
            if len(matching) != 0:
                present.append(i)
        except:
            pass
    return present
