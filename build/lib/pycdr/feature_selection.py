import numpy as np
import anndata as ad
import tqdm
import logging

logger = logging.getLogger(__name__)
logging.disable(logging.DEBUG)

def create_initial_aggregate(adata, nfacs):
    # return the 3-tuple
    count = 0
    mean = adata.uns["Fs"]
    nph = nfacs
    diff = calculate_minmax(mean, nph)
    M2 = np.zeros(diff.shape)
    
    return (count, diff, M2)

def update_welford(existingAggregate, newValue):
    (count, mean, M2) = existingAggregate
    count += 1
    delta = newValue - mean
    mean += delta / count
    delta2 = newValue - mean
    M2 += delta * delta2
    return (count, mean, M2)

# Retrieve the mean, variance and sample variance from an aggregate

def finalize(existingAggregate):
    (count, mean, M2) = existingAggregate
    (mean, variance, sampleVariance) = (mean, M2 / count, M2 / (count - 1))
    return (mean, variance, sampleVariance)


##### now pval mat stuff #####    

def create_initial_pvalmat(adata,nfacs):
    # return the 3-tuple
    count = 0
    mean = adata.uns["Fs"]
    nph = nfacs
    diff = calculate_minmax(mean, nph)
    pval_init = np.zeros(diff.shape)
    return (count, diff, pval_init)


def update_pval(existingAggregate, newValue):
    (count, mean, M2) = existingAggregate
    count += 1
    outcome = np.greater(mean, newValue).astype(np.int64)
    M2 = M2 + outcome
    return (count, mean, M2)

## metric calculation

def calculate_zscore(base, av, var):
    zscore = (base - av)/np.sqrt(var)
    return zscore


def calculate_pvals(agg, nperm, thresh):
    pval = (1 - (agg[2]/nperm)) 
    selection = pval < thresh
    return pval, selection

def select_modules(adata, nperm, thresh, nfacs):

    stup = create_initial_aggregate(adata, nfacs)
    pval = create_initial_pvalmat(adata, nfacs)
    
    Fs_diff = pval[1]

    rng = np.random.default_rng(42) # this needs to exist outisde the loop
    mat_to_permute = adata.uns["Fs"]
    
    for i in tqdm.tqdm(range(nperm)):
        if i % 200 == 0:
            logger.info("num permutations:: %s", i)
        range_extract = rng.permutation(mat_to_permute)
        permuted_Fs_diff = calculate_minmax(range_extract, nfacs)

        stup = update_welford(stup, permuted_Fs_diff)
        pval = update_pval(pval, permuted_Fs_diff)

    
    av, var, _ = finalize(stup)
    
    z_score_mat = calculate_zscore(Fs_diff, av, var)
    
    pval_mat, selection = calculate_pvals(pval, nperm, thresh)
    
    adata.uns["zscores"] = z_score_mat
    adata.uns["pval_mat"] = pval_mat
    adata.uns["Fs_diff"] = Fs_diff
    adata.uns["selection"] = selection
    
    return selection, pval_mat, z_score_mat

def get_significant_genes(adata, nfacs, permnum = 2000, thres = 0.05):
    
    selection, pmat, zcore = select_modules(adata, permnum, thres, nfacs)
    
    factor_loadings = {}

    num_loadings = selection.shape[1]

    for i in range(num_loadings):
        sigs = adata.var.index[selection[:,i]].to_list()
        #zscore = adata.uns["zscores"][:,i].tolist()
        #pvals = adata.uns["pval_mat"][:,i].tolist()
        nameoffl = "factor." + str(i)
        factor_loadings[nameoffl] = sigs
    
    adata.uns["factor_loadings"] = factor_loadings
    
def calculate_minmax(Fs, splits):
    """the metric"""
    mats_divided = np.vsplit(Fs,splits)
    mats_divided_l = np.dstack(mats_divided)
    Fs_diff = np.abs(np.amax(mats_divided_l, axis = 2) - np.amin(mats_divided_l, axis = 2))

    return Fs_diff
