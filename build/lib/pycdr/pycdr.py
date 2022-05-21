#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 13:12:56 2021
@author: weeloongchin
"""

import numpy as np
import scipy.sparse as ss
import logging
import time
import warnings
from .feature_selection import get_significant_genes
from .feature_selection import calculate_minmax

warnings.simplefilter("ignore") 
logging.basicConfig(format='%(process)d - %(levelname)s : %(asctime)s - %(message)s', level=logging.DEBUG)

logger = logging.getLogger(__name__)


def run_CDR_analysis(data, phenotype, capvar = 0.95, pernum = 2000, thres = 0.05):
    """
    input:
    :param data: Anndata object
    phenotype: column in obs df
    backend: dask vs sparse
    num_factor_loadings: choose factor loadings
    num_perm: choose permutation number
    pval_threshold: pval post correction
    gene_quantile: threshold for declaring interesting genes
    
    output:
    Anndata object with results appended to uns 
    
    """
    start = time.time()
    
    gene_num = data.X.shape[0]
    cell_num = data.X.shape[1]
    
    logger.info('processing dataset of %s genes X %s cells', cell_num, gene_num)
    logger.info('target class label:: %s', phenotype)
    logger.info("SVD and threshold selection")
    res = pvalgenerator(data, phenotype, capvar)

    logger.info("completed SVD and varimax")
    logger.info("permutation testing for gene sets:: perms:: %s threshold :: %s", pernum, thres)
    npheno= data.uns["n_pheno"]
    #get_significant_genes_perms(data, npheno, permnum = pernum, thres = thres)
    get_significant_genes(data, npheno, permnum = pernum, thres = thres)
    
    logger.info("computed thresholds for gene selection")

    end = time.time()
    timediff = end - start
    numfact = data.uns["selected_loading"]
    logger.info('N factor loadings:: %s', numfact)
    logger.info('wall clock time in seconds:: %s', timediff)

    
def dask_ver(matrixlist, capvar):
    """provides svd and concatenation with dask"""
    import dask.array as da
    from dask_ml.decomposition import TruncatedSVD

    if ss.issparse(matrixlist[0]):
        list_of_mats_as_dask_arrays = [da.from_array(np.array(d.todense())) for d in matrixlist]
    else:
        list_of_mats_as_dask_arrays = [da.from_array(d) for d in matrixlist]

    list_of_corr_mats = [da.corrcoef(d) for d in list_of_mats_as_dask_arrays]
    X = da.concatenate(list_of_corr_mats, axis=1)
    X[da.isnan(X)] = 0.0

    _, y, Ek, Ss = get_optimal_threshold(X, capvar)
  
    #Ek = svd.components_
    #Ss = svd.singular_values_
    return Ek, Ss, X, y 

def process_svd_to_factors(Ek, Ss, N_k):
    """function for rotation and flips"""
    Ek = Ek.T
    ind = np.argsort(Ss)[::-1]
    Ss  = Ss[ind]
    Ek  = Ek[:, ind]
    
    Lk = Ss**2 # singular values to eigenvalues
    Fk = (Lk[:N_k]**0.5)*Ek[:,:N_k] # factor loadings
    
    # Varimax rotation of the factor loadings
    ROT = classic_orthomax(Fk, gamma=1) # finding rotation (gamma=1 implyes at CLASSIC varimax)
    Fs  = np.dot(Fk,ROT) # rotated factor loadings
    Ls = np.diag(ROT.T@np.diag(Lk[:N_k])@ROT)  # rotated eigenvalues
    
    ind = np.argsort(Ls)[::-1]
    Ls = Ls[ind]
    Fs  = Fs[:, ind] 
    
    Fs = flip_Ek(Fs)
    
    return Fs, Ls, Fk, Lk


### aux functions for matrix extraction

def get_numbers_of_pheno(ad, pheno):
    """return list of nums"""
    vals = ad.obs[pheno].value_counts().tolist()
    return vals

def get_bools_of_pheno(ad, pheno):
    """return list of booleans"""
    phenotypes = ad.obs[pheno].unique()
    bool_list = [ad.obs[pheno] == i for i in phenotypes]
    return bool_list

def extract_matrix_from_anndata(ad, pheno_column):
    ind = get_bools_of_pheno(ad, pheno_column)
    rands = [ad[i,:].X.T for i in ind]
    return rands, len(rands)

#### functions for generating pvals and integrating whole varimax

def _full_Fs(ad, pheno, capvar):
    matlist, numpheno = extract_matrix_from_anndata(ad, pheno)
    Ee, Ss, _, N  = dask_ver(matlist, capvar) # specify algorithm
    Fs, Ls, Fk, Lk = process_svd_to_factors(Ee, Ss, N)
    ad.uns["selected_loading"] = N
    ad.uns["Fs"] = Fs
    ad.uns["Ls"] = Ls
    ad.uns["Fk"] = Fk
    ad.uns["Lk"] = Lk
    ad.uns["n_pheno"] = numpheno
    Fs_diff = calculate_minmax(Fs, numpheno)
    return Fs_diff

def pvalgenerator(ad, pheno, capvar):
    Fs_diff = _full_Fs(ad, pheno, capvar)
    ad.uns["Fs_diff"] = Fs_diff

    return Fs_diff
    
        
# leos' aux functions 

def classic_orthomax(Phi, gamma = 1, q = 20, tol = 1e-6):
    """Returns the orthomax rotation"""
    from numpy import eye, asarray, dot, sum, diag
    from numpy.linalg import svd
    p,k = Phi.shape
    R = eye(k)
    d=0
    for i in range(q):
        d_old = d
        Lambda = dot(Phi, R)
        u,s,vh = svd(dot(Phi.T,asarray(Lambda)**3 - (gamma/p) * dot(Lambda, diag(diag(dot(Lambda.T,Lambda))))))
        R = dot(u,vh)
        d = sum(s)
        if d_old!=0 and d/d_old < 1 + tol: break

    return R

def flip_Ek(Ek):
    """That functions guaranties that the eigenvectors will "point up".
    """
    n, m = Ek.shape
    
    e_k_to_flip = abs(Ek.min(axis=0)) > Ek.max(axis=0)
    
    flip = np.ones(m)
    flip[e_k_to_flip] *= -1

    Ek *= flip
    
    return Ek

### aux functions for detecting factors.

def get_optimal_threshold(num, thres, ncomp = 2000):
    """
    selects number of factors for truncated SVD
    
    """
    from dask_ml.decomposition import TruncatedSVD
    import dask.array as da
    nrows = num.shape[0] # this shows num cells and is required for svd

    numgenes = num.shape[1] # this is to make sure if less 2000
    if numgenes < ncomp:
        ncomp = numgenes - 1
    print(ncomp)
    numm = num.rechunk((nrows, 10))
    svd = TruncatedSVD(n_components=ncomp, n_iter=5, random_state=42)
    svd.fit(numm)
    x = np.cumsum(svd.explained_variance_ratio_)
    y = np.argmax(x>thres)
    if y == 0:
        y = ncomp
    X = svd.components_[0:y]
    v = svd.singular_values_[0:y]    
    return x, y, X, v


