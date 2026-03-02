
import warnings

import numpy as np
import scipy.sparse as ss
import logging
import time
from sklearn.decomposition import TruncatedSVD
from .feature_selection import get_significant_genes
from .feature_selection import calculate_minmax

logger = logging.getLogger(__name__)


def run_CDR_analysis(data, phenotype, capvar=0.95, nperm=2000, thres=0.05,
                     seed=42, correction="fdr_bh", quiet=False, **kwargs):
    """Main CDR-g analysis function

        The key step in CDR-g is an SVD-decomposition on gene
        co-expression matrices. Depending on the sequencing platform,
        this SVD step can produce thousands of factor loadings.
        By default, CDR-g selects number of factor loadings which
        captures 95% of variance in the dataset.

    Args:
        data (anndata): anndata object of interest
        phenotype (str): condition of interest
        capvar (float, optional): factor loadings to examine. Defaults to 0.95.
        nperm (int, optional): nperms to determine importance score. Defaults to 2000.
        thres (float, optional): cut-off for permutation importance to select genes. Defaults to 0.05.
        seed (int, optional): Random seed for reproducibility. Defaults to 42.
        correction (str, optional): Multiple testing correction method.
            ``"fdr_bh"`` for Benjamini-Hochberg FDR, ``"none"`` for no
            correction. Defaults to ``"fdr_bh"``.
        quiet (bool, optional): If True, suppress progress bars. Defaults to False.

    Returns:
        None. Results are stored in-place on *data.uns*: ``Fs``, ``Ls``,
        ``Fk``, ``Lk``, ``Fs_diff``, ``zscores``, ``pval_mat``,
        ``selection``, ``factor_loadings``, and ``selected_loading``.
    """
    # Deprecation shim: accept old 'pernum' kwarg
    if "pernum" in kwargs:
        warnings.warn(
            "Parameter 'pernum' is deprecated; use 'nperm' instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        nperm = kwargs.pop("pernum")
    if kwargs:
        raise TypeError(f"Unexpected keyword arguments: {list(kwargs.keys())}")

    start = time.time()

    n_groups = data.obs[phenotype].nunique()
    if n_groups < 2:
        raise ValueError(
            f"CDR-g requires at least 2 phenotype groups, but '{phenotype}' "
            f"has {n_groups}. Provide a condition column with multiple groups."
        )

    n_cells = data.X.shape[0]
    n_genes = data.X.shape[1]

    logger.info('Processing dataset of %s genes x %s cells', n_genes, n_cells)
    logger.info('Phenotype column: %s', phenotype)
    logger.info("Running SVD and threshold selection")

    t0 = time.time()
    cdr_core(data, phenotype, capvar, seed=seed)
    t1 = time.time()
    logger.info("SVD and varimax rotation: %.1fs", t1 - t0)

    logger.info("Permutation testing: %d permutations, threshold=%.3f", nperm, thres)
    npheno = data.uns["n_pheno"]

    t2 = time.time()
    get_significant_genes(data, npheno, permnum=nperm, thres=thres, seed=seed,
                          correction=correction, quiet=quiet)
    t3 = time.time()
    logger.info("Permutation testing: %.1fs", t3 - t2)

    end = time.time()
    timediff = end - start
    numfact = data.uns["selected_loading"]
    logger.info('Selected %d factor loadings', numfact)
    logger.info('Total time: %.1fs', timediff)

    
def svd_and_concatenate(matrixlist, capvar, seed=42):
    """Concatenate per-condition correlation matrices and perform truncated SVD.

    Args:
        matrixlist (list): Per-condition expression matrices (genes x cells each).
        capvar (float): Cumulative variance threshold for selecting components.
        seed (int, optional): Random seed for SVD. Defaults to 42.

    Returns:
        tuple: (Ek, Ss, X, y) — right singular vectors, singular values,
            concatenated correlation matrix, and number of selected components.
    """
    list_of_dense = [d.toarray() if ss.issparse(d) else d for d in matrixlist]
    list_of_corr_mats = [np.corrcoef(d) for d in list_of_dense]
    X = np.concatenate(list_of_corr_mats, axis=1)
    nan_count = np.isnan(X).sum()
    if nan_count > 0:
        logger.warning("Correlation matrix contains %d NaN values (replaced with 0.0)", nan_count)
    X = np.nan_to_num(X, nan=0.0)
    logger.debug("Correlation matrices computed, shape: %s", X.shape)

    _, y, Ek, Ss = get_optimal_threshold(X, capvar, seed=seed)
    return Ek, Ss, X, y


def process_svd_to_factors(Ek, Ss, N_k):
    """function for rotation and flips"""
    Ek = Ek.T
    ind = np.argsort(Ss)[::-1]
    Ss = Ss[ind]
    Ek = Ek[:, ind]
    Lk = Ss**2  # singular values to eigenvalues
    Fk = (Lk[:N_k]**0.5)*Ek[:, :N_k]  # factor loadings
    # Varimax rotation of the factor loadings
    ROT = classic_orthomax(Fk, gamma=1)  # finding rotation (gamma=1 implyes at CLASSIC varimax)
    Fs = np.dot(Fk, ROT)  # rotated factor loadings
    Ls = np.diag(ROT.T@np.diag(Lk[:N_k])@ROT)  # rotated eigenvalues
    
    ind = np.argsort(Ls)[::-1]
    Ls = Ls[ind]
    Fs = Fs[:, ind] 
    
    Fs = flip_Ek(Fs)
    logger.debug("Factor matrix shape after varimax: %s", Fs.shape)

    return Fs, Ls, Fk, Lk


def extract_matrix_from_anndata(ad, pheno_column):
    """Split the expression matrix by phenotype groups.

    Args:
        ad (anndata.AnnData): AnnData object.
        pheno_column (str): Column in ``ad.obs`` identifying conditions.

    Returns:
        tuple: (list_of_matrices, n_groups, condition_labels) — transposed
            expression matrices per condition, the number of conditions, and
            the condition label strings (in the same order as the matrices).
    """
    phenotypes = ad.obs[pheno_column].unique()
    ind = [ad.obs[pheno_column] == p for p in phenotypes]
    rands = [ad[i,:].X.T for i in ind]
    condition_labels = [str(p) for p in phenotypes]
    return rands, len(rands), condition_labels

# functions for generating pvals and integrating whole varimax


def compute_dominant_condition(Fs, n_conditions, condition_labels):
    """Identify the dominant condition for each factor from the Fs matrix.

    Reshapes Fs to ``(n_conditions, n_genes, n_factors)``, computes the mean
    absolute loading per condition per factor, and returns the condition with
    the highest mean for each factor.

    Args:
        Fs (numpy.ndarray): Factor loading matrix of shape
            ``(n_genes * n_conditions, n_factors)``.
        n_conditions (int): Number of phenotype groups.
        condition_labels (list[str]): Condition names in the same order as the
            blocks of Fs.

    Returns:
        list[str]: Dominant condition label for each factor.
    """
    n_genes = Fs.shape[0] // n_conditions
    # (n_conditions, n_genes, n_factors)
    reshaped = Fs.reshape(n_conditions, n_genes, -1)
    # mean absolute loading per condition per factor -> (n_conditions, n_factors)
    mean_abs = np.abs(reshaped).mean(axis=1)
    dominant_idx = mean_abs.argmax(axis=0)
    result = [condition_labels[i] for i in dominant_idx]
    logger.debug("Dominant conditions: %s", result)
    return result


def cdr_core(ad, pheno, capvar, seed=42):
    """Run the core CDR pipeline: SVD, varimax rotation, and store results in-place.

    Args:
        ad (anndata.AnnData): AnnData object.
        pheno (str): Phenotype column in ``ad.obs``.
        capvar (float): Cumulative variance threshold for component selection.
        seed (int, optional): Random seed for SVD. Defaults to 42.
    """
    matlist, numpheno, condition_labels = extract_matrix_from_anndata(ad, pheno)
    Ee, Ss, _, N  = svd_and_concatenate(matlist, capvar, seed=seed) # specify algorithm
    Fs, Ls, Fk, Lk = process_svd_to_factors(Ee, Ss, N)
    ad.uns["selected_loading"] = N
    ad.uns["Fs"] = Fs
    ad.uns["Ls"] = Ls
    ad.uns["Fk"] = Fk
    ad.uns["Lk"] = Lk
    ad.uns["n_pheno"] = numpheno
    ad.uns["Fs_diff"] = calculate_minmax(Fs, numpheno)
    ad.uns["condition_labels"] = condition_labels
    ad.uns["dominant_condition"] = compute_dominant_condition(Fs, numpheno, condition_labels)
        
# leos' aux functions for performing varimax


def classic_orthomax(Phi, gamma=1, q=200, tol=1e-6):
    """Returns the orthomax rotation"""
    p, k = Phi.shape
    R = np.eye(k)
    d = 0
    converged = False
    for i in range(q):
        d_old = d
        Lambda = np.dot(Phi, R)
        u, s, vh = np.linalg.svd(
            np.dot(Phi.T, np.asarray(Lambda)**3 - (gamma/p) * np.dot(Lambda, np.diag(np.diag(np.dot(Lambda.T, Lambda)))))
        )
        R = np.dot(u, vh)
        d = np.sum(s)
        if d_old != 0 and d/d_old < 1 + tol:
            converged = True
            break

    if not converged:
        logger.warning("Varimax rotation did not converge in %d iterations", q)

    return R


def flip_Ek(Ek):
    """Flip eigenvectors so that the dominant direction points upward.

    Args:
        Ek (numpy.ndarray): Factor loading matrix of shape (n_genes, n_factors).

    Returns:
        numpy.ndarray: Flipped factor loading matrix (same shape).
    """
    n, m = Ek.shape
    e_k_to_flip = abs(Ek.min(axis=0)) > Ek.max(axis=0)
    flip = np.ones(m)
    flip[e_k_to_flip] *= -1

    Ek *= flip
    return Ek


# aux functions for detecting factors.


def get_optimal_threshold(num, thres, ncomp=2000, seed=42):
    """Select the number of SVD components that capture *thres* cumulative variance.

    Args:
        num (numpy.ndarray): Concatenated correlation matrix.
        thres (float): Cumulative explained variance threshold (e.g. 0.95).
        ncomp (int, optional): Maximum components to compute. Defaults to 2000.
        seed (int, optional): Random seed for SVD. Defaults to 42.

    Returns:
        tuple: (x, y, X, v) — cumulative variance array, selected component
            count, right singular vectors for those components, and their
            singular values.
    """
    nrows = num.shape[0]
    numgenes = num.shape[1]
    ncomp = min(ncomp, nrows - 1, numgenes - 1)
    svd = TruncatedSVD(n_components=ncomp, n_iter=5, random_state=seed)
    svd.fit(num)
    x = np.cumsum(svd.explained_variance_ratio_)
    y = np.argmax(x > thres)
    if y == 0:
        y = ncomp
    logger.debug("SVD selected %d components (cumvar=%.4f)", y, x[y-1] if y > 0 else 0)
    X = svd.components_[0:y]
    v = svd.singular_values_[0:y]
    return x, y, X, v


