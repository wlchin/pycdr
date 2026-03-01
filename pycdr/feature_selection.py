import numpy as np
import logging

logger = logging.getLogger(__name__)


def calculate_zscore(base, av, var):
    zscore = (base - av) / np.sqrt(np.maximum(var, 1e-10))
    return zscore


def select_modules(adata, nperm, thresh, nfacs):

    Fs = adata.uns["Fs"]
    Fs_diff = calculate_minmax(Fs, nfacs)

    n_rows = Fs.shape[0]
    rows_per_split = n_rows // nfacs
    rng = np.random.default_rng(42)

    # Pre-generate all permutation indices
    indices = np.array([rng.permutation(n_rows) for _ in range(nperm)])
    # Apply all permutations: (nperm, n_rows, n_cols)
    all_perms = Fs[indices]
    # Reshape and compute minmax in batch: (nperm, nfacs, rows_per_split, n_cols)
    reshaped = all_perms.reshape(nperm, nfacs, rows_per_split, -1)
    all_diffs = np.abs(reshaped.max(axis=1) - reshaped.min(axis=1))

    # Compute stats from the batch
    av = all_diffs.mean(axis=0)
    var = all_diffs.var(axis=0, ddof=1)
    pval_counts = (Fs_diff > all_diffs).sum(axis=0)
    pval_mat = 1 - (pval_counts / nperm)

    z_score_mat = calculate_zscore(Fs_diff, av, var)
    selection = pval_mat < thresh

    adata.uns["zscores"] = z_score_mat
    adata.uns["pval_mat"] = pval_mat
    adata.uns["Fs_diff"] = Fs_diff
    adata.uns["selection"] = selection

    return selection, pval_mat, z_score_mat

def get_significant_genes(adata, nfacs, permnum = 2000, thres = 0.05):
    """Identify significant genes for each factor loading via permutation testing.

    Calls :func:`select_modules` to obtain a boolean selection matrix, then
    stores per-factor gene lists in ``adata.uns["factor_loadings"]``.

    Args:
        adata (anndata.AnnData): AnnData object with ``Fs`` in ``.uns``.
        nfacs (int): Number of phenotype groups (used for min-max splitting).
        permnum (int, optional): Number of permutations. Defaults to 2000.
        thres (float, optional): P-value threshold. Defaults to 0.05.
    """
    selection, pmat, zcore = select_modules(adata, permnum, thres, nfacs)

    factor_loadings = {}

    num_loadings = selection.shape[1]

    for i in range(num_loadings):
        sigs = adata.var.index[selection[:,i]].to_list()
        nameoffl = "factor." + str(i)
        factor_loadings[nameoffl] = sigs

    adata.uns["factor_loadings"] = factor_loadings

def calculate_minmax(Fs, splits):
    """Compute the absolute max-minus-min metric across phenotype splits.

    Reshapes the factor loading matrix into *splits* blocks and computes
    ``abs(max - min)`` along the split axis for each gene and factor.

    Args:
        Fs (numpy.ndarray): Factor loading matrix of shape (n_cells, n_factors).
        splits (int): Number of phenotype groups to split rows into.

    Returns:
        numpy.ndarray: Difference matrix of shape (rows_per_split, n_factors).
    """
    rows_per_split = Fs.shape[0] // splits
    reshaped = Fs.reshape(splits, rows_per_split, -1)
    return np.abs(reshaped.max(axis=0) - reshaped.min(axis=0))
