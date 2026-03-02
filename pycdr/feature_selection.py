import numpy as np
import logging
import tqdm
from statsmodels.stats.multitest import fdrcorrection

logger = logging.getLogger(__name__)

# Default memory budget for batch sizing (256 MB)
_DEFAULT_MEM_BUDGET = 256 * 1024 * 1024
_MAX_BATCH_SIZE = 4096


def _compute_batch_size(n_rows, n_cols, nfacs, mem_budget=_DEFAULT_MEM_BUDGET):
    """Compute batch size that fits within a memory budget.

    Each batch element requires ``n_rows * n_cols * 8`` bytes (float64)
    for the permuted matrix, plus the reshaped and diff arrays.

    Args:
        n_rows (int): Total rows in Fs.
        n_cols (int): Number of columns (factors).
        nfacs (int): Number of phenotype splits.
        mem_budget (int): Memory budget in bytes. Defaults to 256 MB.

    Returns:
        int: Batch size clamped to [1, _MAX_BATCH_SIZE].
    """
    # Main allocation: permuted (batch, n_rows, n_cols) float64
    # Plus reshaped (batch, nfacs, rows_per_split, n_cols) -- same memory
    # Plus diff (batch, rows_per_split, n_cols) -- smaller
    bytes_per_elem = n_rows * n_cols * 8  # permuted array
    rows_per_split = n_rows // nfacs
    bytes_per_elem += rows_per_split * n_cols * 8  # diff array
    if bytes_per_elem <= 0:
        return 1
    batch = max(1, mem_budget // bytes_per_elem)
    return min(batch, _MAX_BATCH_SIZE)


def calculate_zscore(base, av, var):
    zscore = (base - av) / np.sqrt(np.maximum(var, 1e-10))
    return zscore


def select_modules(adata, nperm, thresh, nfacs, seed=42, correction="fdr_bh",
                   quiet=False, batch_size=None, adaptive=False,
                   method="exact"):
    """Permutation testing for gene significance across factor loadings.

    Args:
        adata (anndata.AnnData): AnnData object with ``Fs`` in ``.uns``.
        nperm (int): Number of permutations.
        thresh (float): P-value significance threshold.
        nfacs (int): Number of phenotype groups.
        seed (int): Random seed.
        correction (str): ``"fdr_bh"`` or ``"none"``.
        quiet (bool): Suppress progress bars.
        batch_size (int or None): Permutations per batch. None = auto-size.
        adaptive (bool): Enable adaptive early stopping.
        method (str): ``"exact"``, ``"normal"``, ``"gpd"``, or ``"auto"``.
    """
    Fs = adata.uns["Fs"]
    Fs_diff = calculate_minmax(Fs, nfacs)

    n_rows = Fs.shape[0]
    rows_per_split = n_rows // nfacs
    rng = np.random.default_rng(seed)
    n_cols = Fs.shape[1]

    logger.info("Running %d permutations across %d genes x %d factors",
                nperm, rows_per_split, n_cols)

    # Decide effective method
    if method == "auto":
        effective_method = "exact" if nperm <= 100_000 else "analytical"
    else:
        effective_method = method

    # For analytical methods, cap the pilot run
    pilot_nperm = nperm
    if effective_method in ("normal", "gpd", "analytical"):
        pilot_nperm = min(nperm, 100_000)
        logger.info("Analytical mode: running %d pilot permutations", pilot_nperm)

    # GPD: track top-k null samples per gene-factor pair
    gpd_k = 250
    need_gpd = effective_method in ("gpd", "analytical")
    if need_gpd:
        # Min-heap via sorted insertion; store top-k largest diffs
        top_k_diffs = np.full((rows_per_split, n_cols, gpd_k), -np.inf)

    # Auto batch sizing
    if batch_size is None:
        batch_size = _compute_batch_size(n_rows, n_cols, nfacs)
    batch_size = max(1, min(batch_size, pilot_nperm))

    logger.info("Batch size: %d", batch_size)

    pval_counts = np.zeros((rows_per_split, n_cols))
    sum_diffs = np.zeros((rows_per_split, n_cols))
    sum_sq_diffs = np.zeros((rows_per_split, n_cols))

    n_done = 0
    n_batches = (pilot_nperm + batch_size - 1) // batch_size

    # Adaptive early stopping parameters
    check_interval = max(1000, batch_size * 5)  # check every ~5 batches or 1000 perms
    early_stopped = False

    pbar = tqdm.tqdm(total=pilot_nperm, desc="Permutation testing", disable=quiet)

    for batch_idx in range(n_batches):
        current_batch = min(batch_size, pilot_nperm - n_done)

        if current_batch == 1:
            # Single-perm path (preserves exact RNG stream for batch_size=1)
            perm_idx = rng.permutation(n_rows)
            reshaped = Fs[perm_idx].reshape(nfacs, rows_per_split, -1)
            diff = np.abs(reshaped.max(axis=0) - reshaped.min(axis=0))
            pval_counts += (Fs_diff > diff)
            sum_diffs += diff
            sum_sq_diffs += diff ** 2
            if need_gpd:
                _update_topk(top_k_diffs, diff)
        else:
            # Batched path: generate batch of permutation indices
            base = np.broadcast_to(
                np.arange(n_rows), (current_batch, n_rows)
            ).copy()
            rng.permuted(base, axis=1, out=base)

            # (batch, n_rows, n_cols)
            permuted = Fs[base]
            # (batch, nfacs, rows_per_split, n_cols)
            reshaped = permuted.reshape(current_batch, nfacs, rows_per_split, n_cols)
            # (batch, rows_per_split, n_cols)
            diff = np.abs(reshaped.max(axis=1) - reshaped.min(axis=1))

            pval_counts += (Fs_diff[None] > diff).sum(axis=0)
            sum_diffs += diff.sum(axis=0)
            sum_sq_diffs += (diff ** 2).sum(axis=0)

            if need_gpd:
                for i in range(current_batch):
                    _update_topk(top_k_diffs, diff[i])

        n_done += current_batch
        pbar.update(current_batch)

        # Adaptive early stopping check
        if adaptive and n_done >= 1000 and n_done % check_interval < current_batch:
            p_hat = 1.0 - pval_counts / n_done
            se = np.sqrt(p_hat * (1 - p_hat) / n_done)
            z_crit = 3.29  # 99.9% CI
            decided = ((p_hat - z_crit * se > thresh) |
                       (p_hat + z_crit * se < thresh))
            if decided.all():
                logger.info("Adaptive early stop at %d/%d permutations "
                            "(all gene-factor pairs decided)", n_done, pilot_nperm)
                early_stopped = True
                pbar.close()
                break

    if not early_stopped:
        pbar.close()

    # Compute statistics from the pilot run
    effective_nperm = n_done
    av = sum_diffs / effective_nperm
    var = (sum_sq_diffs - sum_diffs ** 2 / effective_nperm) / max(effective_nperm - 1, 1)

    if effective_method == "exact" or nperm == effective_nperm:
        pval_mat = 1 - (pval_counts / effective_nperm)
    else:
        # Analytical tail approximation
        from ._tail_approx import compute_analytical_pvalues
        pval_mat = compute_analytical_pvalues(
            Fs_diff, pval_counts, effective_nperm, av, var,
            top_k_diffs if need_gpd else None,
            method=effective_method,
        )

    z_score_mat = calculate_zscore(Fs_diff, av, var)

    # Store raw p-values
    adata.uns["pval_mat_raw"] = pval_mat

    # Apply multiple testing correction
    if correction == "fdr_bh":
        flat_pvals = pval_mat.ravel()
        _, flat_corrected = fdrcorrection(flat_pvals)
        corrected_pval_mat = flat_corrected.reshape(pval_mat.shape)
        adata.uns["pval_mat"] = corrected_pval_mat
        selection = corrected_pval_mat < thresh
        logger.info("Applied FDR (Benjamini-Hochberg) correction")
    else:
        adata.uns["pval_mat"] = pval_mat
        selection = pval_mat < thresh

    logger.debug("Permutation complete: %d/%d genes significant",
                 selection.sum(), selection.size)

    adata.uns["zscores"] = z_score_mat
    adata.uns["Fs_diff"] = Fs_diff
    adata.uns["selection"] = selection
    adata.uns["effective_nperm"] = effective_nperm

    return selection, pval_mat, z_score_mat


def _update_topk(top_k, diff):
    """Update the top-k heap with new diff values.

    Args:
        top_k: Array of shape (rows_per_split, n_cols, k) with current top-k.
        diff: Array of shape (rows_per_split, n_cols) with new values.
    """
    # Replace minimums in top_k where diff is larger
    min_idx = top_k.argmin(axis=2)
    min_vals = np.take_along_axis(top_k, min_idx[:, :, None], axis=2).squeeze(axis=2)
    mask = diff > min_vals
    if mask.any():
        rows, cols = np.where(mask)
        for r, c in zip(rows, cols):
            k_idx = min_idx[r, c]
            top_k[r, c, k_idx] = diff[r, c]


def get_significant_genes(adata, nfacs, permnum=10000, thres=0.05, seed=42,
                          correction="fdr_bh", quiet=False, batch_size=None,
                          adaptive=True, method="exact"):
    """Identify significant genes for each factor loading via permutation testing.

    Calls :func:`select_modules` to obtain a boolean selection matrix, then
    stores per-factor gene lists in ``adata.uns["factor_loadings"]``.

    Args:
        adata (anndata.AnnData): AnnData object with ``Fs`` in ``.uns``.
        nfacs (int): Number of phenotype groups (used for min-max splitting).
        permnum (int, optional): Number of permutations. Defaults to 2000.
        thres (float, optional): P-value threshold. Defaults to 0.05.
        seed (int, optional): Random seed for reproducibility. Defaults to 42.
        correction (str, optional): Multiple testing correction method.
            ``"fdr_bh"`` for Benjamini-Hochberg FDR, ``"none"`` for no
            correction. Defaults to ``"fdr_bh"``.
        quiet (bool, optional): If True, suppress progress bars. Defaults to False.
        batch_size (int or None, optional): Permutations per batch. None = auto.
        adaptive (bool, optional): Enable adaptive early stopping. Defaults to False.
        method (str, optional): P-value method: ``"exact"``, ``"normal"``,
            ``"gpd"``, or ``"auto"``. Defaults to ``"exact"``.
    """
    selection, pmat, zcore = select_modules(
        adata, permnum, thres, nfacs, seed=seed,
        correction=correction, quiet=quiet, batch_size=batch_size,
        adaptive=adaptive, method=method,
    )

    factor_loadings = {}

    num_loadings = selection.shape[1]

    for i in range(num_loadings):
        sigs = adata.var.index[selection[:,i]].to_list()
        nameoffl = "factor." + str(i)
        factor_loadings[nameoffl] = sigs
        logger.debug("Factor '%s': %d genes", nameoffl, len(sigs))

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
    if Fs.shape[0] % splits != 0:
        raise ValueError(
            f"Fs row count ({Fs.shape[0]}) is not evenly divisible by splits ({splits})"
        )
    rows_per_split = Fs.shape[0] // splits
    reshaped = Fs.reshape(splits, rows_per_split, -1)
    return np.abs(reshaped.max(axis=0) - reshaped.min(axis=0))
