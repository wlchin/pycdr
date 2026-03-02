"""Analytical tail approximation for permutation p-values.

Provides two methods for estimating p-values beyond the brute-force resolution
of a pilot permutation run:

1. **Normal approximation** — uses moment-matching from Welford accumulators
   (mean and variance of the null distribution). Good for moderate p-values
   (> 10^-4).

2. **GPD (Generalized Pareto Distribution)** — fits the tail of the null
   distribution using the top-k exceedances. Accurate for extreme p-values
   (< 10^-6) where the normal approximation breaks down.

The ``"analytical"`` / ``"auto"`` method applies normal approximation for
non-extreme genes and GPD for genes in the tail.
"""

import logging

import numpy as np
from scipy import stats

logger = logging.getLogger(__name__)


def normal_pvalues(observed, mean, var):
    """Compute p-values via normal approximation (moment-matching).

    Args:
        observed: Array of observed test statistics (Fs_diff).
        mean: Per-element mean of the null distribution.
        var: Per-element variance of the null distribution.

    Returns:
        numpy.ndarray: Estimated p-values (same shape as observed).
    """
    sd = np.sqrt(np.maximum(var, 1e-30))
    z = (observed - mean) / sd
    # P(null >= observed) = 1 - Phi(z)
    pvals = stats.norm.sf(z)
    return np.clip(pvals, 0.0, 1.0)


def gpd_pvalues(observed, top_k_diffs, n_pilot):
    """Compute p-values via GPD (Generalized Pareto Distribution) tail fit.

    For each gene-factor pair, fits a GPD to the top-k exceedances of the null
    distribution and extrapolates the survival probability at the observed value.

    Args:
        observed: Array of observed test statistics, shape (n_genes, n_factors).
        top_k_diffs: Array of top-k null samples, shape (n_genes, n_factors, k).
        n_pilot: Number of pilot permutations completed.

    Returns:
        numpy.ndarray: Estimated p-values (same shape as observed).
    """
    n_genes, n_factors = observed.shape
    k = top_k_diffs.shape[2]
    pvals = np.ones((n_genes, n_factors))

    for g in range(n_genes):
        for f in range(n_factors):
            samples = top_k_diffs[g, f]
            # Remove -inf entries (unfilled slots)
            valid = samples[np.isfinite(samples)]
            if len(valid) < 10:
                # Not enough samples; fall back to empirical
                pvals[g, f] = np.nan
                continue

            threshold = np.min(valid)
            exceedances = valid - threshold

            if np.all(exceedances == 0):
                # All top-k values are identical
                pvals[g, f] = np.nan
                continue

            try:
                c, loc, scale = stats.genpareto.fit(exceedances, floc=0)
                # P(X > observed | X > threshold)
                obs_excess = observed[g, f] - threshold
                if obs_excess <= 0:
                    # Observed is below all top-k; use empirical count
                    pvals[g, f] = np.nan
                    continue
                tail_prob = stats.genpareto.sf(obs_excess, c, loc=0, scale=scale)
                # Scale by fraction of null that exceeds threshold
                pvals[g, f] = (len(valid) / n_pilot) * tail_prob
            except Exception:
                pvals[g, f] = np.nan

    return np.clip(pvals, 0.0, 1.0)


def compute_analytical_pvalues(observed, pval_counts, n_pilot, mean, var,
                               top_k_diffs=None, method="analytical"):
    """Compute p-values using analytical approximation after a pilot run.

    Args:
        observed: Observed test statistics (Fs_diff), shape (n_genes, n_factors).
        pval_counts: Count of permutations where observed > null.
        n_pilot: Number of pilot permutations completed.
        mean: Per-element mean of null distribution.
        var: Per-element variance of null distribution.
        top_k_diffs: Top-k null samples for GPD fitting, or None.
        method: ``"normal"``, ``"gpd"``, or ``"analytical"`` (auto-blend).

    Returns:
        numpy.ndarray: Estimated p-values.
    """
    empirical = 1.0 - (pval_counts / n_pilot)

    if method == "normal":
        pvals = normal_pvalues(observed, mean, var)
        logger.info("Normal approximation: computed %d p-values", pvals.size)
        return pvals

    if method == "gpd":
        if top_k_diffs is None:
            logger.warning("GPD requested but no top-k samples; falling back to normal")
            return normal_pvalues(observed, mean, var)
        gpd_p = gpd_pvalues(observed, top_k_diffs, n_pilot)
        # Fill NaN with normal approximation
        normal_p = normal_pvalues(observed, mean, var)
        mask = np.isnan(gpd_p)
        gpd_p[mask] = normal_p[mask]
        logger.info("GPD tail fit: %d/%d fell back to normal",
                    mask.sum(), gpd_p.size)
        return gpd_p

    # "analytical" = auto-blend: use empirical for clearly non-significant,
    # normal for moderate tail, GPD for extreme tail
    pvals = empirical.copy()

    # Use normal approximation where empirical resolution is limited
    # (p < 10/n_pilot suggests we're near the resolution limit)
    resolution_limit = 10.0 / n_pilot
    in_tail = empirical < resolution_limit

    if in_tail.any():
        normal_p = normal_pvalues(observed, mean, var)

        if top_k_diffs is not None:
            gpd_p = gpd_pvalues(observed, top_k_diffs, n_pilot)
            # Use GPD where available, normal as fallback
            tail_pvals = np.where(np.isnan(gpd_p), normal_p, gpd_p)
        else:
            tail_pvals = normal_p

        pvals[in_tail] = tail_pvals[in_tail]
        logger.info("Analytical blend: %d/%d used tail approximation",
                    in_tail.sum(), pvals.size)

    return np.clip(pvals, 0.0, 1.0)
