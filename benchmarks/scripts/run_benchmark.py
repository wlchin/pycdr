"""Run timed CDR analysis with phase breakdown."""

import json
import platform
import resource
import time
import tracemalloc
import sys
import os

import anndata as ad

# Add project root to path so pycdr is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from pycdr.pycdr import cdr_core
from pycdr.feature_selection import get_significant_genes


def get_peak_rss_mb():
    """Get current peak RSS in MB (handles macOS vs Linux)."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if platform.system() == "Darwin":
        return peak / 1e6  # macOS: bytes -> MB
    return peak / 1e3  # Linux: KB -> MB


def run_benchmark(h5ad_path, n_genes, n_cells, rep, capvar, nperm, thres, output_path):
    adata = ad.read_h5ad(h5ad_path)
    phenotype_col = "condition"

    rss_baseline = get_peak_rss_mb()

    # Phase 1: SVD
    tracemalloc.start()
    t0 = time.perf_counter()
    cdr_core(adata, phenotype_col, capvar)
    t1 = time.perf_counter()
    svd_time = t1 - t0
    _, svd_tracemalloc_peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    rss_after_svd = get_peak_rss_mb()
    svd_rss_delta_mb = rss_after_svd - rss_baseline

    n_factors = adata.uns["selected_loading"]
    npheno = adata.uns["n_pheno"]

    # Phase 2: Permutation testing
    tracemalloc.start()
    t2 = time.perf_counter()
    get_significant_genes(adata, npheno, permnum=nperm, thres=thres)
    t3 = time.perf_counter()
    perm_time = t3 - t2
    _, perm_tracemalloc_peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    rss_after_perm = get_peak_rss_mb()
    perm_rss_delta_mb = rss_after_perm - rss_after_svd

    result = {
        "n_genes": int(n_genes),
        "n_cells": int(n_cells),
        "nperm": int(nperm),
        "rep": int(rep),
        "svd_time_s": round(svd_time, 4),
        "perm_time_s": round(perm_time, 4),
        "total_time_s": round(svd_time + perm_time, 4),
        "n_factors": int(n_factors),
        "peak_rss_mb": round(rss_after_perm, 2),
        "svd_rss_delta_mb": round(svd_rss_delta_mb, 2),
        "perm_rss_delta_mb": round(max(perm_rss_delta_mb, 0), 2),
        "svd_tracemalloc_mb": round(svd_tracemalloc_peak / 1e6, 2),
        "perm_tracemalloc_mb": round(perm_tracemalloc_peak / 1e6, 2),
    }

    with open(output_path, "w") as f:
        json.dump(result, f, indent=2)


if __name__ == "__main__":
    run_benchmark(
        h5ad_path=str(snakemake.input[0]),
        n_genes=snakemake.params.n_genes,
        n_cells=snakemake.params.n_cells,
        rep=snakemake.wildcards.rep,
        capvar=float(snakemake.params.capvar),
        nperm=int(snakemake.params.nperm),
        thres=float(snakemake.params.thres),
        output_path=str(snakemake.output[0]),
    )
