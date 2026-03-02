"""Standalone permutation-step benchmark.

Generates synthetic data inline, runs SVD once per gene count, then sweeps
nperm values to measure how permutation cost scales.

Includes comparison of batch sizes and adaptive early stopping.

Usage:
    cd benchmarks && python scripts/bench_permutation.py
"""

import json
import os
import sys
import time
import tracemalloc

import anndata as ad
import numpy as np
import pandas as pd

# Add project root to path so pycdr is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from pycdr.pycdr import cdr_core
from pycdr.feature_selection import select_modules


# --- Benchmark grid ---
N_GENES_LIST = [500, 2000, 5000, 10000]
N_PERM_LIST = [100, 500, 2000, 5000, 10000, 50000, 100000]
N_CELLS = 1000
N_PHENOTYPES = 2
CAPVAR = 0.95
SEED = 42

# Memory limit for batch sizing (8 GB machine constraint)
MEM_BUDGET = 2 * 1024 * 1024 * 1024  # 2 GB working set (conservative for 8 GB system)


def generate_synthetic_adata(n_genes, n_cells, n_phenotypes, seed):
    """Generate a synthetic AnnData object inline (no disk I/O)."""
    rng = np.random.default_rng(seed)

    gene_log_means = np.log(rng.uniform(1, 20, size=n_genes))

    n_latent = 30
    cell_factors = rng.normal(0, 1, size=(n_cells, n_latent))
    gene_loadings = rng.normal(0, 1, size=(n_latent, n_genes))

    labels = np.repeat(
        [f"condition_{i}" for i in range(n_phenotypes)],
        n_cells // n_phenotypes,
    )
    if len(labels) < n_cells:
        labels = np.append(labels, ["condition_0"] * (n_cells - len(labels)))

    for c in range(n_phenotypes):
        mask = labels == f"condition_{c}"
        shift = rng.normal(0, 0.3, size=n_latent)
        cell_factors[mask] += shift

    log_rates = gene_log_means + (cell_factors @ gene_loadings) * 0.3
    rates = np.exp(np.clip(log_rates, -2, 6))
    counts = rng.poisson(lam=rates).astype(np.float32)
    X = np.log1p(counts)

    obs = pd.DataFrame({"condition": labels}, index=[f"cell_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])

    return ad.AnnData(X=X, obs=obs, var=var)


def measure_memory(func, *args, **kwargs):
    """Run func while tracking peak memory usage. Returns (result, peak_mb)."""
    tracemalloc.start()
    result = func(*args, **kwargs)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return result, peak / (1024 * 1024)


def run_bench():
    results = []

    for n_genes in N_GENES_LIST:
        print(f"\n--- n_genes={n_genes}, n_cells={N_CELLS} ---")

        # Generate data and run SVD once
        adata = generate_synthetic_adata(n_genes, N_CELLS, N_PHENOTYPES, SEED)

        t0 = time.perf_counter()
        cdr_core(adata, "condition", CAPVAR, seed=SEED)
        svd_time = time.perf_counter() - t0

        n_factors = adata.uns["selected_loading"]
        npheno = adata.uns["n_pheno"]
        print(f"  SVD: {svd_time:.3f}s  (n_factors={n_factors})")

        # Sweep nperm values
        for nperm in N_PERM_LIST:
            # Skip very large nperm for small gene counts (would be too slow)
            if nperm > 10000 and n_genes < 2000:
                continue

            configs = [
                {"label": "batch_auto", "batch_size": None, "adaptive": False, "method": "exact"},
                {"label": "batch_1", "batch_size": 1, "adaptive": False, "method": "exact"},
            ]
            # Only test adaptive for larger nperm
            if nperm >= 5000:
                configs.append(
                    {"label": "adaptive", "batch_size": None, "adaptive": True, "method": "exact"}
                )
            # Test analytical for very large nperm
            if nperm >= 50000:
                configs.append(
                    {"label": "auto_method", "batch_size": None, "adaptive": False, "method": "auto"}
                )

            for cfg in configs:
                adata_copy = adata.copy()

                tracemalloc.start()
                t1 = time.perf_counter()
                select_modules(
                    adata_copy, nperm, 0.05, npheno, seed=SEED, quiet=True,
                    batch_size=cfg["batch_size"],
                    adaptive=cfg["adaptive"],
                    method=cfg["method"],
                )
                perm_time = time.perf_counter() - t1
                _, peak_bytes = tracemalloc.get_traced_memory()
                tracemalloc.stop()
                peak_mb = peak_bytes / (1024 * 1024)

                effective_nperm = adata_copy.uns.get("effective_nperm", nperm)
                total_time = svd_time + perm_time
                perm_fraction = perm_time / total_time if total_time > 0 else 0
                per_perm_ms = (perm_time / max(effective_nperm, 1)) * 1000

                rec = {
                    "n_genes": int(n_genes),
                    "n_cells": int(N_CELLS),
                    "nperm": int(nperm),
                    "effective_nperm": int(effective_nperm),
                    "config": cfg["label"],
                    "svd_time": round(svd_time, 4),
                    "perm_time": round(perm_time, 4),
                    "total_time": round(total_time, 4),
                    "perm_fraction": round(perm_fraction, 4),
                    "per_perm_ms": round(per_perm_ms, 4),
                    "peak_mem_mb": round(peak_mb, 1),
                    "n_factors": int(n_factors),
                }
                results.append(rec)
                print(f"  nperm={nperm:7d}  {cfg['label']:12s}  perm={perm_time:.3f}s  "
                      f"eff={effective_nperm:7d}  "
                      f"frac={perm_fraction:.1%}  per_perm={per_perm_ms:.3f}ms  "
                      f"mem={peak_mb:.1f}MB")

    return results


def print_table(results):
    """Print a formatted summary table to stdout."""
    df = pd.DataFrame(results)
    print("\n" + "=" * 120)
    print("PERMUTATION STEP BENCHMARK SUMMARY")
    print("=" * 120)
    print(f"{'n_genes':>8} {'nperm':>8} {'config':>12} {'svd_s':>8} {'perm_s':>8} "
          f"{'eff_nperm':>9} {'total_s':>8} "
          f"{'perm%':>7} {'ms/perm':>8} {'mem_MB':>8} {'n_fac':>6}")
    print("-" * 120)
    for _, row in df.iterrows():
        print(f"{int(row['n_genes']):8d} {int(row['nperm']):8d} {row['config']:>12s} "
              f"{row['svd_time']:8.3f} "
              f"{row['perm_time']:8.3f} {int(row['effective_nperm']):9d} "
              f"{row['total_time']:8.3f} "
              f"{row['perm_fraction']:6.1%} {row['per_perm_ms']:8.3f} "
              f"{row['peak_mem_mb']:8.1f} "
              f"{int(row['n_factors']):6d}")
    print("=" * 120)


def main():
    results = run_bench()
    print_table(results)

    # Save JSON
    out_dir = os.path.join(os.path.dirname(__file__), "..", "results")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "bench_permutation.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nJSON saved to {out_path}")


if __name__ == "__main__":
    main()
