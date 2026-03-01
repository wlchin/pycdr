# pycdr Benchmark Suite

Scaling benchmarks for the pycdr CDR analysis pipeline across varying gene counts and cell counts, using a Snakemake workflow.

## Quick Start

```bash
cd benchmarks/
uv pip install snakemake matplotlib    # or: pip install snakemake matplotlib
snakemake -n                           # dry run — verify the DAG
snakemake --cores 1                    # full run (sequential, memory-safe)
```

Outputs:
- `figures/benchmark_scaling.png` — 3-panel scaling figure
- `results/benchmark_summary.tsv` — full timing and memory data

## Benchmark Grid

| Parameter      | Values                                    |
|----------------|-------------------------------------------|
| `n_genes`      | 500, 1000, 2000, 3000, 5000, 7000, 10000 |
| `n_cells`      | 200, 500, 1000, 2000, 5000, 10000        |
| `n_perms`      | 500                                       |
| `n_phenotypes` | 2 (balanced split)                        |
| `n_repeats`    | 3 (median used for reporting)             |

42 grid points x 3 repeats = 126 benchmark runs.

## Workflow DAG

```
generate_data  -->  run_benchmark  -->  collect_results  -->  plot_results
 (42 jobs)          (126 jobs)           (1 job)              (1 job)
```

- **generate_data**: Creates synthetic h5ad files (one per grid point, no repeat wildcard).
- **run_benchmark**: Runs timed CDR analysis with phase breakdown. Uses Snakemake's `benchmark:` directive for wall-clock and max RSS tracking.
- **collect_results**: Aggregates per-run JSONs and Snakemake benchmark TSVs into `benchmark_summary.tsv`.
- **plot_results**: Produces a 3-panel matplotlib figure (heatmap, line plot, stacked bars).

## Synthetic Data

Data is generated from a latent factor model that mimics real scRNA-seq co-expression structure:

1. **30 latent gene modules** — cell-level factors and gene-level loadings create realistic gene-gene correlations
2. **Condition-specific perturbations** — latent factors are shifted per phenotype group
3. **Poisson counts + log1p transform** — matches CDR's expected input format
4. **Gene-specific base rates** — drawn from Uniform(1, 20) in log-space

This structure ensures the SVD captures 95% variance in ~25-28 factors (matching typical HVG-filtered scRNA-seq data), keeping the batched permutation array memory-bounded.

## Results Summary

### Runtime (median total seconds)

| n_genes |  200 cells |  500 cells | 1000 cells | 2000 cells | 5000 cells | 10000 cells |
|--------:|-----------:|-----------:|-----------:|-----------:|-----------:|------------:|
|     500 |        0.2 |        0.2 |        0.2 |        0.2 |        0.3 |         0.3 |
|    1000 |        0.8 |        0.8 |        0.8 |        0.8 |        0.8 |         0.9 |
|    2000 |        4.3 |        5.4 |        4.3 |        4.5 |        4.5 |         4.5 |
|    3000 |        7.7 |        8.1 |        7.8 |        8.1 |        8.1 |        11.5 |
|    5000 |       15.1 |       15.2 |       15.5 |       15.4 |       16.2 |        16.2 |
|    7000 |       25.7 |       25.9 |       25.2 |       25.7 |       26.1 |        26.5 |
|   10000 |       46.3 |       46.9 |       46.9 |       47.8 |       48.7 |        50.7 |

### Peak Memory (MB)

| n_genes |  200 cells | 10000 cells | Increase |
|--------:|-----------:|------------:|---------:|
|     500 |        197 |         213 |      +8% |
|    1000 |        410 |         442 |      +8% |
|    3000 |      1,230 |       1,326 |      +8% |
|    5000 |      2,049 |       2,210 |      +8% |
|    7000 |      2,869 |       3,122 |      +9% |
|   10000 |      5,808 |       6,200 |      +7% |

### Key Findings

1. **Gene count dominates both time and memory** — runtime scales roughly as O(n_genes^2) due to gene-gene correlation matrix computation and SVD. Memory scales similarly (correlation matrices are n_genes x n_genes).

2. **Cell count has minimal impact** — going from 200 to 10,000 cells (50x increase) adds only ~7-8% memory and <10% runtime. This is because:
   - `np.corrcoef` output is always (n_genes x n_genes) regardless of cell count
   - TruncatedSVD operates on the correlation matrix, not the raw data
   - The permutation array depends on n_genes and n_factors, not n_cells

3. **SVD phase dominates runtime** — at all gene counts, the SVD/varimax phase accounts for >95% of wall-clock time. The batched permutation testing phase is comparatively cheap.

4. **8 GB RAM is sufficient up to ~10,000 genes** — the largest grid point (10,000 genes x 10,000 cells) peaks at 6.2 GB. Use `--cores 1` to avoid parallel memory stacking.

## File Structure

```
benchmarks/
  Snakefile                     # Workflow definition
  config.yaml                   # Grid parameters
  README.md                     # This file
  scripts/
    generate_synthetic.py       # Synthetic h5ad generator
    run_benchmark.py            # Timed CDR with phase breakdown
    collect_results.py          # Aggregate JSONs into TSV
    plot_benchmark.py           # 3-panel matplotlib figure
  CDR_workflows/                # Existing submodule (untouched)
```

Runtime-generated (gitignored):
```
  data/                         # .h5ad files
  results/                      # Per-run JSONs + benchmark_summary.tsv
  benchmarks_sm/                # Snakemake benchmark TSVs
  figures/                      # Final plots
  .snakemake/                   # Snakemake internal state
```

## Memory Considerations

The dominant memory consumers in the CDR pipeline are:

1. **Gene-gene correlation matrices** (n_genes x n_genes x 8 bytes, per condition) — quadratic in gene count, independent of cell count
2. **Concatenated correlation matrix** (n_genes x 2*n_genes x 8 bytes) — fed to TruncatedSVD
3. **Batched permutation array** (n_perms x n_genes x n_factors x 8 bytes) — depends on how many SVD factors are selected

For the synthetic data with 30 latent modules, the SVD selects ~25-28 factors at capvar=0.95, keeping the permutation array small (~0.3-1.1 GB across the grid).

To estimate peak memory for a given gene count:

```
peak_MB ≈ 3 * n_genes^2 * 8 / 1e6  (correlation matrices)
         + n_perms * n_genes * n_factors * 8 / 1e6  (permutation array)
```

For 10,000 genes: ~2,400 MB (correlation) + 1,080 MB (permutations) + overhead ≈ 6 GB.
