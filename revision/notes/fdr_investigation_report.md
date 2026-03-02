# FDR Correction Investigation Report

## Problem

PyCDR's permutation-based gene selection (`feature_selection.py::select_modules()`) applies **no multiple testing correction** to its p-value matrix. On real data (monocyte CD14 dataset), this means:

- **19.6M hypothesis tests** (9,815 genes × 2,000 factors)
- At raw p < 0.05, **100% of genes** are flagged as significant in at least one factor
- Average gene appears in **73.9 factors** — a 73.9× redundancy ratio
- Factor gene sets are meaninglessly large and overlapping

This makes PyCDR's factor gene sets unsuitable for comparison with cNMF/PyWGCNA/Hotspot, and undermines any downstream benchmark.

## Evidence

### P-value distribution

The raw p-value distribution from the permutation test shows ~3.7% of tests at p < 0.05, close to the uniform null expectation of 5%. This is consistent with a valid permutation test producing well-calibrated p-values — the problem is not with the test itself, but with the absence of multiple testing correction.

### Root cause

In `feature_selection.py`, line 42 (prior to this fix):

```python
selection = pval_mat < thresh
```

Raw p-values are thresholded directly with no correction. Each gene is tested across all factors, and each factor tests all genes, but no adjustment is made for the multiplicity.

### Inconsistency with other modules

Both enrichment modules in PyCDR already apply FDR correction:

- `kruskal.py:106`: `_, fdr_values = fdrcorrection(valid)` — corrects across factors
- `perm.py:160`: `whole["fdr"] = fdrcorrection(ps)[1]` — corrects across factors

Only `feature_selection.py` was missing correction, despite running far more tests (n_genes per factor vs. 1 per factor in enrichment).

### Two-stage testing overview

PyCDR has two separate stages of statistical testing:

1. **Stage 1 — Gene selection** (`feature_selection.py`): "Which genes are significant for each factor?" One permutation test per gene × factor. This is where the false positive problem lives.

2. **Stage 2 — Enrichment** (`kruskal.py` / `perm.py`): "Is this factor's gene set differentially enriched across conditions?" One test per factor. Already FDR-corrected.

The fix only touches Stage 1.

## Recommendation: Per-factor Benjamini-Hochberg correction

### Why per-factor (not global)

- Each factor is an independent co-expression program (orthogonal by varimax rotation)
- The biological question per factor is "which genes matter for THIS factor?" — a family of `n_genes` tests
- Global correction across all 19.6M tests is biologically meaningless and yields near-zero significant genes
- BH is valid under independence and PRDS (positive regression dependence on a subset), which holds within factors since genes within a factor are positively correlated by construction
- Consistent with how enrichment tests in `kruskal.py` and `perm.py` correct per-factor

### Implementation

```python
from statsmodels.stats.multitest import fdrcorrection

# Per-factor BH correction
fdr_mat = np.empty_like(pval_mat)
for col in range(pval_mat.shape[1]):
    _, fdr_mat[:, col] = fdrcorrection(pval_mat[:, col], alpha=thresh)
selection = fdr_mat < thresh
```

### Storage

- `adata.uns["pval_mat"]` — always raw p-values (unchanged)
- `adata.uns["fdr_mat"]` — FDR-corrected q-values (new, only when correction applied)

### Expected impact

Dramatic reduction in false positives:
- From 100% gene coverage to a biologically meaningful subset
- Factor gene sets become specific and non-overlapping
- Enables meaningful comparison with cNMF/PyWGCNA/Hotspot
- Backward compatibility preserved via `correction="none"` option
