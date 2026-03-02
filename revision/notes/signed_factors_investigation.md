# Signed Factor Loadings: Investigation and Implementation Plan

## 1. The Problem

PyCDR's gene selection currently treats all significant genes within a factor as a single undifferentiated set. But factor loadings from SVD + varimax have **sign** — a gene can be positively or negatively loaded on a factor, and its loading can be higher in one condition vs another. This directional information is biologically critical:

- **Positive loading, up in STIM**: gene is activated by stimulation as part of this co-expression program
- **Negative loading, up in CTRL**: gene is repressed by stimulation as part of this program
- These represent **different biological mechanisms** (activation vs repression) that should not be lumped together

Currently, all of this information is discarded at `calculate_minmax()` via `abs(max - min)`.

## 2. Where Sign Information Exists and Where It's Lost

### Data flow through the pipeline

```
Fs matrix (SIGNED)                    ← SVD + varimax, shape (n_cond * n_genes, n_factors)
    │
    ├── flip_Ek()                     ← canonicalizes sign direction per factor (preserves relative signs)
    │
    ├── calculate_minmax()            ← abs(max - min) across conditions → Fs_diff (UNSIGNED)  ★ SIGN LOST HERE
    │
    ├── permutation test on Fs_diff   ← compares unsigned magnitudes → p-values (UNSIGNED)
    │
    ├── FDR correction on p-values    ← adjusts p-values → q-values (UNSIGNED, orthogonal to sign)
    │
    └── selection = q < thresh        ← boolean mask → gene lists (NO SIGN)
```

### The Fs matrix structure

Fs has shape `(n_conditions × n_genes, n_factors)`. After reshaping to `(n_conditions, n_genes, n_factors)`, each gene-factor pair has a **vector of loadings across conditions**:

```
Gene ISG15, Factor 3:
  CTRL loading:  0.02
  STIM loading:  0.85
  → Fs_diff = |0.85 - 0.02| = 0.83  (significant, but we don't know direction)
  → Mean loading = 0.44 (positive overall)
  → Direction: UP in STIM
```

### Empirical evidence (monocyte CD14, FDR q<0.05)

| Metric | Value |
|--------|-------|
| Total gene-factor pairs surviving FDR | 1,331 |
| Positive mean loading | 1,286 (96.6%) |
| Negative mean loading | 45 (3.4%) |
| Factors with ALL positive genes | 1,104 |
| Factors with ALL negative genes | 7 |
| Factors with MIXED signs | 6 |

The dominance of positive loadings (96.6%) is an artifact of `flip_Ek()`, which canonicalizes each factor so that its maximum absolute loading is positive. The 3.4% negative genes are genuinely anti-correlated with the factor's dominant program.

### Direction of change (which condition is higher)

| Direction | Gene-factor pairs |
|-----------|-------------------|
| Higher in CTRL | 721 (54%) |
| Higher in STIM | 610 (46%) |

This is roughly balanced — the pipeline selects genes that differ across conditions regardless of which direction. But this information is not exposed to the user.

## 3. Does FDR Correction Impact Directionality?

**No.** The FDR correction is orthogonal to the sign problem. It operates on p-values that were already computed from unsigned Fs_diff values:

```
Fs (signed) → Fs_diff (unsigned) → p-values (unsigned) → FDR correction (unsigned) → selection
                                                           ↑
                                                    FDR sits here — no sign to lose
```

FDR neither helps nor hurts directionality. It correctly adjusts for multiple testing on the **two-sided** question "does this gene's loading differ across conditions?" — which is the right question to ask before splitting by direction.

## 4. Proposed Implementation: `--direction` CLI Parameter

### Design

Add a `direction` parameter with three modes:

| Mode | Behavior | Gene set keys | Default? |
|------|----------|---------------|----------|
| `"unsigned"` | Current behavior. One gene list per factor. | `factor.0`, `factor.1`, ... | **Yes** (backward compatible) |
| `"split"` | After selection, split each factor into up/down gene sets based on direction of loading change across conditions. | `factor.0.up`, `factor.0.down`, ... | |

### How it works

1. **Gene selection is unchanged** — the two-sided permutation test and FDR correction remain the same. The `abs(max - min)` statistic is correct for a two-sided test of "does this gene differ across conditions?"

2. **Post-selection annotation** — after determining which genes are significant, look up the sign information from the original Fs matrix:
   ```python
   reshaped = Fs.reshape(n_conditions, n_genes, n_factors)
   # For each significant gene in factor fi:
   per_cond = reshaped[:, gene_idx, fi]
   dominant_cond_idx = np.argmax(np.abs(per_cond))  # or use pre-computed dominant_condition
   if per_cond[dominant_cond_idx] > per_cond.mean():
       direction = "up"
   else:
       direction = "down"
   ```

3. **FDR is calculated identically** — the statistical test doesn't change. We're just annotating the results post-hoc.

### Why NOT two one-sided tests?

An alternative would be to run separate one-sided permutation tests ("is this gene significantly UP in condition A?" and "is this gene significantly DOWN?"). I don't recommend this because:

- **Halves statistical power** — each one-sided test uses only half the tail
- **Doubles multiple testing burden** — 2× as many tests to correct
- **Unnecessary** — the sign information is already available in Fs without needing a statistical test. The permutation test establishes that the loading **differs** more than chance; the Fs matrix tells us **which direction**.

### Default preference

**`unsigned` should be the default** because:
- Backward compatible with all existing results
- More statistical power (two-sided test)
- Simpler gene sets for initial exploration
- Users who need direction can opt in with `--direction split`

## 5. Files to Modify

| File | Change |
|------|--------|
| `pycdr/feature_selection.py` | Add `direction` parameter to `get_significant_genes()`. When `"split"`, compute per-gene direction from Fs and store as `factor.X.up` / `factor.X.down`. Store direction matrix in `adata.uns["gene_direction"]`. |
| `pycdr/pycdr.py` | Thread `direction` through `run_CDR_analysis()` |
| `pycdr/cli.py` | Add `--direction` option (`type=click.Choice(["unsigned", "split"])`, default `"unsigned"`) |
| `pycdr/utils.py` | `get_top_genes()`: add `direction` column when `gene_direction` matrix exists. `output_results()`: handle `.up`/`.down` factor keys. |
| `pycdr/reporting.py` | Display direction mode in params. Update guide text. Factor summary table: show up/down counts when direction is split. |
| `pycdr/test/test_pycdr.py` | New tests for split mode, direction correctness. |

## 6. Detailed Implementation

### `feature_selection.py` changes

```python
def get_significant_genes(adata, nfacs, permnum=2000, thres=0.05, seed=42,
                          correction="fdr_bh", direction="unsigned", quiet=False):
    selection, pmat, zcore = select_modules(adata, permnum, thres, nfacs, seed=seed,
                                            correction=correction, quiet=quiet)

    Fs = adata.uns["Fs"]
    n_genes = Fs.shape[0] // nfacs
    reshaped = Fs.reshape(nfacs, n_genes, -1)  # (n_conditions, n_genes, n_factors)

    # Compute per-gene direction: +1 if loading increases toward dominant condition, -1 otherwise
    # Use mean loading across conditions as sign indicator
    mean_loading = reshaped.mean(axis=0)  # (n_genes, n_factors)
    # Direction of change: which condition has higher loading?
    # cond_diff[g, f] > 0 means loading is higher in last condition
    cond_diff = reshaped[-1] - reshaped[0]  # (n_genes, n_factors) for 2 conditions

    factor_loadings = {}
    num_loadings = selection.shape[1]

    for i in range(num_loadings):
        sigs_mask = selection[:, i]
        sigs = adata.var.index[sigs_mask].to_list()

        if direction == "split" and len(sigs) > 0:
            # Split into up (positive cond_diff) and down (negative cond_diff)
            up_genes = [g for g, m in zip(sigs, cond_diff[sigs_mask, i]) if m > 0]
            down_genes = [g for g, m in zip(sigs, cond_diff[sigs_mask, i]) if m <= 0]
            factor_loadings[f"factor.{i}.up"] = up_genes
            factor_loadings[f"factor.{i}.down"] = down_genes
        else:
            factor_loadings[f"factor.{i}"] = sigs

    adata.uns["factor_loadings"] = factor_loadings

    if direction == "split":
        adata.uns["gene_direction"] = np.sign(cond_diff)  # +1/-1 matrix
```

### What "up" and "down" mean

For a dataset with conditions `[CTRL, STIM]`:
- `factor.3.up` = genes whose loading is **higher in STIM** than CTRL (activated by stimulation)
- `factor.3.down` = genes whose loading is **higher in CTRL** than STIM (repressed by stimulation)

For >2 conditions, the direction is relative to the condition with the highest absolute mean loading (the dominant condition).

### Impact on enrichment

When `direction="split"`, the enrichment modules (`kruskal.py`, `perm.py`) would receive smaller, directionally coherent gene sets. This could **improve enrichment power** because:
- ssGSEA scores are more interpretable when all genes in a set move in the same direction
- A gene set mixing up- and down-regulated genes may show weak enrichment even when both directions are individually strong
- Pathway databases (GO, KEGG) often distinguish between activation and repression

### Impact on reporting

The HTML report and terminal summary would show:
- Factor count includes both `.up` and `.down` sub-factors
- Gene counts shown as "X up / Y down" instead of a single number
- Top genes table annotated with direction
- The guide text explains what up/down means in the context of the phenotype conditions

## 7. Proposed Experiments

### Experiment 1: Direction coherence validation

**Question:** Are the "up" and "down" gene sets within a factor biologically coherent?

**Method:**
1. Run PyCDR on monocyte data with `direction="split"`
2. For the top factors (by gene count), run GO enrichment separately on up vs down gene sets
3. Compare: do up-genes enrich for activation pathways and down-genes for repression?

**Expected:** In the IFN-stimulated factor, up-genes should be ISGs (interferon-stimulated genes), while down-genes should be housekeeping/baseline genes repressed during IFN response.

### Experiment 2: Enrichment power comparison

**Question:** Does splitting by direction improve enrichment significance?

**Method:**
1. Run Kruskal-Wallis enrichment on unsigned factor gene sets
2. Run Kruskal-Wallis enrichment on direction-split gene sets
3. Compare FDR values: are split sets more often significant?

**Expected:** Split sets should have equal or better enrichment FDR because the ssGSEA scores are computed on directionally coherent gene sets.

### Experiment 3: Comparison with DE gene lists

**Question:** Do PyCDR's directional gene sets overlap with classical differential expression (CTRL vs STIM)?

**Method:**
1. Run Wilcoxon rank-sum DE between CTRL and STIM
2. Classify DE genes as up-in-STIM or down-in-STIM
3. Compute overlap with PyCDR's `factor.X.up` and `factor.X.down` gene sets
4. Measure: are up-genes enriched in DE-up, and down-genes in DE-down?

**Expected:** Strong overlap, confirming that PyCDR's directional annotation recovers known biology.

### Experiment 4: Impact on method comparison (cNMF/PyWGCNA/Hotspot)

**Question:** Does direction-aware PyCDR produce more comparable results to other methods?

**Method:**
1. Run PyCDR (unsigned and split), cNMF, PyWGCNA, Hotspot on monocyte data
2. Compute gene set overlap (Jaccard index) between methods
3. Compare: is the overlap better with directional PyCDR?

**Expected:** cNMF and Hotspot also produce directional gene programs, so PyCDR's split mode should improve cross-method concordance.

## 8. Varimax Convergence and NaN Values

These are pre-existing issues unrelated to FDR/directionality, but worth noting:

### NaN values in correlation matrix

**Cause:** Genes with zero variance in one condition (e.g., not expressed in any cell of that condition) produce `NaN` in `np.corrcoef()`. On the test dataset: 29 zero-variance genes in condition 0, 27 in condition 24 → 9,630 NaN values.

**Current handling:** `np.nan_to_num(X, nan=0.0)` — replaces NaN with 0 (no correlation). This is reasonable but could be improved by filtering zero-variance genes before correlation computation.

**Impact on directionality:** None directly, but genes with zero variance in one condition will have extreme loading differences (loading only in one condition), which could bias direction-split gene sets.

### Varimax non-convergence

**Cause:** `classic_orthomax()` uses `q=20` max iterations. With 1,672 factors on the monocyte data, 20 iterations may be insufficient for convergence.

**Impact:** Non-converged varimax means factor loadings may not be fully rotated to simple structure. This affects ALL downstream analysis (gene selection, direction, enrichment) but is not specifically related to FDR correction.

**Fix (separate PR):** Increase `q` to 100 or 200, or add a convergence check that raises a warning with suggested action.
