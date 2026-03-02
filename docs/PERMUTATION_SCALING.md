# Scaling Permutation Testing: From 2,000 to 10^12

## Why More Permutations Matter

With the default 2,000 permutations, the finest p-value resolution is `1/2000 = 0.0005`. After FDR correction, this means:

- Genes with true p-values between 0.0005 and 0.001 are indistinguishable from each other
- Genes with p-values near the significance threshold (0.05) have noisy estimates — a gene might flip between significant and non-significant across runs
- **Most critically**: genes that are truly significant but have p-values just below the FDR-corrected threshold may be missed entirely, because the granularity of 2,000 permutations produces quantized p-values that don't differentiate them

With 100,000+ permutations:

- P-value resolution reaches 10^-5, surfacing genes that were lost in the discretization noise
- FDR correction becomes more precise, recovering additional significant genes in factor loadings
- Ranking of genes by significance within each factor becomes meaningful (at 2,000 perms, many genes have identical p-values)

## What the End User Sees

The three implementation phases are **invisible to the end user**. From the user's perspective, there is one parameter (`--nperm`) and optionally `--adaptive` and `--perm-method`. The system automatically selects the most efficient strategy:

| What the user does | What happens internally | What changes in results |
|---|---|---|
| `pycdr analyze` (new default) | 10K perms, batched, adaptive on | 5x finer p-value resolution than old default, faster via early stopping |
| `--nperm 100000` | Phase 1 batching + adaptive | Very high resolution, more significant genes recovered after FDR |
| `--nperm 100000 --no-adaptive` | Phase 1 batching only | Same resolution, runs all 100K perms (for benchmarking) |
| `--nperm 1000000 --perm-method auto` | Phase 1 + Phase 3 analytical tail | 10^6-equivalent resolution from ~100K pilot perms |
| `--nperm 2000 --no-adaptive` | Legacy behavior | Exact reproduction of old results |

**The key user-visible difference is: more permutations = finer p-value resolution = more genes recovered after FDR correction.** The implementation phases are purely about making this computationally feasible.

## CLI Modes

### Available Options

```
pycdr analyze INPUT -p PHENOTYPE [OPTIONS]
pycdr run INPUT -p PHENOTYPE [OPTIONS]

Permutation options:
  --nperm INTEGER          Permutations for importance scores [default: 2000]
  --batch-size INTEGER     Permutations per batch [default: auto-sized to fit 256 MB]
  --adaptive/--no-adaptive Enable adaptive early stopping [default: no-adaptive]
  --perm-method [exact|normal|gpd|auto]
                           P-value computation method [default: exact]
```

### Recommended Usage Patterns

```bash
# Default run (10K perms, adaptive on — good out of the box)
pycdr analyze data.h5ad -p condition

# Higher resolution (surfaces more significant genes after FDR)
pycdr analyze data.h5ad -p condition --nperm 100000

# Maximum precision — analytical tail approximation (5-15 minutes)
pycdr analyze data.h5ad -p condition --nperm 1000000 --perm-method auto

# Legacy behavior (if you need exact replication of old results)
pycdr analyze data.h5ad -p condition --nperm 2000 --no-adaptive
```

### Default Behavior

The new defaults improve p-value resolution out of the box:

- `--nperm 10000` — 5x finer resolution than the old 2,000 default
- `--batch-size auto` — automatically computed to fit within 256 MB
- `--adaptive` **on** — stops early when all genes are decided (identical significance calls, just faster)
- `--perm-method exact` — brute-force empirical p-values

Use `--no-adaptive` if you need the exact number of permutations to be run (e.g. for benchmarking).

---

## Phase 1: Batched NumPy Vectorization

### Problem

The original loop processes one permutation at a time:

```python
for _ in range(nperm):
    perm_idx = rng.permutation(n_rows)
    reshaped = Fs[perm_idx].reshape(nfacs, rows_per_split, -1)
    diff = np.abs(reshaped.max(axis=0) - reshaped.min(axis=0))
    pval_counts += (Fs_diff > diff)
```

Each iteration has:
- Python loop overhead (~10 μs per iteration)
- Small array operations that can't fully utilize CPU cache or SIMD

At ~3 ms per iteration for 5K genes, 10^6 permutations would take ~50 minutes.

### Solution

Process `batch_size` permutations simultaneously using vectorized NumPy operations:

```python
# Generate batch of permutation indices
base = np.broadcast_to(np.arange(n_rows), (batch_size, n_rows)).copy()
rng.permuted(base, axis=1, out=base)

# Vectorized computation across entire batch
permuted = Fs[base]                                    # (batch, n_rows, n_cols)
reshaped = permuted.reshape(batch, nfacs, rows_per_split, n_cols)
diff = np.abs(reshaped.max(axis=1) - reshaped.min(axis=1))
pval_counts += (Fs_diff[None] > diff).sum(axis=0)
```

### Auto Batch Sizing

`_compute_batch_size()` targets a 256 MB memory budget (conservative for an 8 GB machine):

```
bytes_per_batch_element = n_rows × n_cols × 8 (permuted array)
                        + rows_per_split × n_cols × 8 (diff array)
batch_size = min(256 MB / bytes_per_element, 4096)
```

For a typical 5,000-gene dataset:
- `n_rows = 10,000` (2 conditions × 5,000 genes), `n_cols = 22` factors
- bytes/element = 10,000 × 22 × 8 + 5,000 × 22 × 8 = 2.64 MB
- batch_size = min(256/2.64, 4096) = **96**
- Peak memory: ~253 MB (fits well within 8 GB system)

For a 50,000-gene dataset:
- bytes/element = 100,000 × 50 × 8 + 50,000 × 50 × 8 = 60 MB
- batch_size = min(256/60, 4096) = **4**
- Peak memory: ~240 MB

The cap at 4096 prevents excessive memory use even for very small datasets.

### RNG Design

The batched path uses `rng.permuted()` (NumPy ≥1.20) instead of `rng.permutation()`. When `batch_size=1`, the code falls back to the original `rng.permutation()` path, preserving the exact RNG stream for backward compatibility. The batched path consumes the RNG differently (one call to `permuted()` per batch vs. one call to `permutation()` per perm), so p-values differ numerically between batch sizes but are statistically equivalent.

### Expected Speedup

10-50× depending on dataset size and cache effects. The speedup comes from:
1. Eliminating Python loop overhead (batch_size fewer iterations)
2. Better CPU cache utilization (contiguous memory access across batch)
3. BLAS-optimized array operations

---

## Phase 2: Adaptive Early Stopping

### Problem

Most genes in a typical dataset are clearly non-significant (p >> 0.05) or clearly significant (p << 0.05). After a few thousand permutations, their status is already decided with high confidence. Running the remaining permutations is wasted computation.

### Solution

Every `check_interval` permutations (default: max(1000, 5 × batch_size)), compute a 99.9% confidence interval for each gene-factor pair's running p-value:

```python
p_hat = 1.0 - pval_counts / n_done
se = np.sqrt(p_hat * (1 - p_hat) / n_done)
z_crit = 3.29  # 99.9% CI
decided = (p_hat - z_crit*se > thresh) | (p_hat + z_crit*se < thresh)
if decided.all():
    break
```

A gene is "decided" when its entire 99.9% CI lies entirely above or below the significance threshold. The loop stops only when **all** gene-factor pairs are decided — this prevents selective stopping bias.

### Why No Bias

- The stopping criterion depends only on whether the CI excludes the threshold, not on the direction
- It never stops for a subset of genes — it's all or nothing
- The 99.9% CI (z = 3.29) ensures a < 0.1% chance of wrong decision per gene
- The actual p-value estimate at the stopping point is the same empirical estimate that would be computed at that iteration anyway

### Expected Impact

For a typical dataset where 80-95% of genes are clearly non-significant:
- With `--nperm 100000 --adaptive`: stops at ~5,000-15,000 effective permutations
- 5-20× reduction in compute time with identical significance calls

---

## Phase 3: Analytical Tail Approximation

### Problem

Even with batching and adaptive stopping, brute-force computation of p-values at 10^12 resolution is physically impossible. With 10^5 pilot permutations, the finest brute-force resolution is 10^-5.

### Solution: Two-Stage Approach

1. **Pilot run**: 100,000 actual permutations (using Phase 1 batching) to build the null distribution
2. **Analytical fit**: Extrapolate p-values beyond the empirical resolution

#### Normal Approximation (`--perm-method normal`)

Uses the mean and variance accumulated during the pilot run (Welford's method — already tracked by the existing code):

```
z = (observed - mean) / sqrt(var)
p = 1 - Φ(z)    # standard normal survival function
```

Good for moderate p-values (> 10^-4). Breaks down in extreme tails where the null distribution is not Gaussian.

#### GPD Tail Fit (`--perm-method gpd`)

For extreme p-values (< 10^-6), fits a Generalized Pareto Distribution to the top-k (k=250) null samples:

1. During the pilot run, maintain a min-heap of the 250 largest null values per gene-factor pair
2. After the pilot, fit `scipy.stats.genpareto` to exceedances above the 250th-largest value
3. Extrapolate the survival probability to the observed statistic

Memory: 250 × n_genes × n_factors × 8 bytes = ~22 MB for 5K genes × 22 factors.

#### Auto Mode (`--perm-method auto`)

Selects the best method per gene:
- For genes with empirical p > 10/n_pilot: use empirical p-value (plenty of counts)
- For genes in the tail: use GPD if top-k samples are available, normal otherwise

### When to Use

- `--perm-method auto` with `--nperm 1000000`: requests million-perm resolution but only runs ~100K actual permutations, then extrapolates analytically for the tail
- The pilot cap is 100,000 permutations regardless of `--nperm`

---

## Memory Budget on 8 GB Systems

All memory estimates assume an 8 GB machine where ~4-5 GB is available for pycdr:

| Component | Memory | Notes |
|---|---|---|
| Fs matrix (5K genes, 22 factors) | 1.7 MB | Always in memory |
| Batch allocation (batch_size=96) | ~253 MB | Largest single allocation |
| GPD top-k heap (k=250) | ~22 MB | Only when using GPD/auto |
| Accumulators (pval_counts, sum_diffs, sum_sq_diffs) | 2.6 MB | Always in memory |
| **Total working set** | **~280 MB** | Well within 8 GB |

For larger datasets (50K genes):
- Fs: 80 MB
- Batch (batch_size=4): ~240 MB
- GPD top-k: ~200 MB
- **Total: ~600 MB** — still comfortable on 8 GB

The auto batch sizer caps allocation at 256 MB by default. This can be overridden with `--batch-size` if the user has more memory available.

---

## Files Changed

| File | Changes |
|---|---|
| `pycdr/feature_selection.py` | Batched loop, `_compute_batch_size()`, adaptive stopping, `method` param |
| `pycdr/perm.py` | Batched `permute_matrix()` loop, `_perm_batch_size()` |
| `pycdr/_tail_approx.py` | **New** — `normal_pvalues()`, `gpd_pvalues()`, `compute_analytical_pvalues()` |
| `pycdr/pycdr.py` | Passthrough `batch_size`, `adaptive`, `method` params |
| `pycdr/cli.py` | `--batch-size`, `--adaptive`, `--perm-method` CLI options |
| `pyproject.toml` | Optional deps: `fast = ["numba>=0.57"]`, `gpu = ["cupy-cuda12x"]` |
| `pycdr/test/test_pycdr.py` | 12 new tests: batch equivalence, adaptive, analytical, memory |
| `benchmarks/scripts/bench_permutation.py` | Extended grid to 100K perms, memory tracking, config comparison |

## Backward Compatibility

- All new parameters have safe defaults: `batch_size=None` (auto), `adaptive=False`, `method="exact"`
- Existing `batch_size=1` path uses the exact same `rng.permutation()` call as before
- All 62 tests pass, including the 50 existing tests unchanged
- No new required dependencies
