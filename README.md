[![Tests](https://github.com/wlchin/pycdr/actions/workflows/test.yml/badge.svg)](https://github.com/wlchin/pycdr/actions/workflows/test.yml) [![Coverage](https://img.shields.io/badge/coverage-100%25-brightgreen)](https://github.com/wlchin/pycdr/actions/workflows/test.yml)


# CDR-g (CDR-genomics)

CDR-g identifies **condition-specific gene expression programs** from multi-condition single-cell RNA-seq data. It uses SVD decomposition on gene co-expression matrices followed by varimax rotation and permutation testing to recover interpretable factors, each representing a set of co-regulated genes that distinguish one condition from others.

CDR-g is described in the manuscript **"Spectral detection of condition-specific biological pathways in single-cell gene expression data"**. It extends the CDR framework introduced in [Portes & Small (2020)](https://doi.org/10.1103/PhysRevE.102.062301).

**Contents:** [How it works](#how-it-works) | [Installation](#installation) | [Preparing your data](#preparing-your-data) | [Python API](#python-api) | [CLI](#command-line-interface) | [Interpreting results](#interpreting-results) | [Example workflows](#example-workflows) | [Citation](#citation)

## How it works

1. **Co-expression matrices** are built per condition from the count matrix (`adata.X`)
2. **Truncated SVD** decomposes the concatenated correlation matrices, retaining factors that capture a target variance threshold (default 95%)
3. **Varimax rotation** produces interpretable factor loadings
4. **Permutation testing** identifies genes with statistically significant loadings on each factor
5. **Enrichment testing** (optional) assesses whether each gene program is differentially active across conditions via chi-square proportions test or Kruskal-Wallis test

## Installation

CDR-g requires Python >= 3.9. Install in a virtual environment:

```bash
pip install cdr-py

# With plotting support (for pycdr plot / pycdr report figures)
pip install cdr-py[plot]
```

## Preparing your data

CDR-g takes an [AnnData](https://anndata.readthedocs.io/) `.h5ad` file. The data must be preprocessed before running CDR-g -- normalization and log-transformation are **not** handled by pycdr.

### What CDR-g expects

- **`adata.X`** -- **log-transformed** count matrix (cells x genes). Sparse or dense.
- **`adata.obs`** -- must contain a column identifying the condition/phenotype of interest (e.g. `"stim"`, `"Hours"`)
- **`adata.var`** -- gene metadata. If running enrichment with the `perm` method, a gene name column (e.g. `"gene_short_name"`) must be present.

### Preprocessing with scanpy (recommended)

Use [scanpy](https://scanpy.readthedocs.io/) or an equivalent tool to prepare your `.h5ad`:

```python
import scanpy as sc

adata = sc.read_10x_h5("raw_data.h5")  # or any loader

# Standard preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Ensure the condition column exists in adata.obs
# e.g. adata.obs["condition"] should contain values like "ctrl", "stim"

adata.write("prepared.h5ad")
```

### Gene filtering with pycdr

After preprocessing, you can optionally filter lowly-expressed genes to reduce computation time. This is the one preprocessing step that pycdr handles directly, either via the Python API or the CLI:

```bash
# Filter by minimum number of expressing cells
pycdr filter prepared.h5ad -m numcells --count-threshold 1 --min-cells 10 -o filtered.h5ad

# Or filter by cell fraction and median count
pycdr filter prepared.h5ad -m percent --cell-fraction 0.05 --median-count 1.0 -o filtered.h5ad
```

Alternatively, `pycdr run` can filter in-line with `--filter-method`:

```bash
pycdr run prepared.h5ad -p condition --filter-method numcells --min-cells 10
```

## Python API

The CLI covers common workflows, but the Python API gives you direct access to the AnnData object at every step. This is useful when you need to integrate CDR-g into a larger analysis (e.g. combining with scanpy clustering or differential expression), inspect intermediate results programmatically, or build custom visualizations from the factor loadings and enrichment scores.

```python
import anndata as ad
from pycdr import run_CDR_analysis, calculate_enrichment
from pycdr.utils import output_results, filter_genecounts_numcells

adata = ad.read_h5ad("data.h5ad")

# Optional: filter lowly expressed genes
adata = filter_genecounts_numcells(adata, count_threshold=1, min_expressed_cells=10)

# Run CDR-g analysis
run_CDR_analysis(adata, "stim")

# Inspect results
print(adata.uns["factor_loadings"])  # dict of factor -> gene lists
print(adata.uns["Fs"].shape)         # rotated factor loading matrix
print(adata.uns["pval_mat"])         # permutation p-values per gene x factor

# Optional: enrichment testing
from pycdr.perm import calculate_enrichment
calculate_enrichment(adata, "stim", list(adata.uns["factor_loadings"].keys()),
                     nperm=100, genecol="gene_short_name", thresh=0.05)

# Export a summary table
df = output_results(adata)
df.to_csv("results.csv", index=False)
```

### Key outputs stored in `adata.uns`

| Key | Description |
|-----|-------------|
| `factor_loadings` | `dict` mapping factor names to lists of significant genes |
| `Fs` | Rotated factor loading matrix (cells x factors) |
| `Fs_diff` | Min-max difference of loadings across conditions |
| `pval_mat` | Permutation p-values (genes x factors) |
| `zscores` | Z-scores for gene significance (genes x factors) |
| `selected_loading` | Number of factors selected by variance threshold |
| `condition_labels` | List of condition group names (same order as Fs blocks) |
| `dominant_condition` | Dominant condition per factor (list of strings) |
| `pval_dict` | Enrichment test statistics and p-values (after enrichment) |

## Command-line interface

After installation, the `pycdr` CLI lets you run the full pipeline from the terminal -- useful for scripting, workflow managers (Snakemake, Nextflow), or quick exploration.

```
pycdr --help
```

### Commands

| Command | Description |
|---------|-------------|
| `pycdr run` | Full pipeline (filter + analyze + enrich + export) |
| `pycdr analyze` | Run CDR-g SVD/varimax analysis only |
| `pycdr filter` | Filter genes from an .h5ad file |
| `pycdr enrich` | Run enrichment on previously analyzed data |
| `pycdr results` | Export results table (CSV, TSV, ASCII table, or markdown) |
| `pycdr info` | Display dataset metadata and analysis summary |
| `pycdr plot` | Generate a summary figure (requires `[plot]` extra) |
| `pycdr report` | Generate a self-contained HTML report |

### Quick start

```bash
# Minimal -- run CDR-g analysis
pycdr run data.h5ad -p stim

# End-to-end: filter, analyze, enrich, and generate HTML report
pycdr run data.h5ad -p stim -o results.h5ad -c results.csv \
  --filter-method numcells --min-cells 25 \
  --enrich --enrich-method kruskal --report report.html

# Subset to a specific cell population before analysis
pycdr run data.h5ad -p stim --subset "celltype=T_cells"

# Multiple subset values (OR) and multiple columns (AND)
pycdr run data.h5ad -p stim \
  --subset "celltype=T_cells,NK_cells" --subset "batch=batch1"

# Inspect a dataset
pycdr info data.h5ad -p stim

# Export top genes for a specific factor
pycdr results analyzed.h5ad --top-genes 20 --factor 3

# View results as a readable table in the terminal
pycdr results analyzed.h5ad -f table

# Re-run enrichment with different parameters
pycdr enrich analyzed.h5ad -p stim -m kruskal -o enriched.h5ad

# Generate a summary figure
pycdr plot analyzed.h5ad -o summary.png

# Generate an HTML report
pycdr report analyzed.h5ad -p stim -o report.html
```

### Common options for `pycdr run`

| Option | Default | Description |
|--------|---------|-------------|
| `-p, --phenotype` | required | Condition column in `adata.obs` |
| `-s, --subset` | - | Subset cells: `COLUMN=VALUE[,VALUE2]`. Repeatable. |
| `-o, --output` | `{stem}_cdr.h5ad` | Output .h5ad path |
| `-c, --csv` | - | Export results CSV |
| `--filter-method` | `none` | `none`, `percent`, or `numcells` |
| `--enrich / --no-enrich` | off | Run enrichment after analysis |
| `--enrich-method` | `perm` | `perm` (chi-square) or `kruskal` (KW test) |
| `--genecol` | - | Gene name column (required for perm) |
| `--report` | - | Generate an HTML report at this path |
| `-v / -vv` | - | INFO / DEBUG logging |

Run `pycdr run --help` for the full list of options.

## Interpreting results

CDR-g identifies **factors** -- each factor is a co-expression program, a group of genes whose expression changes together between your experimental conditions. After running the pipeline, the key question is: *which factors are biologically meaningful and what do they represent?*

### Factors and gene assignments

- **Factor numbering** (factor.0, factor.1, ...) reflects eigenvalue order. Lower-numbered factors capture more variance in the data, but this does not necessarily mean they are the most biologically interesting.
- **Gene count** (`n_genes`): The number of genes assigned to a factor by permutation testing (p < 0.05). Factors with many genes represent broad transcriptional programs; those with few genes are more specific. Factors with **0 genes** had no statistically significant condition-dependent co-expression and can usually be ignored.
- **Mean z-score** (`mean_zscore`): Measures how strongly a factor's genes change across conditions compared to a permutation null. Higher values indicate stronger condition-specific co-expression. A z-score above 2 means a gene's loading difference is well above what would be expected by chance.
- **Top genes**: The genes with the highest z-scores in each factor. These are the most condition-responsive genes in the program and are the best starting point for biological interpretation -- look them up in pathway databases or the literature.
- **Dominant condition** (`dominant_condition`): The phenotype group that drives each factor most strongly. This is always computed from the Fs matrix (the condition block with the highest mean absolute loading). When enrichment is run, it is overridden by a more biologically interpretable measure: the condition with the highest median ssGSEA score (Kruskal) or the highest activation proportion (perm). This label tells you *which* condition a factor is associated with, not just *that* the factor differs between conditions.

### Enrichment testing

If you run enrichment (`--enrich`), each factor's gene set is tested for differential activity across conditions:

- **Kruskal-Wallis** (`--enrich-method kruskal`): a non-parametric test on ssGSEA enrichment scores. Reports an H-statistic -- larger values indicate stronger differences between conditions.
- **Permutation** (`--enrich-method perm`): a chi-square proportions test on binarized gene set activation. Requires `--genecol`.

Both methods report a **p-value**, a Benjamini-Hochberg **FDR**, and a **dominant_condition** column identifying which phenotype group drives each factor. Factors with **FDR < 0.05** are considered significantly differentially active between conditions and are the most important to investigate further.

### Practical tips

- Start with factors that have **many genes, high mean z-scores, and significant FDR**. These represent the dominant condition-specific programs. The **dominant condition** column tells you which phenotype group each factor is associated with.
- To interpret a factor biologically, look up its top genes in pathway databases (e.g. MSigDB, Enrichr, Gene Ontology).
- Factors that share genes may represent overlapping or related biological processes.
- Factors with very few genes (< 5) may reflect noise unless enrichment testing confirms significance -- treat them with caution.
- Use `pycdr report` to generate an HTML report that includes this guidance alongside your results.

## Example workflows

Three Snakemake workflows reproducing the manuscript analyses are provided in [`benchmarks/CDR_workflows`](https://github.com/wlchin/CDR_workflows) (included as a git submodule). They cover monocytes, myocytes, and T-cell datasets.

To run the full workflows:

```bash
pip install scanpy[leiden] bbknn enrichment_utils
```

See the [workflow README](benchmarks/CDR_workflows/README.md) for details.

## Documentation

Full API documentation is available at [cdr-g.readthedocs.io](https://cdr-g.readthedocs.io/).

## Citation

If you use CDR-g, please cite:

> Chin WL et al. "Spectral detection of condition-specific biological pathways in single-cell gene expression data." *(manuscript in preparation)*

CDR-g extends the CDR framework from:

> Portes LL, Small M. "Navigating differential structures in complex networks." *Phys Rev E.* 2020;102(6):062301. [doi:10.1103/PhysRevE.102.062301](https://doi.org/10.1103/PhysRevE.102.062301)
