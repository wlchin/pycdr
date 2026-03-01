[![Tests](https://github.com/wlchin/pycdr/actions/workflows/test.yml/badge.svg)](https://github.com/wlchin/pycdr/actions/workflows/test.yml) [![Coverage](https://img.shields.io/badge/coverage-100%25-brightgreen)](https://github.com/wlchin/pycdr/actions/workflows/test.yml)


# CDR-g (CDR-genomics)

CDR-g identifies **condition-specific gene expression programs** from multi-condition single-cell RNA-seq data. It uses SVD decomposition on gene co-expression matrices followed by varimax rotation and permutation testing to recover interpretable factors — each representing a set of co-regulated genes that distinguish one condition from others.

CDR-g is described in the manuscript **"Spectral detection of condition-specific biological pathways in single-cell gene expression data"**. It extends the CDR framework introduced in [Portes & Small (2020)](https://doi.org/10.1103/PhysRevE.102.062301).

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
```

## Input data requirements

CDR-g takes an [AnnData](https://anndata.readthedocs.io/) `.h5ad` file with:

- **`adata.X`** — log-transformed count matrix (cells x genes). Sparse or dense.
- **`adata.obs`** — must contain a column for the condition/phenotype of interest (e.g. `"stim"`, `"Hours"`)
- **`adata.var`** — gene metadata. If running enrichment with the `perm` method, a gene name column (e.g. `"gene_short_name"`) must be present.

Gene filtering (by expression frequency or count thresholds) is recommended to reduce computation time on large datasets.

## Python API

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
| `pval_dict` | Enrichment test statistics and p-values (after enrichment) |

## Command-line interface

After installation, the `pycdr` CLI lets you run the full pipeline from the terminal — useful for scripting, workflow managers (Snakemake, Nextflow), or quick exploration.

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
| `pycdr results` | Export results table to CSV/TSV or stdout |
| `pycdr info` | Display dataset metadata |

### Quick start

```bash
# Minimal — run CDR-g analysis
pycdr run data.h5ad -p stim

# With gene filtering and enrichment
pycdr run data.h5ad -p stim -o results.h5ad -c results.csv \
  --filter-method numcells --min-cells 25 \
  --enrich --genecol gene_short_name

# Inspect a dataset
pycdr info data.h5ad -p stim

# Export top genes for a specific factor
pycdr results analyzed.h5ad --top-genes 20 --factor 3

# Re-run enrichment with different parameters
pycdr enrich analyzed.h5ad -p stim -m kruskal -o enriched.h5ad
```

### Common options for `pycdr run`

| Option | Default | Description |
|--------|---------|-------------|
| `-p, --phenotype` | required | Condition column in `adata.obs` |
| `-o, --output` | `{stem}_cdr.h5ad` | Output .h5ad path |
| `-c, --csv` | — | Export results CSV |
| `--filter-method` | `none` | `none`, `percent`, or `numcells` |
| `--enrich / --no-enrich` | off | Run enrichment after analysis |
| `--enrich-method` | `perm` | `perm` (chi-square) or `kruskal` (KW test) |
| `--genecol` | — | Gene name column (required for perm) |
| `-v / -vv` | — | INFO / DEBUG logging |

Run `pycdr run --help` for the full list of options.

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
