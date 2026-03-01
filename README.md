[![Tests](https://github.com/wlchin/pycdr/actions/workflows/test.yml/badge.svg)](https://github.com/wlchin/pycdr/actions/workflows/test.yml) [![Coverage](https://img.shields.io/badge/coverage-100%25-brightgreen)](https://github.com/wlchin/pycdr/actions/workflows/test.yml)


# CDR-g (CDR-genomics)

This repository contains the codebase for the CDR-g (CDR-genomics) algorithm described in the manuscript **"Spectral detection of condition-specific biological pathways in single-cell gene expression data"**.

# Installation

CDR-g runs on python 3 and has been tested on python versions 3.5 through 3.9. It is strongly recommended that the installation is performed in a virtual environment. CDR-g is available on pyPI via:
	
	pip install cdr-py

# Usage

The basic workflow is demonstrated below. As input, CDR-g requires a pre-prepared anndata object. Genes can be filtered based on variance or count criteria to reduce computation time. The condition of interest should be a column in the anndata obs. dataframe. CDR-g uses data from the count matrix (X) to construct co-expression matrices. Count data should be log-transformed. The two steps below will (1) run the CDR-g analysis to produce gene expression programs and (2) perform single cell enrichment on each gene expression program recovered by CDR-g.

	from pycdr.pycdr import run_CDR_analysis
	fom pycdr.perm import calculate_enrichment

	run_CDR_analysis(anndata_object, condition_of_interest)
	calculate_enrichment(anndata_object)

# Command-line interface

After installation, the `pycdr` CLI is available. It wraps the full CDR-g pipeline so you can run analyses directly from the terminal without writing Python code.

```
pycdr --help
```

## Commands

| Command | Description |
|---------|-------------|
| `pycdr run` | Full pipeline (filter + analyze + enrich + export) |
| `pycdr analyze` | Run CDR-g SVD/varimax analysis only |
| `pycdr filter` | Filter genes from an .h5ad file |
| `pycdr enrich` | Run enrichment on previously analyzed data |
| `pycdr results` | Export results table to CSV/TSV or stdout |
| `pycdr info` | Display dataset metadata |

## Quick start

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

## `pycdr run` options

```
Usage: pycdr run [OPTIONS] INPUT

Options:
  -p, --phenotype TEXT            Condition column in adata.obs [required]
  -o, --output TEXT               Output .h5ad path [default: {stem}_cdr.h5ad]
  -c, --csv TEXT                  Export results CSV
  --capvar FLOAT                  Variance threshold [default: 0.95]
  --pernum INTEGER                Permutations for importance [default: 2000]
  --thres FLOAT                   P-value threshold [default: 0.05]
  --filter-method [none|percent|numcells]
                                  Gene filter method [default: none]
  --cell-fraction FLOAT           (percent) fraction of cells [default: 0.05]
  --median-count FLOAT            (percent) median count threshold [default: 1.0]
  --count-threshold INTEGER       (numcells) count cutoff [default: 1]
  --min-cells INTEGER             (numcells) min expressed cells [default: 10]
  --enrich / --no-enrich          Run enrichment [default: no-enrich]
  --enrich-method [perm|kruskal]  Enrichment method [default: perm]
  --genecol TEXT                  Gene name column in adata.var
  --nperm INTEGER                 Permutations for ssGSEA [default: 100]
  --enrich-thresh FLOAT           Active gene set threshold [default: 0.05]
  --seed INTEGER                  Random seed [default: 42]
  -v, --verbose                   Increase verbosity (-v INFO, -vv DEBUG)
  -q, --quiet                     Errors only
```

# Example workflows

Three example snakemake workflows are provided in the `benchmarks/CDR_workflows` directory. These workflows generate the results and describe preprocessing steps for each dataset in the manuscript. These CDR-g analyses use the visualisation and preprocessing functions provided in other single cell packages.

To run the full workflows, please install [scanpy](https://scanpy-tutorials.readthedocs.io/en/latest/#), [bbknn](https://bbknn.readthedocs.io/en/latest/) (to allow dataset integration) and enrichment_utils (a simple wrapper around goatools to allow enrichment analysis on anndata objects analysed by CDR-g).

	pip install scanpy[leiden] bbknn enrichment_utils

# Documentation

CDR-g's documentation is available [here](http://cdr-g.readthedocs.io/)
