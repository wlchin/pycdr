[![wlchin](https://circleci.com/gh/wlchin/pycdr.svg?style=svg)](https://circleci.com/gh/circleci/circleci-docs)


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

# Example workflows

Three example snakemake workflows are provided in a separate [repository](https://github.com/wlchin/CDR_workflows). These workflows generate the results and describe preprocessing steps for each dataset in the manuscript. These CDR-g analyses use the visualisation and preprocessing functions provided in other single cell packages. 

To run the full workflows, please install [scanpy](https://scanpy-tutorials.readthedocs.io/en/latest/#), [bbknn](https://bbknn.readthedocs.io/en/latest/) (to allow dataset integration) and enrichment_utils (a simple wrapper around goatools to allow enrichment analysis on anndata objects analysed by CDR-g).

	pip install scanpy[leiden] bbknn enrichment_utils

# Documentation

CDR-g's documentation is available [here](http://cdr-g.readthedocs.io/)
