# CDR

This repository contains the codebase for the CDR algorithm described in the manuscript **"Spectral detection of condition-specific biological pathways in single-cell gene expression data"**.

# Installation

CDR runs on python 3 and has been tested on python versions 3.5 through 3.9. It is strongly recommended that the installation is performed in a virtual environment. The most recent installation can be obtained by cloning the repository, navigating to the folder containing the setup.py file and installing via:

	pip install .
	

# Usage

The basic workflow is demonstrated below. As input, CDR requires a pre-prepared anndata object. Genes can be filtered based on variance or count criteria to reduce computation time. The condition of interest should be a column in the anndata obs. dataframe. CDR uses data from the count matrix (X) to construct co-expression matrices. Count data should be log-transformed. 

	run_CDR_analysis(anndata_object, condition_of_interest)
	calculate_enrichment(anndata_object)

# Example workflows

Three example snakemake workflows are provided in a separate [repository](https://github.com/wlchin/CDR_workflows). These workflows generate the results and describe preprocessing steps for each dataset in the manuscript. These CDR analyses use the visualisation and preprocessing functions provided in other packages. 

To run the full workflows, please install [scanpy](https://scanpy-tutorials.readthedocs.io/en/latest/#), [bbknn](https://bbknn.readthedocs.io/en/latest/) (to allow dataset integration) and enrichment_utils (a simple wrapper around goatools to allow enrichment analysis on anndata objects analysed by CDR).

	pip install scanpy[leiden] bbknn enrichment_utils


