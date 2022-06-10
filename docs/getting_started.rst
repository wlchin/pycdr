.. CDR-g documentation master file, created by
   sphinx-quickstart on Fri May 20 10:42:58 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Usage
=====

Overview
--------

In this section, we provide a high-level overview of the CDR-g workflow. 

As input, CDR requires a pre-prepared anndata object. Genes can be filtered based on variance or count criteria to reduce computation time. The condition of interest should be a column in the anndata.obs dataframe. 

CDR uses data from the count matrix (anndata.X) to construct co-expression matrices. Count data should be scaled and log-transformed. 

The two steps below will (1) run the CDR analysis to produce gene expression programs and (2) perform single cell enrichment on each gene expression program recovered by CDR-g. The output of CDR-g is a dictionary of gene lists, with each list representing a gene expression program which varies between the conditions of interest.

.. code-block:: python

	from pycdr.pycdr import run_CDR_analysis
	from pycdr.perm import calculate_enrichment

	run_CDR_analysis(anndata_object, condition_of_interest)
	calculate_enrichment(anndata_object)


Example (snakemake) workflows
-----------------------------

Three example snakemake workflows are provided in a separate `repository <https://github.com/wlchin/CDR_workflows>`_. These workflows generate the results and describe preprocessing steps for each dataset in the manuscript. These CDR-g analyses use the visualisation and preprocessing functions provided in other single cell packages. 

To run the full workflows, please install `scanpy <https://scanpy-tutorials.readthedocs.io/en/latest/>`_, `bbknn <https://bbknn.readthedocs.io/en/latest/>`_ (to allow dataset integration) and enrichment_utils (a simple wrapper around `goatools <https://github.com/tanghaibao/goatools>`_ to allow enrichment analysis on anndata objects analysed by CDR-g).

.. code-block::

	pip install snakemake scanpy[leiden] bbknn enrichment_utils

.. warning::

    The workflows download large datasets from GEO. 
    

Walkthrough analysis
--------------------

An annotated example is provided here: `Monocyte dataset`_