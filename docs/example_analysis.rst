Example analysis
================

Introduction
------------

Here, we analyse the a subset of the human peripheral blood dataset provided in Kang et al. This is a dataset of human peripheral blood cells, divided into an interferon beta stimulated dataset and one without. Download the dataset.

.. code-block::

	! wget https://github.com/wlchin/CDR_workflows/raw/main/monocytes/resources/raw_monocyte_CD14.h5ad

Data preprossing and cleaning
-----------------------------

.. code-block:: python

    import scanpy as sc
    from pycdr.utils import filter_genecounts_numcells, filter_genecounts_percent

    mono = sc.read("raw_monocyte_CD1.h5ad")

    mono = filter_genecounts_percent(mono, 0.01, 1)
    mono = filter_genecounts_numcells(mono, 0, 100)

    sc.pp.log1p(mono)
    sc.pp.scale(mono, zero_center = False)


Running a CDR analysis
----------------------

.. code-block:: python

    from pycdr.pycdr import run_CDR_analysis
    from pycdr.perm import calculate_enrichment

    mono = ad.read(INPUT)
    run_CDR_analysis(mono, "stim")

Running enrichment
------------------




.. code-block:: python

    from enrichment_utils.ontology_analysis import analyse_adata
    from pycdr.perm import calculate_enrichment

    INPUT = snakemake.input[0]
    INPUT_FILE_GENE2GO = snakemake.input[1]
    INPUT_FILE_ONTOLOGY = snakemake.input[2]

    OUTPUT = snakemake.output[0]

    analyse_adata(mono, INPUT_FILE_ONTOLOGY, INPUT_FILE_GENE2GO, "human", threshold_pvalue = 0.05, ontology_subset = "BP", prop = False)
    factor_list = [i for i in mono.uns["factor_loadings"].keys()]
    calculate_enrichment(mono, "stim", factor_list, 100, "features", 0.05)
