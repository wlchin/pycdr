API
===

Import pycdr in a python session via::

  from pycdr import run_CDR_analysis, calculate_enrichment


Cell filtering
--------------

.. autosummary::
    :toctree: .

    pycdr.utils.filter_genecounts_percent
    pycdr.utils.filter_genecounts_numcells

CDR algorithm
-------------

.. autosummary::
    :toctree: .

    pycdr.pycdr.run_CDR_analysis

Enrichment analysis
-------------------

.. autosummary::
    :toctree: .

    pycdr.experimental.calculate_enrichment
    pycdr.experimental.binarize_gset_on_adata

Results and output
------------------

.. autosummary::
    :toctree: .

    pycdr.utils.output_results
    pycdr.utils.get_top_genes

Legacy (permutation-based enrichment)
-------------------------------------

.. autosummary::
    :toctree: .

    pycdr.perm.calculate_enrichment
    pycdr.perm.get_df_loadings
