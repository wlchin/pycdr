automodule:: pycdr

API
===

Import pycdr in a python session via::

  from pycdr.pycdr import run_CDR_analysis
  from pycdr.perm import calculate_enrichment

  
Cell filtering
--------------

**Filtering anndata objects** 

.. autosummary::
    :toctree: .

    pycdr.utils.filter_genecounts_percent
    pycdr.utils.filter_genecounts_numcells

CDR algorithm
-------------
    
**Decomposition**

SVD function 

.. autosummary::
    :toctree: .
	     
    pycdr.pycdr.run_CDR_analysis
    
**Permutation and thresholding**

Module specific functions for permutations

.. autosummary::
    :toctree: .

    pycdr.perm.create_rank_matrix

Enrichment
----------





Utility
-------

