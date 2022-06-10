Algorithm description
=================================

The CDR-g approach identifies differentially-expressed and differentially co-expressed genes in multi-condition data. 

.. image:: algorithm.png
   :width: 500px

The steps in the CDR-g algorithm involve:

1. concatenation of condition-specific gene co-expression matrices (A, B), 
2. and the decomposition of their product to produce the resulting factor loading matrix using truncated SVD (C),
3. followed by varimax rotation of the resulting factor loading matrix (D), 
4. Gene set selection (E) is performed on each rotated factor loading. Genes which distinguish conditions have high factor loading scores and are selected using a permutation filter. 

The red box describes the original CDR framework as described in [#f1]_ and the blue box describes CDR-g, the extension to multi-condition single cell data.


.. [#f1] Portes LL, Small M. Navigating differential structures in complex networks. Phys Rev E. 2020 Dec;102(6-1):062301. doi: 10.1103/PhysRevE.102.062301. PMID: 33466036.
