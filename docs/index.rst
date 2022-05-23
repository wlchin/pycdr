.. CDR-g documentation master file, created by
   sphinx-quickstart on Fri May 20 10:42:58 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CDR-g's documentation!
=================================

CDR-g (CDR-genomics) is a tool for extracting condition-dependent sources of variation from single cell data. 

CDR-g extends the CDR algorithm [#f1]_ initially proposed by Dr Leonardo Portes des Santos and Dr Michael Small to large single cell datasets.

The source code for this package is available on `github <https://github.com/wlchin/pycdr>`_.

The CDR-g workflow, provided as a separate `repository <https://github.com/wlchin/CDR_workflows>_` demonstrates example use-cases of CDR-g and reproduces the results in the manuscript.

Motivation
==========

In complex single-cell RNA-sequencing experiments, a crucial analysis is to identify conditional sources of transcriptional variation. Differential expression and differential co-expression are well established strategies for analysing RNA seq data. However, the challenge in these types of analyses in single cell datasets is that:

- each condition may be composed of multiple "subconditions", which escape detection
- variation may be shared between conditions. 

CDR-g overcomes these limitations, allowing sensitive, condition-dependent variation to be performed in an unsupervised manner within a single analysis. 

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation.rst
   algorithm.rst
   getting_started.rst
   api.rst
   
Reference
=========

.. [#f1] Portes LL, Small M. Navigating differential structures in complex networks. Phys Rev E. 2020 Dec;102(6-1):062301. doi: 10.1103/PhysRevE.102.062301. PMID: 33466036.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
