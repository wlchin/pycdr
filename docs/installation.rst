.. CDR-g documentation master file, created by
   sphinx-quickstart on Fri May 20 10:42:58 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Installation
============

CDR-g runs on python 3 and has been tested on python versions 3.5 through 3.9. It is strongly recommended that the installation is performed in a virtual environment. CDR is available on pyPI via:
	
.. code-block::python

	pip install cdr-py

For enrichment of gene sets, one can also install the enrichment_utils package:

.. code-block::python

    pip install enrichment_utils

CDR-g does not require scanpy as a dependency but integrates with a scanpy workflow, so a scanpy installation is recommended for analysis.