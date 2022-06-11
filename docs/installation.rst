.. CDR-g documentation master file, created by
   sphinx-quickstart on Fri May 20 10:42:58 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Installation
============

CDR-g runs on python 3 and has been tested on python versions 3.5 through 3.9. It is strongly recommended that the installation is performed in a virtual environment. CDR-g is available on pyPI via:
	
.. code-block::

	pip install cdr-py

For enrichment of gene sets, we provide the enrichment_utils package, a simple wrapper around the python `goatools <https://github.com/tanghaibao/goatools>`_ package.

.. code-block::

    pip install enrichment_utils

CDR-g does not require `scanpy <https://scanpy-tutorials.readthedocs.io/en/latest/>`_ as a dependency but integrates well with it. A `scanpy <https://scanpy-tutorials.readthedocs.io/en/latest/>`_ installation is therefore recommended for analysis.

Docker
------

`Here <https://hub.docker.com/repository/docker/wlc27/pycdr_jupyter>`_ is a link to a image with CDR-g, enrichment_utils, and scanpy. It has jupyterlab installed. It is built directly off the `gcfntnu/scanpy <https://hub.docker.com/r/gcfntnu/scanpy>`_ image. 