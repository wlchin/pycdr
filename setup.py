#!/usr/bin/env python

from setuptools import setup

setup(name='cdr_py',
      version='0.0.2',
      packages=['pycdr',
                'pycdr.test'],
      install_requires = ['anndata==0.7.6', 'scipy', 'dask-ml', 'tqdm', 'pandas==1.3.5', 'numpy==1.20.0', 'statsmodels'],
      python_requires='>3.5',
      author='WL Chin',
      package_data={'': ['data/*.*']},
      author_email="melvin.chin@research.uwa.edu.au",
      url = "https://github.com/wlchin/pycdr",
      description = "CDR algorithm for multi-condition single-cell data"
      )
