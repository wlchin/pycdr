#!/usr/bin/env python

from setuptools import setup

setup(name='pycdr',
      version='0.0.1',
      packages=['pycdr',
                'pycdr.test'],
      install_requires = ['anndata', 'scipy', 'dask-ml', 'tqdm', 'pandas', 'numpy==1.20.0', 'statsmodels'],
      python_requires='>3.5'
      )
