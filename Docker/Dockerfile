FROM gcfntnu/scanpy:1.7.0

RUN mamba install -y -q jupyterlab snakemake bbknn

RUN pip install cdr-py enrichment_utils

RUN mkdir /tmp/work

WORKDIR /tmp/work

CMD jupyter-lab --allow-root --port 9999 --no-browser 