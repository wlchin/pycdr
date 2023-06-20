FROM gcfntnu/scanpy:1.7.0

RUN mamba install -y -q jupyterlab snakemake bbknn

RUN pip install cdr-py enrichment_utils

RUN mkdir /tmp/work

WORKDIR /tmp/work

RUN pip install --upgrade matplotlib

RUN pip install pandas==1.5.0

CMD jupyter-lab --allow-root --port=9999 --ip=0.0.0.0 --no-browser 
