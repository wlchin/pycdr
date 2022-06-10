Monocyte dataset 
================

Introduction
------------

Here, we analyse the a subset of the human peripheral blood dataset provided in Kang et al [fn1]_. This is a monocyte subset of a dataset of human peripheral blood cells, divided into an interferon-beta stimulated vs treatment naive population. We first download the dataset from the github repository:

.. code-block:: shell

	! wget https://github.com/wlchin/CDR_workflows/raw/main/monocytes/resources/raw_monocyte_CD14.h5ad


Data preprocessing
------------------

After loading the dataset, we implement several preprocessing steps: filtering, log-transform, and scaling step prior to running CDR-g. 

.. code-block:: python

    import scanpy as sc
    from pycdr.utils import filter_genecounts_numcells, filter_genecounts_percent

    mono = sc.read("raw_monocyte_CD1.h5ad")
    mono = filter_genecounts_percent(mono, 0.01, 1)
    mono = filter_genecounts_numcells(mono, 0, 100)

    sc.pp.log1p(mono)
    sc.pp.scale(mono, zero_center = False)


Running a CDR analysis
----------------------

This code-block runs the CDR analysis. The condition of interest should be provided as a column in the addndata.obs dataframe. In this case, we are interested in the "stim" column, which labels cells as either interferon-stimulated or treatment naive. 

Briefly, this CDR function wraps three critical steps: (1) The construction and concatenation of coexpression matrices for each condition (2) The singular value decomposition step and (3) the selection of important genes from each factor loading resulting from SVD using a permutation filter. 

.. code-block:: python

    from pycdr.pycdr import run_CDR_analysis

    run_CDR_analysis(mono, "stim")


Running gene set enrichment
---------------------------

The resulting gene expression programs extracted by CDR-g, which show variation between conditions, can be found in the unstructured annotation (anndata.uns) of the anndata object. These results are stored as a dictionary of lists, with each key corresponding to the variable genes found in each corresponding factor loading. 

To better understand these genes, We perform gene set enrichment on these gene sets. We use enrichment_utils, a simple wrapper around the goatools [fn2]_ package. Enrichment using goatools requires an ontology (in this case the PANTHER GO-SLIM ontology) and the NCBIs gene2go mapping, so we download these accordingly. 

.. code-block:: shell

    ! wget http://data.pantherdb.org/PANTHER17.0/ontology/PANTHERGOslim.obo
    ! wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
    ! tar -xzvf gene2go.gz


We run the ontology analysis with the code block below. We examine only enriched terms from the GO-Biological Processes subset of the ontology in humans.  

.. code-block:: python
    
    from enrichment_utils.ontology_analysis import analyse_adata

    INPUT_FILE_GENE2GO = "PANTHERGOslim.obo"
    INPUT_FILE_ONTOLOGY = "gene2go"

    analyse_adata(mono, INPUT_FILE_ONTOLOGY, INPUT_FILE_GENE2GO, "human", ontology_subset = "BP")
    

Comparing gene set activation between condition
-----------------------------------------------

The final stage of the analysis is to identify gene sets which are more activated between conditions of interest. We have implemented a simple test of proportions that compares the number of cells with "activated gene set" in each condition, which we calculate gene set activation using ssGSEA [fn3]_. Below, we provide all factors as a list and calculate whether a gene set is activated based on a permutation test, thresholded at 0.05.

.. code-block:: python

    from pycdr.perm import calculate_enrichment

    factor_list = [i for i in mono.uns["factor_loadings"].keys()]
    calculate_enrichment(mono, "stim", factor_list, 100, "features", 0.05)


References
----------

.. [fn1] Kang, H. M., Subramaniam, M., Targ, S., Nguyen, M., Maliskova, L., McCarthy, E., Wan, E., Wong, S., Byrnes, L., Lanata, C. M., Gate, R. E., Mostafavi, S., Marson, A., Zaitlen, N., Criswell, L. A., & Ye, C. J. (2018). Multiplexed droplet single-cell RNA-sequencing using natural genetic variation. Nature biotechnology, 36(1), 89–94. https://doi.org/10.1038/nbt.4042

.. [fn2] Foroutan, M., Bhuva, D. D., Lyu, R., Horan, K., Cursons, J., & Davis, M. J. (2018). Single sample scoring of molecular phenotypes. BMC bioinformatics, 19(1), 404. https://doi.org/10.1186/s12859-018-2435-4

.. [fn3] Klopfenstein, D. V., Zhang, L., Pedersen, B. S., Ramírez, F., Warwick Vesztrocy, A., Naldi, A., Mungall, C. J., Yunes, J. M., Botvinnik, O., Weigel, M., Dampier, W., Dessimoz, C., Flick, P., & Tang, H. (2018). GOATOOLS: A Python library for Gene Ontology analyses. Scientific reports, 8(1), 10872. https://doi.org/10.1038/s41598-018-28948-z