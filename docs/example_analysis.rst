.. _monocyte-example:

Monocyte dataset 
================

Introduction
------------

Here, we analyse a subset of the human peripheral blood dataset provided in Kang et al [#fn1]_. This is a CD14 monocyte dataset, divided into an interferon-beta stimulated vs. treatment naive population. We first download the dataset from the github repository:

.. code-block:: shell

	wget https://github.com/wlchin/CDR_workflows/raw/main/monocytes/resources/raw_monocyte_CD14.h5ad


Data pre-processing
------------------

After loading the dataset, we implement several pre-processing steps: filtering, a log-transform, and scaling prior to running CDR-g. Below, the filtering step removes genes based on two criteria: count threshold and number of cells in which the gene is expressed. 

.. code-block:: python

    import scanpy as sc
    from pycdr.utils import filter_genecounts_numcells, filter_genecounts_percent

    mono = sc.read("raw_monocyte_CD14.h5ad")
    mono = filter_genecounts_percent(mono, 0.01, 1)
    mono = filter_genecounts_numcells(mono, 0, 100)

    sc.pp.log1p(mono)
    sc.pp.scale(mono, zero_center = False)


Running a CDR-g analysis
----------------------

This code-block runs the CDR-g analysis. The condition of interest should be provided as a column in the anndata.obs dataframe. In this case, we are interested in the "stim" column, which labels cells as either interferon-stimulated or treatment naive. 

Briefly, this CDR-g function wraps three critical steps: (1) The construction and concatenation of coexpression matrices for each condition (2) The singular value decomposition (SVD) step and (3) the selection of important genes from each factor loading resulting from SVD using a permutation filter. 

.. code-block:: python

    from pycdr.pycdr import run_CDR_analysis

    run_CDR_analysis(mono, "stim")


Running gene set enrichment
---------------------------

.. note::
    This analysis is optional. You are welcome to use other gene enrichment analysis strategies.

The resulting gene expression programs extracted by CDR-g, which show variation between conditions, can be found in the unstructured annotation (anndata.uns) of the anndata object. These results are stored as a dictionary of lists, with each key corresponding to the variable genes found in each corresponding factor loading. 

To better understand these genes, We perform gene set enrichment on these gene sets. We provide enrichment_utils, a simple wrapper around the goatools [#fn2]_ package. Enrichment using goatools requires an ontology file (in this case the PANTHER GO-SLIM ontology) and the NCBI gene2go mapping, so we download these accordingly. 

.. code-block:: shell

    wget http://data.pantherdb.org/PANTHER17.0/ontology/PANTHERGOslim.obo
    wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
    tar -xzvf gene2go.gz


We run the ontology analysis with the code block below. We examine only enriched GO-terms from the biological processes subset of the ontology terms in humans.  

.. code-block:: python
    
    from enrichment_utils.ontology_analysis import analyse_adata

    INPUT_FILE_GENE2GO = "PANTHERGOslim.obo"
    INPUT_FILE_ONTOLOGY = "gene2go"

    analyse_adata(mono, INPUT_FILE_ONTOLOGY, INPUT_FILE_GENE2GO, "human", ontology_subset = "BP")
    

Comparing gene set activation between conditions
-----------------------------------------------

The final stage of the analysis is to identify gene sets which are more activated between conditions of interest. We have implemented a `test of proportions <https://www.statsmodels.org/devel/generated/statsmodels.stats.proportion.proportions_chisquare.html>`_ that compares the number of cells with the "activated gene set" in each condition. We calculate gene set activation using ssGSEA [#fn3]_. Below, we test all factors and calculate whether a gene set is activated based on a permutation test, thresholded at a pvalue of =<0.05.

.. code-block:: python

    from pycdr.perm import calculate_enrichment

    factor_list = [i for i in mono.uns["factor_loadings"].keys()]
    calculate_enrichment(mono, "stim", factor_list, 100, "features", 0.05)

Viewing and saving results
--------------------------

The results from the analysis can be viewed as a dataframe and persisted to disk.

.. code-block:: python

    res = output_results(mono)

Example output from the first row of the dataframe:

.. code-block:: python

    res.iloc[0,]


    genes                     RSAD2,IFIT3,IFIT1,ISG20,APOBEC3A,MX1
    terms        defense response to virus,type I interferon si...
    max_P                                                     STIM
    a_max                                                 0.098277
    a_range                                               0.098277
    a_mean                                                0.049138
    statistic                                           228.747876
    pvalue                                                     0.0
    fdr                                                        0.0
    Name: factor.0, dtype: object


.. [#fn1] Kang, H. M., Subramaniam, M., Targ, S., Nguyen, M., Maliskova, L., McCarthy, E., Wan, E., Wong, S., Byrnes, L., Lanata, C. M., Gate, R. E., Mostafavi, S., Marson, A., Zaitlen, N., Criswell, L. A., & Ye, C. J. (2018). Multiplexed droplet single-cell RNA-sequencing using natural genetic variation. Nature biotechnology, 36(1), 89–94. https://doi.org/10.1038/nbt.4042

.. [#fn3] Foroutan, M., Bhuva, D. D., Lyu, R., Horan, K., Cursons, J., & Davis, M. J. (2018). Single sample scoring of molecular phenotypes. BMC bioinformatics, 19(1), 404. https://doi.org/10.1186/s12859-018-2435-4

.. [#fn2] Klopfenstein, D. V., Zhang, L., Pedersen, B. S., Ramírez, F., Warwick Vesztrocy, A., Naldi, A., Mungall, C. J., Yunes, J. M., Botvinnik, O., Weigel, M., Dampier, W., Dessimoz, C., Flick, P., & Tang, H. (2018). GOATOOLS: A Python library for Gene Ontology analyses. Scientific reports, 8(1), 10872. https://doi.org/10.1038/s41598-018-28948-z