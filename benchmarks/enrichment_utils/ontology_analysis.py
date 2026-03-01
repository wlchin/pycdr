from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from . import mouse
from . import human
from goatools.associations import read_ncbi_gene2go

from pathlib import Path
import anndata as ad
import pandas as pd

import os
import sys
import logging
import time
import warnings

warnings.simplefilter("ignore")
logging.basicConfig(format='%(process)d - %(levelname)s : %(asctime)s - %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)


def get_datafile_name():
    """
    return the default data file that comes with the package
    """
    return Path(__file__).parent / "data/gene_index.csv"

def get_datafile_name_mouse():
    """
    return the default data file that comes with the package
    """
    return Path(__file__).parent / "data/mouse_index.csv"

def _create_human_go_object(godag_file, gene2go_file, propcounts):
    obodag = GODag(godag_file)
    objanno = Gene2GoReader(gene2go_file, taxids=True)
    from .human import GENEID2NT as GeneID2nt_human
    ns2assoc_h = objanno.get_ns2assc(9606)
    goeaobj_h = GOEnrichmentStudyNS(
        GeneID2nt_human.keys(), # List of mouse protein-coding genes
        ns2assoc_h,
        obodag, # Ontologies
        propagate_counts = propcounts,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method
    return goeaobj_h

def _create_mouse_go_object(godag_file, gene2go_file, propcounts):
    obodag = GODag(godag_file)
    objanno = Gene2GoReader(gene2go_file, taxids=True)
    ns2assoc_h = objanno.get_ns2assc(10090)
    from .mouse import GENEID2NT as GeneID2nt_mouse
    goeaobj_m = GOEnrichmentStudyNS(
        GeneID2nt_mouse.keys(), # List of mouse protein-coding genes
        ns2assoc_h,
        obodag, # Ontologies
        propagate_counts = propcounts,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method
    return goeaobj_m

def use_id_map(testlist, species):
    import pandas as pd
    if species == "human":
        df = pd.read_csv(get_datafile_name())
    if species == "mouse":
        df = pd.read_csv(get_datafile_name_mouse())
    df.geneid = df.geneid.astype(int)
    c = [x.upper() for x in testlist]
    index_mask = df['symbol'].isin(c)
    converted_list = df['geneid'][index_mask].to_list()

    return converted_list

def reverse_id_map(testlist_ncbi, species):
    import pandas as pd
    if species == "human":
        df = pd.read_csv(get_datafile_name())
    if species == "mouse":
        df = pd.read_csv(get_datafile_name_mouse())
    df.geneid = df.geneid.astype(int)
    index_mask = df['geneid'].isin(testlist_ncbi)
    dfs = df[index_mask]

    return dfs

def test_geneset(goeaobj, testlist, ontology_subset, threshold_pvalue):
    goea_results_all = goeaobj.run_study(testlist, prt=None)
    goea_results_sig = [r.name for r in goea_results_all if r.p_fdr_bh < threshold_pvalue and r.NS == ontology_subset]
    return goea_results_all, goea_results_sig

# functions for GOATOOLS enrichment


def get_terms_one_goea(oneentry, sp):
    goid = oneentry.GO
    goname = oneentry.name
    gotype = oneentry.NS
    hits = reverse_id_map(oneentry.kws["study_items"], sp).symbol.to_list()
    go_pval = oneentry.kws["p_uncorrected"]
    go_fdr = oneentry.p_fdr_bh
    oneid = {"goid" : goid, "goname" : goname, "gotype" : gotype, "hits" : hits, "go_pval" : go_pval, "fdr" : go_fdr}
    return goname, oneid

def create_list_for_append(goaobjlist, sp, ontology_subset, threshold_pvalue):
    # just do so for a single factor loading
    dict_int = {}
    for i in goaobjlist:
        if i.p_fdr_bh < threshold_pvalue and i.NS == ontology_subset:
            a, b = get_terms_one_goea(i, sp)
            dict_int[a] = b ## dict of dicts
    return dict_int

### helper function

def analyse_list(input_list,
                          dag,
                          g2g,
                          species,
                          threshold_pvalue = 0.05,
                          ontology_subset = "BP",
                          prop = False):

    if species == "human":
        goeaobj = _create_human_go_object(dag, g2g, prop)
    if species == "mouse":
        goeaobj = _create_mouse_go_object(dag, g2g, prop)

    genelist = input_list
    testlist = use_id_map(genelist, species)
    goa_obj, terms = test_geneset(goeaobj, testlist, ontology_subset, threshold_pvalue)

    return goa_obj, terms

def analyse_adata(adata, dag, g2g, species, threshold_pvalue = 0.05, ontology_subset = "BP", prop = False):
    # store terms in anndata

    """
    returns enrichment dict of terms and the dictionary of goaterm
    """
    #adata = ad.read(adata_f)

    import tqdm

    if species == "human":
        goeaobj = _create_human_go_object(dag, g2g, prop)
    if species == "mouse":
        goeaobj = _create_mouse_go_object(dag, g2g, prop)

    genelist = adata.var.index.to_list()

    logger.info('initialization complete')

    factor_dict = adata.uns["factor_loadings"]

    #print(factor_dict)

    enrichment_dict = {}

    goatools_dict = {}
    start = time.time()

    for k, v in tqdm.tqdm(factor_dict.items()):
        testlist = use_id_map(v, species)
        goa_obj, terms = test_geneset(goeaobj, testlist, ontology_subset, threshold_pvalue)
        enrichment_dict[k] = terms

        goatools_dict[k] = create_list_for_append(goa_obj, species, ontology_subset, threshold_pvalue)
    #em = get_sig_df(goatools_dict, 0.05)
    end = time.time()
    timediff = end - start

    logger.info('wall clock time in seconds:: %s', timediff)

    adata.uns["enrichment_results"] = enrichment_dict
    adata.uns["enrichment_details"] = goatools_dict
    logger.info('enrichment complete')
