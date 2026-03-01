
# get the data file as a fixture with scope module

import pytest
import numpy as np
import scipy.sparse as ss
import dask.array as da
import pandas as pd
from pycdr import pycdr
import anndata as ad
from pathlib import Path
from pycdr import perm
from pycdr import utils
from pycdr import experimental

@pytest.fixture(scope='module')
def muscleobject():
    """
    return the default data file that comes with the package
    """
    adatapath = Path(__file__).parents[1] / "data/test_muscle.h5ad"
    x = ad.read_h5ad(adatapath)
    return x

@pytest.fixture(scope='module')
def muscleobjectoutput():
    """
    return the default data file that comes with the package
    """
    adatapath = Path(__file__).parents[1] / "data/test_fl_e_s.h5ad"
    x = ad.read_h5ad(adatapath)
    return x

@pytest.fixture(scope='module')
def analyzed_muscle(muscleobject, tmp_path_factory):
    """Run CDR analysis once and return the serialized result."""
    pycdr.run_CDR_analysis(muscleobject, "Hours")
    file = tmp_path_factory.mktemp("data") / "output.h5ad"
    muscleobject.write(str(file))
    x = ad.read_h5ad(str(file))
    return x

@pytest.fixture(scope='module')
def sparse_muscle(muscleobject):
    """Return a copy of muscleobject with sparse X matrix."""
    adata = muscleobject.copy()
    adata.X = ss.csr_matrix(adata.X)
    return adata

def test_output_function(muscleobjectoutput):
    x = utils.output_results(muscleobjectoutput)
    assert x.shape == (1749, 9)

def test_CDR_muscle_adata_shape(analyzed_muscle):
    assert analyzed_muscle.shape[0] == 100

def test_CDR_muscle_Fs(analyzed_muscle):
    assert analyzed_muscle.uns["Fs"].shape == (200,22)

def test_CDR_muscle_loading_selection(analyzed_muscle):
    assert analyzed_muscle.uns["selected_loading"] == 22

def test_CDR_muscle_Fs_mat(analyzed_muscle):
    assert round(analyzed_muscle.uns["Fs"][25,3],4) == 0.1033

def test_CDR_muscle_Fs_shape_pval(analyzed_muscle):
    assert analyzed_muscle.uns["pval_mat"].shape == (100, 22)

def test_CDR_muscle_Fs_pval_values(analyzed_muscle):
    assert analyzed_muscle.uns["pval_mat"][15,10] == 0.4535

def test_CDR_muscle_Fs_factor_loading(analyzed_muscle):
    assert analyzed_muscle.uns['factor_loadings']['factor.2'][0] == 'ARHGAP33'

def test_CDR_muscle_Fs_diff(analyzed_muscle):
    assert round(analyzed_muscle.uns['Fs_diff'][23,6], 4) == 0.235

def test_enrichment(analyzed_muscle):
    a, b = perm.calculate_enrichment(analyzed_muscle, "Hours",  ['factor.9'] , 10, "gene_short_name", 0.1)
    assert a["factor.9"][0] == [5, 0]
    assert b["factor.9"][1] == 0.06659680506323108

def test_df_loading(analyzed_muscle):
    a, b = perm.calculate_enrichment(analyzed_muscle, "Hours",  ['factor.9'] , 10, "gene_short_name", 0.1)
    c = perm.get_df_loadings(analyzed_muscle).fdr[0]
    d = perm.get_df_loadings(analyzed_muscle).a_max[0]
    assert c == 0.06659680506323108
    assert d == 0.08196721311475409

def test_filter_percent(analyzed_muscle):
    y = utils.filter_genecounts_percent(analyzed_muscle, 5, 4)
    assert y.shape == (100, 6)

def test_filter_counts(analyzed_muscle):
    y = utils.filter_genecounts_numcells(analyzed_muscle, 5, 25)
    assert y.shape == (100, 35)

def test_top_genes(analyzed_muscle):
    y = utils.get_top_genes(analyzed_muscle, 1)
    assert round(y["z_score"][0],4) == 2.3099
    assert round(y["pval"][2], 4) == 0.1295


# --- A. Sparse data paths (pycdr:65, utils:25+52, perm:24) ---

def test_CDR_sparse_input(sparse_muscle):
    adata = sparse_muscle.copy()
    adata.X = ss.csr_matrix(adata.X)
    pycdr.run_CDR_analysis(adata, "Hours")
    assert adata.uns["Fs"].shape == (200, 22)

def test_filter_percent_sparse(sparse_muscle):
    adata = sparse_muscle.copy()
    adata.X = ss.csr_matrix(adata.X)
    y = utils.filter_genecounts_percent(adata, 5, 4)
    assert y.shape == (100, 6)

def test_filter_counts_sparse(sparse_muscle):
    adata = sparse_muscle.copy()
    adata.X = ss.csr_matrix(adata.X)
    y = utils.filter_genecounts_numcells(adata, 5, 25)
    assert y.shape == (100, 35)

def test_rank_matrix_sparse(sparse_muscle):
    adata = sparse_muscle.copy()
    adata.X = ss.csr_matrix(adata.X)
    arrrank = perm.create_rank_matrix(adata)
    assert arrrank.shape == (adata.shape[1], adata.shape[0])


# --- B. pycdr.py branches (103-104, 151, 197) ---

def test_get_numbers_of_pheno(analyzed_muscle):
    vals = pycdr.get_numbers_of_pheno(analyzed_muscle, "Hours")
    assert isinstance(vals, list)
    assert all(isinstance(v, (int, np.integer)) for v in vals)
    assert sum(vals) == analyzed_muscle.shape[0]

def test_orthomax_convergence():
    R = pycdr.classic_orthomax(np.eye(5))
    assert R.shape == (5, 5)
    np.testing.assert_allclose(np.abs(R), np.eye(5), atol=1e-5)

def test_optimal_threshold_fallback():
    rng = np.random.default_rng(42)
    arr = rng.random((20, 30))
    darr = da.from_array(arr, chunks=(20, 30))
    x, y, X, v = pycdr.get_optimal_threshold(darr, thres=2.0)
    assert y == 19  # ncomp = min(2000, 19, 29) = 19, fallback since no cumvar > 2.0


# --- C. utils.py output_results branches (109-111, 117-118, 123-124, 127, 130, 133) ---

def test_output_no_analysis():
    adata = ad.AnnData(np.random.default_rng(0).random((5, 5)))
    result = utils.output_results(adata)
    assert result is None

def test_output_genes_only(analyzed_muscle):
    adata = analyzed_muscle.copy()
    for key in ['enrichment_results', 'pval_dict', 'dict_res_prop']:
        adata.uns.pop(key, None)
    result = utils.output_results(adata)
    assert result is not None
    assert 'genes' in result.columns

def test_output_genes_and_stats(analyzed_muscle):
    adata = analyzed_muscle.copy()
    perm.calculate_enrichment(adata, "Hours", ['factor.9'], 10, "gene_short_name", 0.1)
    adata.uns.pop('enrichment_results', None)
    result = utils.output_results(adata)
    assert result is not None
    assert 'genes' in result.columns
    assert 'fdr' in result.columns

def test_output_with_enrichment_no_stats(analyzed_muscle):
    adata = analyzed_muscle.copy()
    adata.uns['enrichment_results'] = {'factor.9': ['GO:0001', 'GO:0002']}
    for key in ['pval_dict', 'dict_res_prop']:
        adata.uns.pop(key, None)
    result = utils.output_results(adata)
    assert result is not None
    assert 'genes' in result.columns
    assert 'terms' in result.columns


# --- D. experimental.py (all 74 lines) ---

def test_exp_create_rank_matrix(analyzed_muscle):
    arrrank = experimental.create_rank_matrix(analyzed_muscle.X)
    assert arrrank.shape == (analyzed_muscle.shape[1], analyzed_muscle.shape[0])

def test_exp_create_rank_matrix_sparse(sparse_muscle):
    arrrank_sparse = experimental.create_rank_matrix(sparse_muscle.X)
    arrrank_dense = experimental.create_rank_matrix(sparse_muscle.X.toarray())
    np.testing.assert_array_equal(arrrank_sparse, arrrank_dense)

def test_exp_single_geneset(analyzed_muscle):
    arrrank = experimental.create_rank_matrix(analyzed_muscle.X)
    arr_index = analyzed_muscle.var.index
    geneset = analyzed_muscle.uns['factor_loadings']['factor.9']
    result = experimental.calculate_enrichment_single_geneset(geneset, arr_index, arrrank)
    assert result.shape == (analyzed_muscle.shape[0],)

def test_exp_calculate_enrichment(analyzed_muscle):
    adata = analyzed_muscle.copy()
    results = experimental.calculate_enrichment(adata, "Hours")
    assert isinstance(results, pd.DataFrame)
    assert 'stat' in results.columns
    assert 'pval' in results.columns
    assert 'fdr' in results.columns
    assert 'rank_matrix' in adata.layers
    assert 'enrichment_score_matrix' in adata.obsm

def test_exp_binarize_gset(analyzed_muscle):
    arrrank = experimental.create_rank_matrix(analyzed_muscle.X)
    arr_index = analyzed_muscle.var.index
    geneset = analyzed_muscle.uns['factor_loadings']['factor.9']
    pmat, matreal, active_cells = experimental.binarize_gset(arrrank, arr_index, geneset, nperm=10)
    assert pmat.shape == (analyzed_muscle.shape[0],)
    assert matreal.shape == (analyzed_muscle.shape[0],)
    assert active_cells.dtype == bool

def test_exp_binarize_on_adata(analyzed_muscle):
    adata = analyzed_muscle.copy()
    factor_list = ['factor.9', 'factor.2']
    result = experimental.binarize_gset_on_adata(adata, factor_list, nperm=10)
    assert result.shape == (adata.shape[0], 2)
    assert result.dtype == bool
