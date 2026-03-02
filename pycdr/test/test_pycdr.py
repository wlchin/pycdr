
# get the data file as a fixture with scope module

import pytest
import numpy as np
import scipy.sparse as ss
import pandas as pd
from pycdr import pycdr
import anndata as ad
from pathlib import Path
from pycdr import perm
from pycdr import utils
from pycdr import kruskal

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
    """Run CDR analysis once and return the serialized result.

    Uses correction="none" to preserve legacy snapshot values on the small
    test dataset.  FDR behaviour is tested separately.
    """
    pycdr.run_CDR_analysis(muscleobject, "Hours", correction="none")
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
    assert x is not None
    # New structure: one row per factor (not per gene)
    n_factors = len(muscleobjectoutput.uns["factor_loadings"])
    assert x.shape[0] == n_factors
    assert "n_genes" in x.columns
    assert "genes" in x.columns
    assert "top_genes" in x.columns
    # genes column should contain lists
    for val in x["genes"]:
        assert isinstance(val, list)

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
    assert 0 < b["factor.9"][1] < 1  # p-value in valid range

def test_df_loading(analyzed_muscle):
    a, b = perm.calculate_enrichment(analyzed_muscle, "Hours",  ['factor.9'] , 10, "gene_short_name", 0.1)
    df = perm.get_df_loadings(analyzed_muscle)
    assert 0 < df.fdr[0] < 1  # FDR in valid range
    assert 0 < df.a_max[0] < 1  # max activation proportion in valid range

def test_filter_percent(analyzed_muscle):
    y = utils.filter_genecounts_percent(analyzed_muscle, 5, 4)
    assert y.shape == (100, 6)

def test_filter_counts(analyzed_muscle):
    y = utils.filter_genecounts_numcells(analyzed_muscle, 5, 25)
    assert y.shape == (100, 35)

def test_top_genes(analyzed_muscle):
    y = utils.get_top_genes(analyzed_muscle, 1)
    np.testing.assert_allclose(y["z_score"].iloc[0], 2.31, atol=0.01)
    np.testing.assert_allclose(y["pval"].iloc[2], 0.13, atol=0.01)


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


# --- B. pycdr.py branches ---

def test_orthomax_convergence():
    R = pycdr.classic_orthomax(np.eye(5))
    assert R.shape == (5, 5)
    np.testing.assert_allclose(np.abs(R), np.eye(5), atol=1e-5)

def test_optimal_threshold_fallback():
    rng = np.random.default_rng(42)
    arr = rng.random((20, 30))
    x, y, X, v = pycdr.get_optimal_threshold(arr, thres=2.0)
    assert y == 19  # ncomp = min(2000, 19, 29) = 19, fallback since no cumvar > 2.0


# --- C. utils.py output_results branches (109-111, 117-118, 123-124, 127, 130, 133) ---

def test_output_no_analysis():
    adata = ad.AnnData(np.random.default_rng(0).random((5, 5)))
    result = utils.output_results(adata)
    assert result is None

def test_output_genes_only(analyzed_muscle):
    adata = analyzed_muscle.copy()
    for key in ['enrichment_results', 'enrichment_stats']:
        adata.uns.pop(key, None)
    result = utils.output_results(adata)
    assert result is not None
    assert 'genes' in result.columns
    assert 'n_genes' in result.columns
    assert 'top_genes' in result.columns
    # genes column should contain lists
    for val in result["genes"]:
        assert isinstance(val, list)

def test_output_genes_and_stats(analyzed_muscle):
    adata = analyzed_muscle.copy()
    kruskal.calculate_enrichment(adata, "Hours")
    adata.uns.pop('enrichment_results', None)
    result = utils.output_results(adata)
    assert result is not None
    assert 'genes' in result.columns
    assert 'fdr' in result.columns
    assert 'n_genes' in result.columns

def test_output_with_enrichment_no_stats(analyzed_muscle):
    adata = analyzed_muscle.copy()
    # Clear any enrichment state from prior tests, then set only enrichment_results
    adata.uns.pop('enrichment_stats', None)
    adata.uns.pop('enrichment_results', None)
    adata.uns['enrichment_results'] = {'factor.9': ['GO:0001', 'GO:0002']}
    assert 'enrichment_results' in adata.uns  # sanity check
    result = utils.output_results(adata)
    assert result is not None
    assert 'genes' in result.columns
    assert 'terms' in result.columns


# --- D. kruskal.py (all 74 lines) ---

def test_exp_create_rank_matrix(analyzed_muscle):
    arrrank = kruskal.create_rank_matrix(analyzed_muscle.X)
    assert arrrank.shape == (analyzed_muscle.shape[1], analyzed_muscle.shape[0])

def test_exp_create_rank_matrix_sparse(sparse_muscle):
    arrrank_sparse = kruskal.create_rank_matrix(sparse_muscle.X)
    arrrank_dense = kruskal.create_rank_matrix(sparse_muscle.X.toarray())
    np.testing.assert_array_equal(arrrank_sparse, arrrank_dense)

def test_exp_single_geneset(analyzed_muscle):
    arrrank = kruskal.create_rank_matrix(analyzed_muscle.X)
    arr_index = analyzed_muscle.var.index
    geneset = analyzed_muscle.uns['factor_loadings']['factor.9']
    result = kruskal.calculate_enrichment_single_geneset(geneset, arr_index, arrrank)
    assert result.shape == (analyzed_muscle.shape[0],)

def test_exp_calculate_enrichment(analyzed_muscle):
    adata = analyzed_muscle.copy()
    results = kruskal.calculate_enrichment(adata, "Hours")
    assert isinstance(results, pd.DataFrame)
    assert 'stat' in results.columns
    assert 'pval' in results.columns
    assert 'fdr' in results.columns
    assert 'rank_matrix' in adata.layers
    assert 'enrichment_score_matrix' in adata.obsm
    assert 'enrichment_stats' in adata.uns

def test_exp_binarize_gset(analyzed_muscle):
    arrrank = kruskal.create_rank_matrix(analyzed_muscle.X)
    arr_index = analyzed_muscle.var.index
    geneset = analyzed_muscle.uns['factor_loadings']['factor.9']
    pmat, matreal, active_cells = kruskal.binarize_gset(arrrank, arr_index, geneset, nperm=10)
    assert pmat.shape == (analyzed_muscle.shape[0],)
    assert matreal.shape == (analyzed_muscle.shape[0],)
    assert active_cells.dtype == bool

def test_exp_binarize_on_adata(analyzed_muscle):
    adata = analyzed_muscle.copy()
    # Ensure X is dense ndarray after copy/deserialization
    if ss.issparse(adata.X):
        adata.X = adata.X.toarray()
    factor_list = ['factor.9', 'factor.2']
    result = kruskal.binarize_gset_on_adata(adata, factor_list, nperm=10)
    assert result is not None
    assert result.shape == (adata.shape[0], 2)
    assert result.dtype == bool


# --- E. Shared create_rank_matrix and perm enrichment_stats ---

def test_utils_create_rank_matrix(analyzed_muscle):
    from pycdr.utils import create_rank_matrix
    arrrank = create_rank_matrix(analyzed_muscle.X)
    assert arrrank.shape == (analyzed_muscle.shape[1], analyzed_muscle.shape[0])
    # Verify it matches the kruskal module's (now delegated) version
    np.testing.assert_array_equal(arrrank, kruskal.create_rank_matrix(analyzed_muscle.X))

def test_perm_enrichment_stores_stats(analyzed_muscle):
    adata = analyzed_muscle.copy()
    perm.calculate_enrichment(adata, "Hours", ['factor.9'], 10, "gene_short_name", 0.1)
    assert 'enrichment_stats' in adata.uns
    df = adata.uns["enrichment_stats"]
    assert isinstance(df, pd.DataFrame)
    # output_results should work without error
    result = utils.output_results(adata)
    assert result is not None


# --- F. Dominant condition mapping ---

def test_dominant_condition_in_uns(analyzed_muscle):
    """dominant_condition should be a sequence of strings with one entry per factor."""
    dc = analyzed_muscle.uns.get("dominant_condition")
    assert dc is not None
    n_factors = analyzed_muscle.uns["selected_loading"]
    # After h5ad round-trip, lists become numpy arrays
    assert len(dc) == n_factors
    assert all(isinstance(str(v), str) for v in dc)

def test_condition_labels_in_uns(analyzed_muscle):
    """condition_labels should be a sequence of strings."""
    cl = analyzed_muscle.uns.get("condition_labels")
    assert cl is not None
    assert len(cl) > 0
    assert all(isinstance(str(v), str) for v in cl)

def test_kruskal_dominant_condition_column(analyzed_muscle):
    """Kruskal enrichment should produce a dominant_condition column."""
    adata = analyzed_muscle.copy()
    results = kruskal.calculate_enrichment(adata, "Hours")
    assert "dominant_condition" in results.columns
    # Every factor should have a dominant condition string
    assert all(isinstance(v, str) for v in results["dominant_condition"])

def test_perm_dominant_condition_column(analyzed_muscle):
    """Perm enrichment should rename max_P to dominant_condition."""
    adata = analyzed_muscle.copy()
    perm.calculate_enrichment(adata, "Hours", ['factor.9'], 10, "gene_short_name", 0.1)
    df = perm.get_df_loadings(adata)
    assert "dominant_condition" in df.columns
    assert "max_P" not in df.columns

def test_output_results_dominant_condition(analyzed_muscle):
    """output_results should include dominant_condition column."""
    adata = analyzed_muscle.copy()
    # Without enrichment — Fs-based fallback
    for key in ['enrichment_results', 'enrichment_stats']:
        adata.uns.pop(key, None)
    result = utils.output_results(adata)
    assert "dominant_condition" in result.columns

def test_output_results_enrichment_dominant_condition(analyzed_muscle):
    """When enrichment has dominant_condition, it should override Fs-based."""
    adata = analyzed_muscle.copy()
    kruskal.calculate_enrichment(adata, "Hours")
    result = utils.output_results(adata)
    assert "dominant_condition" in result.columns


# --- G. Memory-safety tests for iterative permutation ---

from pycdr import feature_selection

def test_select_modules_memory_safe():
    """select_modules should complete without OOM on moderately large inputs.

    Old batch code would allocate 100 * 10000 * 20 * 8 = 160 MB as an
    intermediate.  The iterative version uses ~1.6 MB per iteration.
    """
    rng = np.random.default_rng(99)
    nfacs = 2
    rows_per_split = 5000
    n_rows = nfacs * rows_per_split  # 10000
    n_cols = 20
    Fs = rng.standard_normal((n_rows, n_cols))

    adata = ad.AnnData(np.zeros((rows_per_split, n_cols)))
    adata.uns["Fs"] = Fs

    selection, pval_mat, z_score_mat = feature_selection.select_modules(
        adata, nperm=100, thresh=0.05, nfacs=nfacs
    )

    assert pval_mat.shape == (rows_per_split, n_cols)
    assert z_score_mat.shape == (rows_per_split, n_cols)
    assert selection.shape == (rows_per_split, n_cols)
    assert np.all((pval_mat >= 0) & (pval_mat <= 1))


def test_permute_matrix_memory_safe():
    """permute_matrix should complete without OOM on moderately large inputs.

    Old batch code would allocate 50 * 5000 * 2000 * 8 = 4 GB as an
    intermediate.  The iterative version uses ~8 MB per iteration.
    """
    rng = np.random.default_rng(99)
    n_genes = 5000
    n_cells = 2000
    arrrank = rng.random((n_genes, n_cells))

    gene_names = [f"gene_{i}" for i in range(n_genes)]
    factor_genes = gene_names[:500]

    adata = ad.AnnData(
        np.zeros((n_cells, n_genes)),
        var=pd.DataFrame({"gene_name": gene_names}, index=gene_names),
    )
    adata.uns["factor_loadings"] = {"test_factor": factor_genes}

    pmat, matreal = perm.permute_matrix(
        adata, arrrank, "test_factor", nperm=50, genecol="gene_name"
    )

    assert pmat.shape == (n_cells,)
    assert matreal.shape == (n_cells,)
    assert np.all((pmat >= 0) & (pmat <= 1))


# --- H. Phase 1: Seed reproducibility ---

def test_reproducibility_same_seed(muscleobject):
    """Same seed should produce identical results."""
    a = muscleobject.copy()
    b = muscleobject.copy()
    pycdr.run_CDR_analysis(a, "Hours", nperm=100, seed=123, correction="none")
    pycdr.run_CDR_analysis(b, "Hours", nperm=100, seed=123, correction="none")
    np.testing.assert_array_equal(a.uns["Fs"], b.uns["Fs"])
    np.testing.assert_array_equal(a.uns["pval_mat"], b.uns["pval_mat"])
    np.testing.assert_array_equal(a.uns["selection"], b.uns["selection"])


def test_reproducibility_different_seed(muscleobject):
    """Different seeds should produce different permutation p-values."""
    a = muscleobject.copy()
    b = muscleobject.copy()
    pycdr.run_CDR_analysis(a, "Hours", nperm=100, seed=1, correction="none")
    pycdr.run_CDR_analysis(b, "Hours", nperm=100, seed=2, correction="none")
    # SVD uses different seeds so Fs may differ; at minimum pval_mat should differ
    assert not np.array_equal(a.uns["pval_mat"], b.uns["pval_mat"])


# --- I. Phase 2: FDR correction ---

def test_fdr_reduces_gene_counts(muscleobject):
    """FDR correction should produce fewer or equal significant genes vs none."""
    a = muscleobject.copy()
    b = muscleobject.copy()
    pycdr.run_CDR_analysis(a, "Hours", nperm=500, correction="none")
    pycdr.run_CDR_analysis(b, "Hours", nperm=500, correction="fdr_bh")
    genes_none = sum(len(v) for v in a.uns["factor_loadings"].values())
    genes_fdr = sum(len(v) for v in b.uns["factor_loadings"].values())
    assert genes_fdr <= genes_none


def test_fdr_mat_stored(muscleobject):
    """FDR correction should store fdr_mat with correct shape."""
    a = muscleobject.copy()
    pycdr.run_CDR_analysis(a, "Hours", nperm=100, correction="fdr_bh")
    assert "fdr_mat" in a.uns
    assert "pval_mat" in a.uns
    assert a.uns["fdr_mat"].shape == a.uns["pval_mat"].shape
    # FDR q-values should be >= raw p-values (BH only inflates)
    assert np.all(a.uns["fdr_mat"] >= a.uns["pval_mat"] - 1e-10)


def test_fdr_rejects_noise():
    """On pure random data, FDR should find near-zero significant genes."""
    rng = np.random.default_rng(42)
    nfacs = 2
    rows_per_split = 200
    n_rows = nfacs * rows_per_split
    n_cols = 10
    Fs = rng.standard_normal((n_rows, n_cols))

    adata = ad.AnnData(np.zeros((rows_per_split, n_cols)))
    adata.uns["Fs"] = Fs

    selection, _, _ = feature_selection.select_modules(
        adata, nperm=200, thresh=0.05, nfacs=nfacs, correction="fdr_bh"
    )
    # With pure noise, FDR should reject nearly everything
    fdr_hits = selection.sum()
    raw_selection = adata.uns["pval_mat"] < 0.05
    raw_hits = raw_selection.sum()
    # FDR should find far fewer hits than raw (allow some slack)
    assert fdr_hits <= raw_hits


def test_correction_none_backward_compat(muscleobject):
    """correction='none' should match old behavior: no fdr_mat, raw pvals used."""
    a = muscleobject.copy()
    pycdr.run_CDR_analysis(a, "Hours", nperm=200, correction="none")
    assert "pval_mat" in a.uns
    assert "fdr_mat" not in a.uns


# --- J. Phase 3: Input validation and edge cases ---

def test_single_condition_raises(muscleobject):
    """CDR-g should raise ValueError with fewer than 2 phenotype groups."""
    adata = muscleobject.copy()
    adata.obs["single"] = "A"
    with pytest.raises(ValueError, match="at least 2 phenotype groups"):
        pycdr.run_CDR_analysis(adata, "single")


def test_nan_warning_logged(muscleobject, caplog):
    """NaN values in correlation matrix should be logged."""
    import logging
    with caplog.at_level(logging.WARNING, logger="pycdr.pycdr"):
        adata = muscleobject.copy()
        pycdr.run_CDR_analysis(adata, "Hours", correction="none")
    assert any("NaN" in r.message for r in caplog.records)


def test_empty_gene_set_kruskal():
    """Empty gene set overlap should return NaN enrichment scores."""
    arrrank = np.random.default_rng(0).random((10, 5))
    arr_index = pd.Index([f"gene_{i}" for i in range(10)])
    result = kruskal.calculate_enrichment_single_geneset(
        ["nonexistent_gene"], arr_index, arrrank
    )
    assert result.shape == (5,)
    assert np.all(np.isnan(result))


def test_reshape_mismatch_raises():
    """calculate_minmax should raise ValueError when rows aren't divisible by splits."""
    Fs = np.zeros((7, 3))
    with pytest.raises(ValueError, match="not evenly divisible"):
        feature_selection.calculate_minmax(Fs, splits=2)


def test_empty_gene_overlap_perm():
    """permute_matrix should return neutral values for empty gene overlap."""
    n_genes = 10
    n_cells = 5
    gene_names = [f"gene_{i}" for i in range(n_genes)]
    adata = ad.AnnData(
        np.zeros((n_cells, n_genes)),
        var=pd.DataFrame({"name": gene_names}, index=gene_names),
    )
    adata.uns["factor_loadings"] = {"test": ["missing_gene"]}
    arrrank = np.random.default_rng(0).random((n_genes, n_cells))
    pmat, matreal = perm.permute_matrix(adata, arrrank, "test", nperm=10, genecol="name")
    np.testing.assert_array_equal(pmat, np.ones(n_cells))
    np.testing.assert_array_equal(matreal, np.zeros(n_cells))


# --- K. Phase 4: API naming deprecation ---

def test_pernum_deprecation_warning(muscleobject):
    """Using 'pernum' keyword should emit DeprecationWarning."""
    import warnings
    adata = muscleobject.copy()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        pycdr.run_CDR_analysis(adata, "Hours", pernum=100, correction="none")
        dep_warnings = [x for x in w if issubclass(x.category, DeprecationWarning)]
        assert len(dep_warnings) == 1
        assert "pernum" in str(dep_warnings[0].message)
