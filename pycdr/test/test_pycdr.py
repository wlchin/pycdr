
# get the data file as a fixture with scope module

# there is no mock
import pytest
from pycdr import pycdr
import anndata as ad
from pathlib import Path
from pycdr import perm
from pycdr import utils

@pytest.fixture(scope='module')
def muscleobject():
    """
    return the default data file that comes with the package
    """
    adatapath = Path(__file__).parents[1] / "data/test_muscle.h5ad"
    x = ad.read(adatapath)
    return x
    
def test_CDR_muscle_adata_shape(muscleobject, tmpdir):
    pycdr.run_CDR_analysis(muscleobject, "Hours")
    file = tmpdir.join('output.h5ad')
    muscleobject.write(file.strpath)  # or use str(file)
    x = ad.read(file.strpath)
    assert x.shape[0] == 100

def test_CDR_muscle_Fs(muscleobject, tmpdir):
    pycdr.run_CDR_analysis(muscleobject, "Hours")
    file = tmpdir.join('output.h5ad')
    muscleobject.write(file.strpath)  # or use str(file)
    x = ad.read(file.strpath)
    assert x.uns["Fs"].shape == (200,22)

    
def test_CDR_muscle_loading_selection(muscleobject, tmpdir):
    pycdr.run_CDR_analysis(muscleobject, "Hours")
    file = tmpdir.join('output.h5ad')
    muscleobject.write(file.strpath)  # or use str(file)
    x = ad.read(file.strpath)
    assert x.uns["selected_loading"] == 22

def test_CDR_muscle_Fs_mat(muscleobject, tmpdir):
    pycdr.run_CDR_analysis(muscleobject, "Hours")
    file = tmpdir.join('output.h5ad')
    muscleobject.write(file.strpath)  # or use str(file)
    x = ad.read(file.strpath)
    assert round(x.uns["Fs"][25,3],4) == 0.1033

def test_CDR_muscle_Fs_shape_pval(muscleobject, tmpdir):
    pycdr.run_CDR_analysis(muscleobject, "Hours")
    file = tmpdir.join('output.h5ad')
    muscleobject.write(file.strpath)  # or use str(file)
    x = ad.read(file.strpath)
    assert x.uns["pval_mat"].shape == (100, 22)

def test_CDR_muscle_Fs_pval_values(muscleobject, tmpdir):
    pycdr.run_CDR_analysis(muscleobject, "Hours")
    file = tmpdir.join('output.h5ad')
    muscleobject.write(file.strpath)  # or use str(file)
    x = ad.read(file.strpath)
    assert x.uns["pval_mat"][15,10] == 0.4535

def test_CDR_muscle_Fs_factor_loading(muscleobject, tmpdir):
    pycdr.run_CDR_analysis(muscleobject, "Hours")
    file = tmpdir.join('output.h5ad')
    muscleobject.write(file.strpath)  # or use str(file)
    x = ad.read(file.strpath)
    assert x.uns['factor_loadings']['factor.2'][0] == 'ARHGAP33'
    assert x.uns['factor_loadings']['factor.2'][0] == 'ARHGAP33'
    assert x.uns['factor_loadings']['factor.2'][0] == 'ARHGAP33'

def test_CDR_muscle_Fs_diff(muscleobject, tmpdir):
    pycdr.run_CDR_analysis(muscleobject, "Hours")
    file = tmpdir.join('output.h5ad')
    muscleobject.write(file.strpath)  # or use str(file)
    x = ad.read(file.strpath)
    assert round(x.uns['Fs_diff'][23,6], 4) == 0.235

def test_enrichment(muscleobject, tmpdir):
    pycdr.run_CDR_analysis(muscleobject, "Hours")
    file = tmpdir.join('output.h5ad')
    muscleobject.write(file.strpath)  # or use str(file)
    x = ad.read(file.strpath)
    a, b = perm.calculate_enrichment(x, "Hours",  ['factor.9'] , 10, "gene_short_name", 0.1)
    assert a["factor.9"][0] == [9,0]
    assert b["factor.9"][1] == 0.011916939343153464

def test_df_loading(muscleobject, tmpdir):
    pycdr.run_CDR_analysis(muscleobject, "Hours")
    file = tmpdir.join('output.h5ad')
    muscleobject.write(file.strpath)  # or use str(file)
    x = ad.read(file.strpath)
    a, b = perm.calculate_enrichment(x, "Hours",  ['factor.9'] , 10, "gene_short_name", 0.1)
    c = perm.get_df_loadings(x).fdr[0]
    d = perm.get_df_loadings(x).a_max[0]
    assert c == 0.011916939343153464
    assert d == 0.14754098360655737
    
def test_filter_percent(muscleobject, tmpdir):
    pycdr.run_CDR_analysis(muscleobject, "Hours")
    file = tmpdir.join('output.h5ad')
    muscleobject.write(file.strpath)  # or use str(file)
    x = ad.read(file.strpath)
    y = utils.filter_genecounts_percent(x, 5, 4)
    assert y.shape == (100, 6)

def test_filter_counts(muscleobject, tmpdir):
    pycdr.run_CDR_analysis(muscleobject, "Hours")
    file = tmpdir.join('output.h5ad')
    muscleobject.write(file.strpath)  # or use str(file)
    x = ad.read(file.strpath)
    y = utils.filter_genecounts_numcells(x, 5, 25)
    assert y.shape == (100, 35)

def test_top_genes(muscleobject, tmpdir):
    pycdr.run_CDR_analysis(muscleobject, "Hours")
    file = tmpdir.join('output.h5ad')
    muscleobject.write(file.strpath)  # or use str(file)
    x = ad.read(file.strpath)
    y = utils.get_top_genes(x, 1)
    assert round(y["z_score"][0],4) == 2.3099
    assert round(y["pval"][2], 4) == 0.1295


