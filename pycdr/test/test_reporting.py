"""Tests for the pycdr reporting module."""

import pytest
import numpy as np
import pandas as pd
import anndata as ad
from pathlib import Path

from pycdr import pycdr as pycdr_mod
from pycdr import kruskal
from pycdr.reporting import (
    get_factor_summary,
    format_run_summary,
    format_info_summary,
    format_results_table,
    generate_html_report,
)


@pytest.fixture(scope="module")
def analyzed_adata():
    """Run CDR analysis on test data and return the AnnData."""
    adatapath = Path(__file__).parents[1] / "data/test_muscle.h5ad"
    adata = ad.read_h5ad(str(adatapath))
    pycdr_mod.run_CDR_analysis(adata, "Hours", correction="none")
    return adata


@pytest.fixture(scope="module")
def enriched_adata(analyzed_adata):
    """Run Kruskal enrichment on analyzed data."""
    adata = analyzed_adata.copy()
    kruskal.calculate_enrichment(adata, "Hours")
    return adata


# ---------------------------------------------------------------------------
# get_factor_summary
# ---------------------------------------------------------------------------

def test_factor_summary_columns(analyzed_adata):
    df = get_factor_summary(analyzed_adata)
    assert "factor" in df.columns
    assert "n_genes" in df.columns
    assert "top_genes" in df.columns
    assert "mean_zscore" in df.columns


def test_factor_summary_row_count(analyzed_adata):
    df = get_factor_summary(analyzed_adata)
    n_factors = len(analyzed_adata.uns["factor_loadings"])
    assert len(df) == n_factors


def test_factor_summary_with_enrichment(enriched_adata):
    df = get_factor_summary(enriched_adata)
    assert "enrich_fdr" in df.columns


# ---------------------------------------------------------------------------
# format_run_summary
# ---------------------------------------------------------------------------

def test_run_summary_contains_header(analyzed_adata):
    out = format_run_summary(analyzed_adata, "Hours")
    assert "CDR-g Analysis Summary" in out


def test_run_summary_contains_factors(analyzed_adata):
    out = format_run_summary(analyzed_adata, "Hours")
    # Should mention factor names and gene counts
    assert "factor." in out
    assert "Genes" in out


def test_run_summary_with_enrichment(enriched_adata):
    out = format_run_summary(enriched_adata, "Hours", enriched=True)
    assert "FDR" in out


# ---------------------------------------------------------------------------
# format_info_summary
# ---------------------------------------------------------------------------

def test_info_summary_analyzed(analyzed_adata):
    out = format_info_summary(analyzed_adata)
    assert out is not None
    assert "CDR-g Analysis Results" in out
    assert "unique" in out


def test_info_summary_no_analysis():
    adata = ad.AnnData(np.random.default_rng(0).random((5, 5)))
    out = format_info_summary(adata)
    assert out is None


# ---------------------------------------------------------------------------
# format_results_table
# ---------------------------------------------------------------------------

def test_results_table_ascii():
    df = pd.DataFrame({
        "n_genes": [10, 5],
        "genes": [["A", "B"], ["C"]],
        "top_genes": ["A, B", "C"],
    })
    out = format_results_table(df, fmt="table")
    # Should have aligned columns with separator line
    lines = out.strip().split("\n")
    assert len(lines) >= 3  # header, separator, rows
    assert "---" in lines[1]


def test_results_table_markdown():
    df = pd.DataFrame({
        "n_genes": [10, 5],
        "genes": [["A", "B"], ["C"]],
        "top_genes": ["A, B", "C"],
    })
    out = format_results_table(df, fmt="markdown")
    assert "|" in out
    lines = out.strip().split("\n")
    assert all(line.startswith("|") for line in lines)


def test_results_table_csv():
    df = pd.DataFrame({
        "n_genes": [10],
        "genes": [["A", "B", "C"]],
        "top_genes": ["A, B, C"],
    })
    out = format_results_table(df, fmt="csv")
    assert "," in out
    # genes column should be serialized
    assert "A,B,C" in out


# ---------------------------------------------------------------------------
# generate_html_report
# ---------------------------------------------------------------------------

def test_html_report_basic(analyzed_adata, tmp_path):
    out = tmp_path / "report.html"
    generate_html_report(analyzed_adata, str(out), "Hours")
    content = out.read_text()
    assert "<html" in content
    assert "factor." in content


def test_html_report_with_params(analyzed_adata, tmp_path):
    adata = analyzed_adata.copy()
    adata.uns["cdr_params"] = {
        "phenotype": "Hours",
        "capvar": 0.95,
        "pernum": 2000,
        "thres": 0.05,
        "filter_method": "numcells",
        "count_threshold": 1,
        "min_cells": 10,
        "enrich_method": "kruskal",
        "subset": [],
    }
    out = tmp_path / "report_params.html"
    generate_html_report(adata, str(out), "Hours")
    content = out.read_text()
    assert "Run Parameters" in content
    assert "capvar" in content
    assert "numcells" in content
    assert "kruskal" in content


def test_html_report_no_params(analyzed_adata, tmp_path):
    adata = analyzed_adata.copy()
    adata.uns.pop("cdr_params", None)
    out = tmp_path / "report_no_params.html"
    generate_html_report(adata, str(out), "Hours")
    content = out.read_text()
    assert "Run Parameters" not in content


def test_html_report_without_matplotlib(analyzed_adata, tmp_path, monkeypatch):
    """When matplotlib is not available, report should include install note."""
    import pycdr.reporting as rep_mod
    # Monkey-patch _try_embed_figure to simulate missing matplotlib
    original = rep_mod._try_embed_figure
    def _mock_no_mpl(adata):
        return (
            '<p class="note">Install matplotlib for embedded plots: '
            '<code>pip install cdr-py[plot]</code></p>'
        )
    monkeypatch.setattr(rep_mod, "_try_embed_figure", _mock_no_mpl)
    out = tmp_path / "report_no_mpl.html"
    generate_html_report(analyzed_adata, str(out), "Hours")
    content = out.read_text()
    assert "Install matplotlib" in content
