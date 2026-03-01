"""Tests for the pycdr CLI."""

import pytest
import numpy as np
import anndata as ad
from pathlib import Path
from click.testing import CliRunner

from pycdr.cli import cli
from pycdr import pycdr as pycdr_mod


@pytest.fixture(scope="module")
def runner():
    return CliRunner()


@pytest.fixture(scope="module")
def muscle_path():
    return str(Path(__file__).parents[1] / "data/test_muscle.h5ad")


@pytest.fixture(scope="module")
def analyzed_path(muscle_path, tmp_path_factory):
    """Run CDR analysis and write to a temp file."""
    adata = ad.read_h5ad(muscle_path)
    pycdr_mod.run_CDR_analysis(adata, "Hours")
    out = tmp_path_factory.mktemp("cli") / "analyzed.h5ad"
    adata.write(str(out))
    return str(out)


@pytest.fixture(scope="module")
def output_path():
    return str(Path(__file__).parents[1] / "data/test_fl_e_s.h5ad")


# ---------------------------------------------------------------------------
# Top-level
# ---------------------------------------------------------------------------

def test_help(runner):
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "pycdr" in result.output


def test_version(runner):
    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0


# ---------------------------------------------------------------------------
# info
# ---------------------------------------------------------------------------

def test_info_basic(runner, muscle_path):
    result = runner.invoke(cli, ["info", muscle_path])
    assert result.exit_code == 0
    assert "cells" in result.output
    assert "genes" in result.output


def test_info_with_phenotype(runner, muscle_path):
    result = runner.invoke(cli, ["info", muscle_path, "-p", "Hours"])
    assert result.exit_code == 0
    assert "Phenotype" in result.output


def test_info_bad_phenotype(runner, muscle_path):
    result = runner.invoke(cli, ["info", muscle_path, "-p", "NONEXISTENT"])
    assert result.exit_code != 0


def test_info_analyzed(runner, analyzed_path):
    result = runner.invoke(cli, ["info", analyzed_path])
    assert result.exit_code == 0
    assert "CDR-g Analysis Results" in result.output
    assert "unique" in result.output


# ---------------------------------------------------------------------------
# analyze
# ---------------------------------------------------------------------------

def test_analyze(runner, muscle_path, tmp_path_factory):
    out = str(tmp_path_factory.mktemp("cli_analyze") / "out.h5ad")
    result = runner.invoke(cli, [
        "analyze", muscle_path, "-p", "Hours", "-o", out, "-v",
    ])
    assert result.exit_code == 0
    assert "Wrote" in result.output
    adata = ad.read_h5ad(out)
    assert "factor_loadings" in adata.uns


def test_analyze_bad_phenotype(runner, muscle_path):
    result = runner.invoke(cli, ["analyze", muscle_path, "-p", "BAD"])
    assert result.exit_code != 0


# ---------------------------------------------------------------------------
# filter
# ---------------------------------------------------------------------------

def test_filter_percent(runner, muscle_path, tmp_path_factory):
    out = str(tmp_path_factory.mktemp("cli_filter") / "filt.h5ad")
    result = runner.invoke(cli, [
        "filter", muscle_path, "-m", "percent",
        "--cell-fraction", "5", "--median-count", "4", "-o", out,
    ])
    assert result.exit_code == 0
    assert "After filtering" in result.output
    adata = ad.read_h5ad(out)
    assert adata.shape[1] < 100  # fewer genes after filtering


def test_filter_numcells(runner, muscle_path, tmp_path_factory):
    out = str(tmp_path_factory.mktemp("cli_filter2") / "filt.h5ad")
    result = runner.invoke(cli, [
        "filter", muscle_path, "-m", "numcells",
        "--count-threshold", "5", "--min-cells", "25", "-o", out,
    ])
    assert result.exit_code == 0
    assert "After filtering" in result.output


# ---------------------------------------------------------------------------
# enrich (perm method)
# ---------------------------------------------------------------------------

def test_enrich_perm(runner, analyzed_path, tmp_path_factory):
    out = str(tmp_path_factory.mktemp("cli_enrich") / "enriched.h5ad")
    result = runner.invoke(cli, [
        "enrich", analyzed_path, "-p", "Hours", "-m", "perm",
        "--genecol", "gene_short_name", "--nperm", "10",
        "--factors", "factor.9", "-o", out,
    ])
    assert result.exit_code == 0
    assert "Wrote" in result.output


def test_enrich_perm_missing_genecol(runner, analyzed_path):
    result = runner.invoke(cli, [
        "enrich", analyzed_path, "-p", "Hours", "-m", "perm",
    ])
    assert result.exit_code != 0


def test_enrich_kruskal(runner, analyzed_path, tmp_path_factory):
    out = str(tmp_path_factory.mktemp("cli_enrich_kw") / "enriched.h5ad")
    result = runner.invoke(cli, [
        "enrich", analyzed_path, "-p", "Hours", "-m", "kruskal", "-o", out,
    ])
    assert result.exit_code == 0


# ---------------------------------------------------------------------------
# results
# ---------------------------------------------------------------------------

def test_results_stdout(runner, output_path):
    result = runner.invoke(cli, ["results", output_path])
    assert result.exit_code == 0
    assert len(result.output) > 0


def test_results_to_file(runner, output_path, tmp_path_factory):
    out = str(tmp_path_factory.mktemp("cli_results") / "res.csv")
    result = runner.invoke(cli, ["results", output_path, "-o", out])
    assert result.exit_code == 0
    assert "Wrote" in result.output
    assert Path(out).exists()


def test_results_tsv(runner, output_path):
    result = runner.invoke(cli, ["results", output_path, "-f", "tsv"])
    assert result.exit_code == 0
    assert "\t" in result.output


def test_results_top_genes(runner, analyzed_path):
    result = runner.invoke(cli, [
        "results", analyzed_path, "--top-genes", "5", "--factor", "1",
    ])
    assert result.exit_code == 0
    lines = result.output.strip().split("\n")
    assert len(lines) == 6  # header + 5 rows


def test_results_no_analysis(runner, muscle_path):
    result = runner.invoke(cli, ["results", muscle_path])
    assert result.exit_code != 0


def test_results_table_format(runner, analyzed_path):
    result = runner.invoke(cli, ["results", analyzed_path, "-f", "table"])
    assert result.exit_code == 0
    assert "---" in result.output  # separator line


def test_results_markdown_format(runner, analyzed_path):
    result = runner.invoke(cli, ["results", analyzed_path, "-f", "markdown"])
    assert result.exit_code == 0
    assert "|" in result.output


# ---------------------------------------------------------------------------
# plot
# ---------------------------------------------------------------------------

def test_plot_command(runner, analyzed_path, tmp_path_factory):
    pytest.importorskip("matplotlib")
    out = str(tmp_path_factory.mktemp("cli_plot") / "summary.png")
    result = runner.invoke(cli, ["plot", analyzed_path, "-o", out])
    assert result.exit_code == 0
    assert "Wrote" in result.output
    assert Path(out).exists()


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

def test_report_command(runner, analyzed_path, tmp_path_factory):
    out = str(tmp_path_factory.mktemp("cli_report") / "report.html")
    result = runner.invoke(cli, [
        "report", analyzed_path, "-p", "Hours", "-o", out,
    ])
    assert result.exit_code == 0
    assert "Wrote" in result.output
    content = Path(out).read_text()
    assert "<html" in content


# ---------------------------------------------------------------------------
# run (full pipeline)
# ---------------------------------------------------------------------------

def test_run_minimal(runner, muscle_path, tmp_path_factory):
    out = str(tmp_path_factory.mktemp("cli_run") / "out.h5ad")
    result = runner.invoke(cli, [
        "run", muscle_path, "-p", "Hours", "-o", out,
    ])
    assert result.exit_code == 0
    assert "CDR-g analysis complete" in result.output
    assert "CDR-g Analysis Summary" in result.output
    assert "Wrote" in result.output
    adata = ad.read_h5ad(out)
    assert "factor_loadings" in adata.uns


def test_run_with_filter(runner, muscle_path, tmp_path_factory):
    out = str(tmp_path_factory.mktemp("cli_run2") / "out.h5ad")
    result = runner.invoke(cli, [
        "run", muscle_path, "-p", "Hours", "-o", out,
        "--filter-method", "numcells", "--min-cells", "25",
    ])
    assert result.exit_code == 0
    assert "After filtering" in result.output


def test_run_with_enrich(runner, muscle_path, tmp_path_factory):
    outdir = tmp_path_factory.mktemp("cli_run3")
    out = str(outdir / "out.h5ad")
    csv = str(outdir / "out.csv")
    result = runner.invoke(cli, [
        "run", muscle_path, "-p", "Hours", "-o", out, "-c", csv,
        "--enrich", "--enrich-method", "perm",
        "--genecol", "gene_short_name", "--nperm", "10",
    ])
    assert result.exit_code == 0
    assert "Enrichment complete" in result.output
    assert Path(csv).exists()


def test_run_with_kruskal(runner, muscle_path, tmp_path_factory):
    out = str(tmp_path_factory.mktemp("cli_run4") / "out.h5ad")
    result = runner.invoke(cli, [
        "run", muscle_path, "-p", "Hours", "-o", out,
        "--enrich", "--enrich-method", "kruskal",
    ])
    assert result.exit_code == 0
    assert "Enrichment complete" in result.output
    assert "FDR" in result.output


def test_run_bad_phenotype(runner, muscle_path):
    result = runner.invoke(cli, ["run", muscle_path, "-p", "NONEXISTENT"])
    assert result.exit_code != 0


def test_run_enrich_missing_genecol(runner, muscle_path):
    result = runner.invoke(cli, [
        "run", muscle_path, "-p", "Hours", "--enrich",
    ])
    assert result.exit_code != 0
