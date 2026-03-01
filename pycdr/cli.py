"""pycdr command-line interface."""

import logging
import sys
from pathlib import Path

import click

logger = logging.getLogger("pycdr")


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------

def validate_phenotype(adata, col):
    """Check that *col* exists in adata.obs; raise click.BadParameter otherwise."""
    if col not in adata.obs.columns:
        available = ", ".join(adata.obs.columns.tolist())
        raise click.BadParameter(
            f"Column '{col}' not found in adata.obs. Available: {available}",
            param_hint="'-p'",
        )


def validate_analyzed(adata):
    """Check that the dataset contains CDR-g analysis results."""
    if "factor_loadings" not in adata.uns:
        raise click.ClickException(
            "No CDR-g analysis found in this dataset. Run 'pycdr analyze' first."
        )


def validate_genecol(adata, col):
    """Check that *col* exists in adata.var."""
    if col not in adata.var.columns:
        available = ", ".join(adata.var.columns.tolist())
        raise click.BadParameter(
            f"Column '{col}' not found in adata.var. Available: {available}",
            param_hint="'--genecol'",
        )


def default_output(input_path, suffix="_cdr"):
    """Generate an output path from the input stem."""
    p = Path(input_path)
    return str(p.with_name(f"{p.stem}{suffix}{p.suffix}"))


def _setup_logging(verbose, quiet):
    """Configure the pycdr root logger based on -v/-q flags."""
    if quiet:
        level = logging.ERROR
    elif verbose >= 2:
        level = logging.DEBUG
    elif verbose == 1:
        level = logging.INFO
    else:
        level = logging.WARNING

    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    root = logging.getLogger("pycdr")
    root.handlers.clear()
    root.addHandler(handler)
    root.setLevel(level)
    # Also set level on library sub-loggers so their messages propagate
    for name in ("pycdr.pycdr", "pycdr.perm", "pycdr.utils",
                 "pycdr.feature_selection", "pycdr.experimental"):
        logging.getLogger(name).setLevel(level)


def _read_h5ad(path):
    """Read an .h5ad file, raising a nice error on failure."""
    import anndata as ad

    p = Path(path)
    if not p.exists():
        raise click.BadParameter(f"File not found: {path}", param_hint="'INPUT'")
    if p.suffix != ".h5ad":
        raise click.BadParameter(
            f"Expected .h5ad file, got '{p.suffix}'", param_hint="'INPUT'"
        )
    return ad.read_h5ad(str(p))


# ---------------------------------------------------------------------------
# CLI group
# ---------------------------------------------------------------------------

@click.group()
@click.version_option(package_name="cdr-py")
def cli():
    """pycdr — CDR-g analysis for single-cell RNA-seq data."""


# ---------------------------------------------------------------------------
# info
# ---------------------------------------------------------------------------

@cli.command()
@click.argument("input", type=click.Path(exists=True))
@click.option("-p", "--phenotype", default=None, help="Show value counts for this obs column.")
def info(input, phenotype):
    """Display dataset metadata."""
    adata = _read_h5ad(input)

    click.echo(f"Shape: {adata.shape[0]} cells x {adata.shape[1]} genes")

    if "factor_loadings" in adata.uns:
        n = len(adata.uns["factor_loadings"])
        click.echo(f"CDR-g analysis: found ({n} factors)")
    else:
        click.echo("CDR-g analysis: not found")

    if phenotype is not None:
        validate_phenotype(adata, phenotype)
        counts = adata.obs[phenotype].value_counts()
        parts = ", ".join(f"{k}={v}" for k, v in counts.items())
        click.echo(f"Phenotype '{phenotype}': {parts}")


# ---------------------------------------------------------------------------
# analyze
# ---------------------------------------------------------------------------

@cli.command()
@click.argument("input", type=click.Path(exists=True))
@click.option("-p", "--phenotype", required=True, help="Condition column in adata.obs.")
@click.option("-o", "--output", default=None, help="Output .h5ad path.")
@click.option("--capvar", default=0.95, show_default=True, help="Variance threshold.")
@click.option("--pernum", default=2000, show_default=True, help="Permutations for importance scores.")
@click.option("--thres", default=0.05, show_default=True, help="P-value threshold for gene selection.")
@click.option("-v", "--verbose", count=True, help="Increase verbosity (-v INFO, -vv DEBUG).")
@click.option("-q", "--quiet", is_flag=True, help="Errors only.")
def analyze(input, phenotype, output, capvar, pernum, thres, verbose, quiet):
    """Run CDR-g SVD/varimax analysis."""
    _setup_logging(verbose, quiet)
    from .pycdr import run_CDR_analysis

    adata = _read_h5ad(input)
    validate_phenotype(adata, phenotype)

    logger.info("Running CDR-g analysis on %s", input)
    run_CDR_analysis(adata, phenotype, capvar=capvar, pernum=pernum, thres=thres)

    out = output or default_output(input)
    adata.write(out)
    click.echo(f"Wrote {out}")


# ---------------------------------------------------------------------------
# filter
# ---------------------------------------------------------------------------

@cli.command("filter")
@click.argument("input", type=click.Path(exists=True))
@click.option("-m", "--method", "filter_method", required=True,
              type=click.Choice(["percent", "numcells"]), help="Filter method.")
@click.option("--cell-fraction", default=0.05, show_default=True,
              help="(percent) Fraction of cells.")
@click.option("--median-count", default=1.0, show_default=True,
              help="(percent) Median count threshold.")
@click.option("--count-threshold", default=1, show_default=True,
              help="(numcells) Count cutoff.")
@click.option("--min-cells", default=10, show_default=True,
              help="(numcells) Min expressed cells.")
@click.option("-o", "--output", default=None, help="Output .h5ad path.")
@click.option("-v", "--verbose", count=True)
@click.option("-q", "--quiet", is_flag=True)
def filter_cmd(input, filter_method, cell_fraction, median_count,
               count_threshold, min_cells, output, verbose, quiet):
    """Filter genes from an .h5ad file."""
    _setup_logging(verbose, quiet)
    from .utils import filter_genecounts_percent, filter_genecounts_numcells

    adata = _read_h5ad(input)

    if filter_method == "percent":
        logger.info("Filtering with percent method (fraction=%.3f, median_count=%.1f)",
                     cell_fraction, median_count)
        adata = filter_genecounts_percent(adata, cell_fraction, median_count)
    else:
        logger.info("Filtering with numcells method (count_threshold=%d, min_cells=%d)",
                     count_threshold, min_cells)
        adata = filter_genecounts_numcells(adata, count_threshold, min_cells)

    click.echo(f"After filtering: {adata.shape[0]} cells x {adata.shape[1]} genes")

    out = output or default_output(input, suffix="_filtered")
    adata.write(out)
    click.echo(f"Wrote {out}")


# ---------------------------------------------------------------------------
# enrich
# ---------------------------------------------------------------------------

@cli.command()
@click.argument("input", type=click.Path(exists=True))
@click.option("-p", "--phenotype", required=True, help="Condition column in adata.obs.")
@click.option("-m", "--method", "enrich_method", default="perm", show_default=True,
              type=click.Choice(["perm", "kruskal"]), help="Enrichment method.")
@click.option("--genecol", default=None,
              help="Gene name column in adata.var (required for perm method).")
@click.option("--nperm", default=100, show_default=True, help="Permutations for ssGSEA.")
@click.option("--enrich-thresh", default=0.05, show_default=True,
              help="Active gene set threshold.")
@click.option("--seed", default=42, show_default=True, help="Random seed.")
@click.option("--factors", default=None,
              help="Comma-separated factor names to test (default: all).")
@click.option("-o", "--output", default=None, help="Output .h5ad path.")
@click.option("-c", "--csv", default=None, help="Export results CSV.")
@click.option("-v", "--verbose", count=True)
@click.option("-q", "--quiet", is_flag=True)
def enrich(input, phenotype, enrich_method, genecol, nperm, enrich_thresh,
           seed, factors, output, csv, verbose, quiet):
    """Run enrichment on previously analyzed data."""
    _setup_logging(verbose, quiet)

    adata = _read_h5ad(input)
    validate_phenotype(adata, phenotype)
    validate_analyzed(adata)

    factor_list = list(adata.uns["factor_loadings"].keys())
    if factors is not None:
        factor_list = [f.strip() for f in factors.split(",")]

    if enrich_method == "perm":
        from .perm import calculate_enrichment

        if genecol is None:
            raise click.BadParameter(
                "--genecol is required for the 'perm' enrichment method.",
                param_hint="'--genecol'",
            )
        validate_genecol(adata, genecol)

        logger.info("Running perm enrichment (%d factors, nperm=%d)", len(factor_list), nperm)
        calculate_enrichment(adata, phenotype, factor_list, nperm, genecol, enrich_thresh, seed=seed)
    else:
        from .experimental import calculate_enrichment as calc_kruskal

        logger.info("Running Kruskal-Wallis enrichment (%d factors)", len(factor_list))
        calc_kruskal(adata, phenotype)

    out = output or default_output(input, suffix="_enriched")
    adata.write(out)
    click.echo(f"Wrote {out}")

    if csv:
        from .utils import output_results
        df = output_results(adata)
        if df is not None:
            df.to_csv(csv, index=False)
            click.echo(f"Wrote {csv}")


# ---------------------------------------------------------------------------
# results
# ---------------------------------------------------------------------------

@cli.command()
@click.argument("input", type=click.Path(exists=True))
@click.option("-o", "--output", default=None, help="Output file (CSV/TSV). Stdout if omitted.")
@click.option("-f", "--format", "fmt", default="csv", show_default=True,
              type=click.Choice(["csv", "tsv"]), help="Output format.")
@click.option("--top-genes", default=None, type=int, help="Show top N genes per factor.")
@click.option("--factor", default=None, type=int, help="Show results for a single factor index.")
@click.option("-v", "--verbose", count=True)
@click.option("-q", "--quiet", is_flag=True)
def results(input, output, fmt, top_genes, factor, verbose, quiet):
    """Export results table to CSV/TSV or stdout."""
    _setup_logging(verbose, quiet)
    from .utils import output_results, get_top_genes

    adata = _read_h5ad(input)
    validate_analyzed(adata)

    if factor is not None:
        df = get_top_genes(adata, factor)
        if top_genes is not None:
            df = df.head(top_genes)
    else:
        df = output_results(adata)
        if df is None:
            raise click.ClickException("No results to export.")

    sep = "\t" if fmt == "tsv" else ","
    if output:
        df.to_csv(output, sep=sep, index=(factor is not None))
        click.echo(f"Wrote {output}")
    else:
        click.echo(df.to_csv(sep=sep, index=(factor is not None)), nl=False)


# ---------------------------------------------------------------------------
# run  (full pipeline)
# ---------------------------------------------------------------------------

@cli.command()
@click.argument("input", type=click.Path(exists=True))
@click.option("-p", "--phenotype", required=True, help="Condition column in adata.obs.")
@click.option("-o", "--output", default=None, help="Output .h5ad path.")
@click.option("-c", "--csv", default=None, help="Export results CSV.")
@click.option("--capvar", default=0.95, show_default=True, help="Variance threshold.")
@click.option("--pernum", default=2000, show_default=True, help="Permutations for importance.")
@click.option("--thres", default=0.05, show_default=True, help="P-value threshold.")
@click.option("--filter-method", default="none", show_default=True,
              type=click.Choice(["none", "percent", "numcells"]), help="Gene filter method.")
@click.option("--cell-fraction", default=0.05, show_default=True)
@click.option("--median-count", default=1.0, show_default=True)
@click.option("--count-threshold", default=1, show_default=True)
@click.option("--min-cells", default=10, show_default=True)
@click.option("--enrich/--no-enrich", default=False, show_default=True,
              help="Run enrichment after analysis.")
@click.option("--enrich-method", default="perm", show_default=True,
              type=click.Choice(["perm", "kruskal"]))
@click.option("--genecol", default=None, help="Gene name column in adata.var.")
@click.option("--nperm", default=100, show_default=True)
@click.option("--enrich-thresh", default=0.05, show_default=True)
@click.option("--seed", default=42, show_default=True, help="Random seed.")
@click.option("-v", "--verbose", count=True)
@click.option("-q", "--quiet", is_flag=True)
def run(input, phenotype, output, csv, capvar, pernum, thres,
        filter_method, cell_fraction, median_count, count_threshold, min_cells,
        enrich, enrich_method, genecol, nperm, enrich_thresh, seed, verbose, quiet):
    """Full CDR-g pipeline: filter, analyze, enrich, export."""
    _setup_logging(verbose, quiet)
    from .pycdr import run_CDR_analysis
    from .utils import (filter_genecounts_percent, filter_genecounts_numcells,
                        output_results)

    adata = _read_h5ad(input)
    validate_phenotype(adata, phenotype)

    # --- filter ---
    if filter_method == "percent":
        logger.info("Filtering genes (percent method)")
        adata = filter_genecounts_percent(adata, cell_fraction, median_count)
        click.echo(f"After filtering: {adata.shape[0]} cells x {adata.shape[1]} genes")
    elif filter_method == "numcells":
        logger.info("Filtering genes (numcells method)")
        adata = filter_genecounts_numcells(adata, count_threshold, min_cells)
        click.echo(f"After filtering: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # --- analyze ---
    logger.info("Running CDR-g analysis")
    run_CDR_analysis(adata, phenotype, capvar=capvar, pernum=pernum, thres=thres)
    n_factors = len(adata.uns["factor_loadings"])
    click.echo(f"CDR-g analysis complete: {n_factors} factors")

    # --- enrich ---
    if enrich:
        factor_list = list(adata.uns["factor_loadings"].keys())

        if enrich_method == "perm":
            from .perm import calculate_enrichment

            if genecol is None:
                raise click.BadParameter(
                    "--genecol is required when --enrich is used with perm method.",
                    param_hint="'--genecol'",
                )
            validate_genecol(adata, genecol)
            logger.info("Running perm enrichment")
            calculate_enrichment(adata, phenotype, factor_list, nperm, genecol, enrich_thresh, seed=seed)
        else:
            from .experimental import calculate_enrichment as calc_kruskal

            logger.info("Running Kruskal-Wallis enrichment")
            calc_kruskal(adata, phenotype)

        click.echo("Enrichment complete")

    # --- output ---
    out = output or default_output(input)
    adata.write(out)
    click.echo(f"Wrote {out}")

    if csv:
        df = output_results(adata)
        if df is not None:
            df.to_csv(csv, index=False)
            click.echo(f"Wrote {csv}")
