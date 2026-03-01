"""pycdr command-line interface."""

import logging
import sys
import warnings
from pathlib import Path

import click

logger = logging.getLogger("pycdr")

# Module-level quiet flag for _echo helper
_quiet_mode = False


def _echo(message, **kwargs):
    """click.echo wrapper that respects quiet mode."""
    if not _quiet_mode:
        click.echo(message, **kwargs)


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


def subset_cells(adata, subset_specs):
    """Subset *adata* by one or more ``column=value[,value2]`` specs.

    Multiple specs are combined with AND logic.  Within a single spec,
    comma-separated values are combined with OR logic.
    """
    if not subset_specs:
        return adata

    for spec in subset_specs:
        if "=" not in spec:
            raise click.BadParameter(
                f"Invalid subset format '{spec}'. Expected COLUMN=VALUE[,VALUE2,...]",
                param_hint="'--subset'",
            )
        col, raw_values = spec.split("=", 1)
        if col not in adata.obs.columns:
            available = ", ".join(adata.obs.columns.tolist())
            raise click.BadParameter(
                f"Column '{col}' not found in adata.obs. Available: {available}",
                param_hint="'--subset'",
            )
        values = [v.strip() for v in raw_values.split(",")]
        mask = adata.obs[col].astype(str).isin(values)
        if not mask.any():
            raise click.ClickException(
                f"No cells match --subset '{spec}'. "
                f"Unique values in '{col}': {adata.obs[col].unique().tolist()}"
            )
        adata = adata[mask].copy()

    logger.info("After subsetting: %d cells x %d genes", adata.shape[0], adata.shape[1])
    click.echo(f"Subset: {adata.shape[0]} cells x {adata.shape[1]} genes")
    return adata


def _validate_phenotype_after_subset(adata, phenotype):
    """Warn if any phenotype group has 0 cells after subsetting."""
    counts = adata.obs[phenotype].value_counts()
    empty = counts[counts == 0]
    if len(empty) > 0:
        raise click.ClickException(
            f"After subsetting, phenotype group(s) {empty.index.tolist()} "
            f"have 0 cells. CDR-g requires cells in every condition."
        )


def default_output(input_path, suffix="_cdr"):
    """Generate an output path from the input stem."""
    p = Path(input_path)
    return str(p.with_name(f"{p.stem}{suffix}{p.suffix}"))


def _setup_logging(verbose, quiet):
    """Configure the pycdr root logger based on -v/-q flags.

    Returns:
        bool: True if quiet mode is active (level >= ERROR).
    """
    global _quiet_mode

    if quiet:
        level = logging.ERROR
    elif verbose >= 2:
        level = logging.DEBUG
    elif verbose == 1:
        level = logging.INFO
    else:
        level = logging.WARNING

    if level <= logging.DEBUG:
        fmt = "%(asctime)s %(levelname)s: %(message)s"
    else:
        fmt = "%(levelname)s: %(message)s"

    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter(fmt))
    root = logging.getLogger("pycdr")
    root.handlers.clear()
    root.addHandler(handler)
    root.setLevel(level)
    # Also set level on library sub-loggers so their messages propagate
    for name in ("pycdr.pycdr", "pycdr.perm", "pycdr.utils",
                 "pycdr.feature_selection", "pycdr.kruskal"):
        logging.getLogger(name).setLevel(level)

    # Suppress noisy Python warnings at all levels except DEBUG
    if level > logging.DEBUG:
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        warnings.filterwarnings("ignore", category=UserWarning)
        warnings.filterwarnings("ignore", category=FutureWarning)

    _quiet_mode = (level >= logging.ERROR)
    return _quiet_mode


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


def _serialize_genes_column(df):
    """Serialize the ``genes`` list column to comma-joined strings for CSV export."""
    out = df.copy()
    if "genes" in out.columns:
        out["genes"] = out["genes"].apply(
            lambda x: ",".join(x) if isinstance(x, list) else str(x)
        )
    return out


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
@click.option("-s", "--subset", multiple=True,
              help="Subset cells: COLUMN=VALUE[,VALUE2]. Repeatable.")
@click.option("-v", "--verbose", count=True, help="Increase verbosity (-v INFO, -vv DEBUG).")
@click.option("-q", "--quiet", is_flag=True, help="Errors only.")
def info(input, phenotype, subset, verbose, quiet):
    """Display dataset metadata."""
    _setup_logging(verbose, quiet)
    from .reporting import format_info_summary

    adata = _read_h5ad(input)
    adata = subset_cells(adata, subset)

    click.echo(f"Shape: {adata.shape[0]} cells x {adata.shape[1]} genes")

    if phenotype is not None:
        validate_phenotype(adata, phenotype)
        counts = adata.obs[phenotype].value_counts()
        parts = ", ".join(f"{k}={v}" for k, v in counts.items())
        click.echo(f"Phenotype '{phenotype}': {parts}")

    summary = format_info_summary(adata)
    if summary is not None:
        click.echo(summary)
    else:
        click.echo("CDR-g analysis: not found")


# ---------------------------------------------------------------------------
# analyze
# ---------------------------------------------------------------------------

@cli.command()
@click.argument("input", type=click.Path(exists=True))
@click.option("-p", "--phenotype", required=True, help="Condition column in adata.obs.")
@click.option("-o", "--output", default=None, help="Output .h5ad path.")
@click.option("--capvar", default=0.95, show_default=True, help="Variance threshold.")
@click.option("--nperm", default=2000, show_default=True, help="Permutations for importance scores.")
@click.option("--pernum", default=None, type=int, hidden=True, help="Deprecated alias for --nperm.")
@click.option("--thres", default=0.05, show_default=True, help="P-value threshold for gene selection.")
@click.option("--seed", default=42, show_default=True, help="Random seed for reproducibility.")
@click.option("--correction", default="fdr_bh", show_default=True,
              type=click.Choice(["fdr_bh", "none"]), help="Multiple testing correction method.")
@click.option("-s", "--subset", multiple=True,
              help="Subset cells: COLUMN=VALUE[,VALUE2]. Repeatable.")
@click.option("-v", "--verbose", count=True, help="Increase verbosity (-v INFO, -vv DEBUG).")
@click.option("-q", "--quiet", is_flag=True, help="Errors only.")
def analyze(input, phenotype, output, capvar, nperm, pernum, thres, seed, correction, subset, verbose, quiet):
    """Run CDR-g SVD/varimax analysis."""
    _setup_logging(verbose, quiet)
    from .pycdr import run_CDR_analysis
    from .reporting import format_run_summary
    import warnings

    if pernum is not None:
        warnings.warn("--pernum is deprecated; use --nperm instead.", DeprecationWarning, stacklevel=2)
        nperm = pernum

    adata = _read_h5ad(input)
    adata = subset_cells(adata, subset)
    validate_phenotype(adata, phenotype)
    _validate_phenotype_after_subset(adata, phenotype)

    logger.info("Running CDR-g analysis on %s", input)
    run_CDR_analysis(adata, phenotype, capvar=capvar, nperm=nperm, thres=thres,
                     seed=seed, correction=correction, quiet=quiet)

    adata.uns["cdr_params"] = {
        "phenotype": phenotype,
        "capvar": capvar,
        "nperm": nperm,
        "thres": thres,
        "seed": seed,
        "correction": correction,
        "subset": list(subset) if subset else [],
    }

    n_factors = len(adata.uns["factor_loadings"])
    _echo(f"CDR-g analysis complete: {n_factors} factors")
    _echo(format_run_summary(adata, phenotype))

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
@click.option("-s", "--subset", multiple=True,
              help="Subset cells: COLUMN=VALUE[,VALUE2]. Repeatable.")
@click.option("-o", "--output", default=None, help="Output .h5ad path.")
@click.option("-v", "--verbose", count=True)
@click.option("-q", "--quiet", is_flag=True)
def filter_cmd(input, filter_method, cell_fraction, median_count,
               count_threshold, min_cells, subset, output, verbose, quiet):
    """Filter genes from an .h5ad file."""
    _setup_logging(verbose, quiet)
    from .utils import filter_genecounts_percent, filter_genecounts_numcells

    adata = _read_h5ad(input)
    adata = subset_cells(adata, subset)

    if filter_method == "percent":
        logger.info("Filtering with percent method (fraction=%.3f, median_count=%.1f)",
                     cell_fraction, median_count)
        adata = filter_genecounts_percent(adata, cell_fraction, median_count)
    else:
        logger.info("Filtering with numcells method (count_threshold=%d, min_cells=%d)",
                     count_threshold, min_cells)
        adata = filter_genecounts_numcells(adata, count_threshold, min_cells)

    _echo(f"After filtering: {adata.shape[0]} cells x {adata.shape[1]} genes")

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
@click.option("-s", "--subset", multiple=True,
              help="Subset cells: COLUMN=VALUE[,VALUE2]. Repeatable.")
@click.option("-o", "--output", default=None, help="Output .h5ad path.")
@click.option("-c", "--csv", default=None, help="Export results CSV.")
@click.option("-v", "--verbose", count=True)
@click.option("-q", "--quiet", is_flag=True)
def enrich(input, phenotype, enrich_method, genecol, nperm, enrich_thresh,
           seed, factors, subset, output, csv, verbose, quiet):
    """Run enrichment on previously analyzed data."""
    _setup_logging(verbose, quiet)

    adata = _read_h5ad(input)
    adata = subset_cells(adata, subset)
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
        calculate_enrichment(adata, phenotype, factor_list, nperm, genecol, enrich_thresh, seed=seed, quiet=quiet)
    else:
        from .kruskal import calculate_enrichment as calc_kruskal

        logger.info("Running Kruskal-Wallis enrichment (%d factors)", len(factor_list))
        calc_kruskal(adata, phenotype, quiet=quiet)

    out = output or default_output(input, suffix="_enriched")
    adata.write(out)
    click.echo(f"Wrote {out}")

    if csv:
        from .utils import output_results
        df = output_results(adata)
        if df is not None:
            _serialize_genes_column(df).to_csv(csv, index=False)
            click.echo(f"Wrote {csv}")


# ---------------------------------------------------------------------------
# results
# ---------------------------------------------------------------------------

@cli.command()
@click.argument("input", type=click.Path(exists=True))
@click.option("-o", "--output", default=None, help="Output file (CSV/TSV). Stdout if omitted.")
@click.option("-f", "--format", "fmt", default="csv", show_default=True,
              type=click.Choice(["csv", "tsv", "table", "markdown"]), help="Output format.")
@click.option("--top-genes", default=None, type=int, help="Show top N genes per factor.")
@click.option("--factor", default=None, type=int, help="Show results for a single factor index.")
@click.option("-v", "--verbose", count=True)
@click.option("-q", "--quiet", is_flag=True)
def results(input, output, fmt, top_genes, factor, verbose, quiet):
    """Export results table to CSV/TSV/table/markdown or stdout."""
    _setup_logging(verbose, quiet)
    from .utils import output_results, get_top_genes
    from .reporting import format_results_table

    adata = _read_h5ad(input)
    validate_analyzed(adata)

    if factor is not None:
        df = get_top_genes(adata, factor)
        if top_genes is not None:
            df = df.head(top_genes)
        # Single-factor view always uses csv/tsv (no table/markdown)
        sep = "\t" if fmt == "tsv" else ","
        if output:
            df.to_csv(output, sep=sep, index=True)
            click.echo(f"Wrote {output}")
        else:
            click.echo(df.to_csv(sep=sep, index=True), nl=False)
        return

    df = output_results(adata)
    if df is None:
        raise click.ClickException("No results to export.")

    if fmt in ("table", "markdown"):
        click.echo(format_results_table(df, fmt=fmt), nl=False)
    else:
        out_df = _serialize_genes_column(df)
        sep = "\t" if fmt == "tsv" else ","
        if output:
            out_df.to_csv(output, sep=sep, index=False)
            click.echo(f"Wrote {output}")
        else:
            click.echo(out_df.to_csv(sep=sep, index=False), nl=False)


# ---------------------------------------------------------------------------
# plot
# ---------------------------------------------------------------------------

@cli.command()
@click.argument("input", type=click.Path(exists=True))
@click.option("-o", "--output", default=None, help="Output image path (default: cdr_summary.png).")
@click.option("--dpi", default=150, show_default=True, help="Figure DPI.")
@click.option("--max-factors", default=30, show_default=True,
              help="Max factors to display. Use 0 for all.")
@click.option("-v", "--verbose", count=True, help="Increase verbosity (-v INFO, -vv DEBUG).")
@click.option("-q", "--quiet", is_flag=True, help="Errors only.")
def plot(input, output, dpi, max_factors, verbose, quiet):
    """Generate a summary figure from CDR-g results."""
    _setup_logging(verbose, quiet)
    try:
        from .plotting import plot_summary
    except ImportError:
        raise click.ClickException(
            "matplotlib is required for plotting. Install with: pip install cdr-py[plot]"
        )

    adata = _read_h5ad(input)
    validate_analyzed(adata)

    limit = None if max_factors == 0 else max_factors
    out = output or default_output(input, suffix="_summary").replace(".h5ad", ".png")
    n_omitted = plot_summary(adata, out, dpi=dpi, max_factors=limit)
    if n_omitted:
        click.echo(
            f"Note: {n_omitted} factors omitted from plot "
            f"(showing top {max_factors}, ranked by FDR then gene count). "
            f"Use --max-factors 0 to show all."
        )
    click.echo(f"Wrote {out}")


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

@cli.command()
@click.argument("input", type=click.Path(exists=True))
@click.option("-p", "--phenotype", required=True, help="Condition column in adata.obs.")
@click.option("-o", "--output", default=None, help="Output HTML path.")
@click.option("-v", "--verbose", count=True, help="Increase verbosity (-v INFO, -vv DEBUG).")
@click.option("-q", "--quiet", is_flag=True, help="Errors only.")
def report(input, phenotype, output, verbose, quiet):
    """Generate an HTML report from CDR-g results."""
    _setup_logging(verbose, quiet)
    from .reporting import generate_html_report

    adata = _read_h5ad(input)
    validate_phenotype(adata, phenotype)
    validate_analyzed(adata)

    out = output or default_output(input, suffix="_report").replace(".h5ad", ".html")
    generate_html_report(adata, out, phenotype)
    click.echo(f"Wrote {out}")


# ---------------------------------------------------------------------------
# run  (full pipeline)
# ---------------------------------------------------------------------------

@cli.command()
@click.argument("input", type=click.Path(exists=True))
@click.option("-p", "--phenotype", required=True, help="Condition column in adata.obs.")
@click.option("-o", "--output", default=None, help="Output .h5ad path.")
@click.option("-c", "--csv", default=None, help="Export results CSV.")
@click.option("--capvar", default=0.95, show_default=True, help="Variance threshold.")
@click.option("--nperm-analysis", "nperm_analysis", default=2000, show_default=True,
              help="Permutations for importance scores.")
@click.option("--pernum", default=None, type=int, hidden=True, help="Deprecated alias for --nperm-analysis.")
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
@click.option("--correction", default="fdr_bh", show_default=True,
              type=click.Choice(["fdr_bh", "none"]), help="Multiple testing correction method.")
@click.option("-s", "--subset", multiple=True,
              help="Subset cells: COLUMN=VALUE[,VALUE2]. Repeatable.")
@click.option("--report", "report_path", default=None,
              help="Generate an HTML report at this path.")
@click.option("-v", "--verbose", count=True)
@click.option("-q", "--quiet", is_flag=True)
def run(input, phenotype, output, csv, capvar, nperm_analysis, pernum, thres,
        filter_method, cell_fraction, median_count, count_threshold, min_cells,
        enrich, enrich_method, genecol, nperm, enrich_thresh, seed, correction,
        subset, report_path, verbose, quiet):
    """Full CDR-g pipeline: filter, analyze, enrich, export."""
    _setup_logging(verbose, quiet)
    from .pycdr import run_CDR_analysis
    from .utils import (filter_genecounts_percent, filter_genecounts_numcells,
                        output_results)
    from .reporting import format_run_summary
    import warnings

    if pernum is not None:
        warnings.warn("--pernum is deprecated; use --nperm-analysis instead.", DeprecationWarning, stacklevel=2)
        nperm_analysis = pernum

    adata = _read_h5ad(input)
    adata = subset_cells(adata, subset)
    validate_phenotype(adata, phenotype)
    _validate_phenotype_after_subset(adata, phenotype)

    # --- filter ---
    if filter_method == "percent":
        logger.info("Filtering genes (percent method)")
        adata = filter_genecounts_percent(adata, cell_fraction, median_count)
        _echo(f"After filtering: {adata.shape[0]} cells x {adata.shape[1]} genes")
    elif filter_method == "numcells":
        logger.info("Filtering genes (numcells method)")
        adata = filter_genecounts_numcells(adata, count_threshold, min_cells)
        _echo(f"After filtering: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # --- analyze ---
    logger.info("Running CDR-g analysis")
    run_CDR_analysis(adata, phenotype, capvar=capvar, nperm=nperm_analysis, thres=thres,
                     seed=seed, correction=correction, quiet=quiet)
    n_factors = len(adata.uns["factor_loadings"])
    _echo(f"CDR-g analysis complete: {n_factors} factors")

    # --- store run parameters ---
    params = {
        "phenotype": phenotype,
        "capvar": capvar,
        "nperm": nperm_analysis,
        "thres": thres,
        "seed": seed,
        "correction": correction,
        "filter_method": filter_method,
        "subset": list(subset) if subset else [],
    }
    if filter_method == "percent":
        params["cell_fraction"] = cell_fraction
        params["median_count"] = median_count
    elif filter_method == "numcells":
        params["count_threshold"] = count_threshold
        params["min_cells"] = min_cells

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
            calculate_enrichment(adata, phenotype, factor_list, nperm, genecol, enrich_thresh, seed=seed, quiet=quiet)
        else:
            from .kruskal import calculate_enrichment as calc_kruskal

            logger.info("Running Kruskal-Wallis enrichment")
            calc_kruskal(adata, phenotype, quiet=quiet)

        _echo("Enrichment complete")

        params["enrich_method"] = enrich_method
        if enrich_method == "perm":
            params["enrich_nperm"] = nperm
            params["enrich_thresh"] = enrich_thresh
            params["genecol"] = genecol

    adata.uns["cdr_params"] = params

    # --- summary ---
    _echo(format_run_summary(adata, phenotype, enriched=enrich))

    # --- output ---
    out = output or default_output(input)
    adata.write(out)
    click.echo(f"Wrote {out}")

    if csv:
        df = output_results(adata)
        if df is not None:
            _serialize_genes_column(df).to_csv(csv, index=False)
            click.echo(f"Wrote {csv}")

    if report_path:
        from .reporting import generate_html_report

        generate_html_report(adata, report_path, phenotype)
        click.echo(f"Wrote {report_path}")


# ---------------------------------------------------------------------------
# demo
# ---------------------------------------------------------------------------

@cli.command()
@click.option("-o", "--output-dir", default="./pycdr_demo",
              show_default=True, help="Directory for demo output files.")
@click.option("-v", "--verbose", count=True, help="Increase verbosity (-v INFO, -vv DEBUG).")
@click.option("-q", "--quiet", is_flag=True, help="Errors only.")
def demo(output_dir, verbose, quiet):
    """Run an end-to-end example on bundled test data."""
    _setup_logging(verbose, quiet)

    import anndata as ad

    from .utils import filter_genecounts_numcells, output_results
    from .pycdr import run_CDR_analysis
    from .kruskal import calculate_enrichment as calc_kruskal
    from .reporting import generate_html_report, format_run_summary

    click.echo("pycdr demo: running end-to-end example...\n")

    # 1. Create output directory
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 2. Load bundled test data
    data_path = Path(__file__).parent / "data" / "test_muscle.h5ad"
    adata = ad.read_h5ad(str(data_path))
    phenotype = "Hours"

    n_cells, n_genes = adata.shape
    click.echo(f"Dataset: {n_cells} cells x {n_genes} genes (muscle differentiation time-course)")
    groups = sorted(adata.obs[phenotype].unique().tolist(), key=str)
    click.echo(f"Phenotype: '{phenotype}' ({', '.join(str(g) for g in groups)})\n")

    # 3. Filter genes
    click.echo("Filtering genes...")
    adata = filter_genecounts_numcells(adata, count_threshold=1, min_expressed_cells=10)
    click.echo(f"  After filtering: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # 4. CDR-g analysis
    click.echo("Running CDR-g analysis...")
    run_CDR_analysis(adata, phenotype, nperm=500, quiet=True)
    n_factors = len(adata.uns["factor_loadings"])
    click.echo(f"  {n_factors} factors identified")

    # 5. Kruskal-Wallis enrichment
    click.echo("Running Kruskal-Wallis enrichment...")
    calc_kruskal(adata, phenotype, quiet=True)
    click.echo("  Enrichment complete\n")

    # Store params for the report
    adata.uns["cdr_params"] = {
        "phenotype": phenotype,
        "capvar": 0.95,
        "nperm": 500,
        "thres": 0.05,
        "seed": 42,
        "filter_method": "numcells",
        "count_threshold": 1,
        "min_cells": 10,
        "enrich_method": "kruskal",
    }

    # 6. Save analyzed.h5ad + results.csv
    h5ad_path = outdir / "analyzed.h5ad"
    adata.write(str(h5ad_path))

    csv_path = outdir / "results.csv"
    df = output_results(adata)
    if df is not None:
        _serialize_genes_column(df).to_csv(str(csv_path), index=False)

    # 7. HTML report
    report_path = outdir / "report.html"
    generate_html_report(adata, str(report_path), phenotype)

    # 8. Summary plot (optional)
    png_path = outdir / "summary.png"
    try:
        from .plotting import plot_summary
        plot_summary(adata, str(png_path))
        plot_ok = True
    except ImportError:
        plot_ok = False

    # 9. Print summary
    click.echo(f"Results written to {output_dir}/")
    click.echo(f"  analyzed.h5ad    — annotated dataset with CDR-g results")
    click.echo(f"  results.csv      — factor summary table")
    click.echo(f"  report.html      — HTML report")
    if plot_ok:
        click.echo(f"  summary.png      — summary figure")
    click.echo(f"\nOpen report.html in a browser to explore the results.")
