"""Terminal formatting, factor summaries, and HTML report generation for pycdr."""

import html as _html
import io
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def _pretty_factor_name(name):
    """Convert internal factor key to a display name for HTML reports.

    ``"factor.0"`` becomes ``"Factor 0"``, ``"factor.12"`` becomes ``"Factor 12"``.
    Non-matching names are returned as-is.
    """
    if isinstance(name, str) and name.startswith("factor."):
        try:
            return f"Factor {int(name.split('.')[-1])}"
        except ValueError:
            pass
    return name


def _separator(char="=", width=60):
    """Return a horizontal rule string."""
    return char * width


def _format_table(headers, rows, min_col_width=6):
    """Return a simple ASCII-aligned table string.

    Parameters
    ----------
    headers : list[str]
    rows : list[list[str]]
    min_col_width : int
        Minimum column width (default 6).

    Returns
    -------
    str
        Multi-line table with header underline.
    """
    n_cols = len(headers)
    widths = [max(min_col_width, len(h)) for h in headers]
    for row in rows:
        for i, cell in enumerate(row[:n_cols]):
            widths[i] = max(widths[i], len(str(cell)))

    def _row(cells):
        parts = []
        for i, cell in enumerate(cells[:n_cols]):
            parts.append(str(cell).ljust(widths[i]))
        return "  ".join(parts)

    lines = [_row(headers), _row(["-" * w for w in widths])]
    for row in rows:
        lines.append(_row(row))
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Factor summary
# ---------------------------------------------------------------------------

def get_factor_summary(adata):
    """Build a summary DataFrame from CDR-g results stored in *adata.uns*.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object after ``run_CDR_analysis`` (and optionally enrichment).

    Returns
    -------
    pandas.DataFrame
        Columns: ``factor``, ``n_genes``, ``top_genes``, ``mean_zscore``,
        and optionally ``enrich_fdr``.
    """
    factor_loadings = adata.uns.get("factor_loadings", {})
    zscores = adata.uns.get("zscores", None)
    selection = adata.uns.get("selection", None)
    gene_names = adata.var_names.tolist()

    # Build a {gene: index} lookup for O(1) access
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}

    rows = []
    for fname, genes in factor_loadings.items():
        n_genes = len(genes)

        # Compute mean z-score for genes in this factor
        mean_z = np.nan
        top = []
        if zscores is not None and n_genes > 0:
            # Factor index from name like "factor.3"
            try:
                fi = int(fname.split(".")[-1])
            except (ValueError, IndexError):
                fi = None

            if fi is not None and fi < zscores.shape[1]:
                # Get z-scores for genes in this factor
                idxs = [gene_to_idx[g] for g in genes if g in gene_to_idx]
                if idxs:
                    z_vals = zscores[idxs, fi]
                    mean_z = float(np.mean(z_vals))
                    # Top 3 by z-score
                    order = np.argsort(z_vals)[::-1]
                    top = [genes[j] for j in order[:3] if j < len(genes)]

        row = {
            "factor": fname,
            "n_genes": n_genes,
            "top_genes": ", ".join(top) if top else "",
            "mean_zscore": round(mean_z, 2) if not np.isnan(mean_z) else np.nan,
        }

        # Add enrichment FDR and dominant condition if available
        enrich_stats = adata.uns.get("enrichment_stats", None)
        enrich_dominant = None
        if enrich_stats is not None and isinstance(enrich_stats, pd.DataFrame):
            if fname in enrich_stats.index:
                fdr_col = "fdr" if "fdr" in enrich_stats.columns else None
                if fdr_col:
                    row["enrich_fdr"] = enrich_stats.loc[fname, fdr_col]
                if "dominant_condition" in enrich_stats.columns:
                    enrich_dominant = enrich_stats.loc[fname, "dominant_condition"]

        # Dominant condition: prefer enrichment-derived, fall back to Fs-based
        fs_dominant = adata.uns.get("dominant_condition")
        if enrich_dominant is not None:
            row["dominant_condition"] = str(enrich_dominant)
        elif fs_dominant is not None:
            try:
                fi = int(fname.split(".")[-1])
                if fi < len(fs_dominant):
                    row["dominant_condition"] = str(fs_dominant[fi])
            except (ValueError, IndexError):
                pass

        rows.append(row)

    df = pd.DataFrame(rows)
    return df


# ---------------------------------------------------------------------------
# Terminal summary after `pycdr run`
# ---------------------------------------------------------------------------

def format_run_summary(adata, phenotype, enriched=False):
    """Return a multi-line terminal summary string for display after ``pycdr run``.

    Parameters
    ----------
    adata : anndata.AnnData
    phenotype : str
    enriched : bool
        Whether enrichment was run.

    Returns
    -------
    str
    """
    lines = []
    sep = _separator("=", 60)
    lines.append("")
    lines.append(sep)
    lines.append("  CDR-g Analysis Summary")
    lines.append(sep)
    lines.append("")

    n_cells, n_genes = adata.shape
    lines.append(f"  Dataset:    {n_cells} cells x {n_genes} genes")

    # Phenotype info
    counts = adata.obs[phenotype].value_counts()
    n_groups = len(counts)
    lines.append(f"  Phenotype:  '{phenotype}' ({n_groups} groups)")

    # Factor info
    fl = adata.uns.get("factor_loadings", {})
    n_factors_total = adata.uns.get("selected_loading", len(fl))
    factors_with_genes = sum(1 for v in fl.values() if len(v) > 0)
    lines.append(f"  Factors:    {n_factors_total} ({factors_with_genes} with genes)")
    lines.append("")

    # Summary table
    df = get_factor_summary(adata)
    # Only show factors with genes
    df = df[df["n_genes"] > 0]

    if len(df) > 0:
        headers = ["Factor", "Genes", "Mean Z"]
        has_fdr = "enrich_fdr" in df.columns
        has_condition = "dominant_condition" in df.columns
        if has_condition:
            headers.append("Condition")
        if has_fdr:
            headers.append("FDR")
        headers.append("Top genes")

        rows = []
        for _, r in df.iterrows():
            row = [r["factor"], str(r["n_genes"]), f"{r['mean_zscore']:.2f}"]
            if has_condition:
                row.append(str(r.get("dominant_condition", "")))
            if has_fdr:
                fdr = r.get("enrich_fdr", np.nan)
                row.append(f"{fdr:.2e}" if not pd.isna(fdr) else "n/a")
            row.append(r["top_genes"])
            rows.append(row)

        table = _format_table(headers, rows)
        for line in table.split("\n"):
            lines.append(f"  {line}")
    else:
        lines.append("  No factors with significant genes.")

    lines.append("")
    lines.append(sep)
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Enhanced info summary
# ---------------------------------------------------------------------------

def format_info_summary(adata):
    """Return a summary string for the ``pycdr info`` command when analysis exists.

    Parameters
    ----------
    adata : anndata.AnnData

    Returns
    -------
    str or None
        ``None`` if no CDR-g analysis found.
    """
    fl = adata.uns.get("factor_loadings")
    if fl is None:
        return None

    lines = []
    lines.append("")
    lines.append("  CDR-g Analysis Results")
    lines.append("  " + _separator("-", 40))

    n_factors = adata.uns.get("selected_loading", len(fl))
    lines.append(f"  Factors:         {n_factors}")

    # Gene hit counts
    all_genes = []
    for genes in fl.values():
        all_genes.extend(genes)
    total_hits = len(all_genes)
    unique_hits = len(set(all_genes))
    lines.append(f"  Total gene hits: {total_hits} ({unique_hits} unique)")
    lines.append("")

    # Per-factor gene count table
    headers = ["Factor", "Genes"]
    rows = []
    for fname, genes in fl.items():
        rows.append([fname, str(len(genes))])
    table = _format_table(headers, rows)
    for line in table.split("\n"):
        lines.append(f"  {line}")

    # Enrichment significance count
    enrich_stats = adata.uns.get("enrichment_stats")
    if enrich_stats is not None and isinstance(enrich_stats, pd.DataFrame):
        fdr_col = "fdr" if "fdr" in enrich_stats.columns else None
        if fdr_col:
            n_sig = (enrich_stats[fdr_col] < 0.05).sum()
            n_total = len(enrich_stats)
            lines.append("")
            lines.append(f"  Enrichment: {n_sig}/{n_total} factors significant (FDR<0.05)")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Results table formatting
# ---------------------------------------------------------------------------

def format_results_table(df, fmt="csv"):
    """Format a results DataFrame for output.

    Parameters
    ----------
    df : pandas.DataFrame
        From ``output_results()``.
    fmt : str
        One of ``"csv"``, ``"tsv"``, ``"table"``, ``"markdown"``.

    Returns
    -------
    str
    """
    if fmt in ("csv", "tsv"):
        sep = "\t" if fmt == "tsv" else ","
        # Serialize genes column if present as list
        out_df = df.copy()
        if "genes" in out_df.columns:
            out_df["genes"] = out_df["genes"].apply(
                lambda x: ",".join(x) if isinstance(x, list) else str(x)
            )
        return out_df.to_csv(sep=sep, index=False)

    # For table/markdown, show summary columns only (drop full genes list)
    display_cols = [c for c in df.columns if c != "genes"]
    out_df = df[display_cols].copy()

    if fmt == "markdown":
        return _format_markdown(out_df)
    else:
        return _format_ascii_table(out_df)


def _format_markdown(df):
    """Format DataFrame as a markdown table."""
    cols = df.columns.tolist()
    lines = []
    lines.append("| " + " | ".join(cols) + " |")
    lines.append("| " + " | ".join(["---"] * len(cols)) + " |")
    for _, row in df.iterrows():
        cells = []
        for c in cols:
            val = row[c]
            if isinstance(val, float):
                if abs(val) < 0.01 and val != 0:
                    cells.append(f"{val:.2e}")
                else:
                    cells.append(f"{val:.4f}" if not pd.isna(val) else "")
            else:
                cells.append(str(val))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def _format_ascii_table(df):
    """Format DataFrame as an ASCII aligned table."""
    cols = df.columns.tolist()
    rows = []
    for _, row in df.iterrows():
        cells = []
        for c in cols:
            val = row[c]
            if isinstance(val, float):
                if abs(val) < 0.01 and val != 0:
                    cells.append(f"{val:.2e}")
                else:
                    cells.append(f"{val:.4f}" if not pd.isna(val) else "")
            else:
                cells.append(str(val))
        rows.append(cells)
    return _format_table(cols, rows)


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------

_HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>CDR-g Report</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
         max-width: 900px; margin: 40px auto; padding: 0 20px; color: #333; }}
  h1 {{ border-bottom: 2px solid #2563eb; padding-bottom: 8px; }}
  h2 {{ color: #2563eb; margin-top: 2em; }}
  h3 {{ color: #475569; margin-top: 1.5em; }}
  table {{ border-collapse: collapse; width: 100%; margin: 1em 0; }}
  th, td {{ border: 1px solid #ddd; padding: 8px 12px; text-align: left; }}
  th {{ background: #f8fafc; font-weight: 600; }}
  tr:nth-child(even) {{ background: #f8fafc; }}
  .meta {{ color: #666; }}
  .note {{ background: #fffbeb; border: 1px solid #fbbf24; padding: 12px; border-radius: 4px; }}
  .guide {{ background: #f0f9ff; border: 1px solid #bae6fd; border-radius: 6px;
            padding: 16px 20px; margin: 1.5em 0; line-height: 1.6; }}
  .guide p {{ margin: 0.5em 0; }}
  .guide dl {{ margin: 0.5em 0 0.5em 0; }}
  .guide dt {{ font-weight: 600; margin-top: 0.5em; }}
  .guide dd {{ margin: 0.1em 0 0 1.5em; }}
  img.figure {{ max-width: 100%; border: 1px solid #ddd; border-radius: 4px; margin: 1em 0; }}
</style>
</head>
<body>
<h1>CDR-g Analysis Report</h1>

<h2>Dataset Overview</h2>
<p class="meta">
  <strong>Shape:</strong> {n_cells} cells &times; {n_genes} genes<br>
  <strong>Phenotype:</strong> '{phenotype}' ({n_groups} groups: {group_detail})<br>
  <strong>Factors:</strong> {n_factors} ({factors_with_genes} with genes)
</p>

{params_section}

<h2>How to read this report</h2>
<div class="guide">
<p>
CDR-g decomposes gene co-expression patterns across your experimental conditions
using SVD and varimax rotation. Each <strong>factor</strong> represents a
co-expression program &mdash; a group of genes whose expression changes together
between conditions.
</p>

<h3>Factor summary table</h3>
<dl>
  <dt>Factor</dt>
  <dd>A numbered identifier (e.g. Factor 0, Factor 1). Factors are ordered by
  eigenvalue, so lower-numbered factors capture more variance. A factor with
  <strong>0 genes</strong> showed no statistically significant condition-dependent
  co-expression and can usually be ignored. In the underlying AnnData object,
  factors are stored as <code>factor.0</code>, <code>factor.1</code>, etc.</dd>

  <dt>n_genes</dt>
  <dd>The number of genes assigned to this factor by permutation testing
  (p&nbsp;&lt;&nbsp;0.05). More genes suggests a broader transcriptional program;
  fewer genes suggests a more specific one.</dd>

  <dt>mean_zscore</dt>
  <dd>The average z-score of gene loading differences across conditions for
  genes in the factor. Higher values indicate stronger condition-dependent
  co-expression. A z-score &gt;&nbsp;2 indicates a gene&rsquo;s loading
  difference is well above what would be expected by chance.</dd>

  <dt>top_genes</dt>
  <dd>The three genes with the highest z-scores in each factor. These are the
  most condition-responsive genes in the program and are a good starting point
  for biological interpretation (e.g. pathway lookup, literature search).</dd>
</dl>

{enrichment_guide}

<h3>Practical interpretation</h3>
<p>
Factors with many genes and high mean z-scores represent the dominant
condition-specific programs in your data. To interpret a factor biologically,
examine its top genes for known pathway membership or functional annotations.
Factors that share genes may represent overlapping or related programs.
Factors with few genes (&lt;&nbsp;5) may reflect noise or very specific
regulatory events &mdash; consider them with caution unless the enrichment
test confirms significance.
</p>
</div>

<h2>Factor Summary</h2>
{factor_table_html}

{enrichment_section}

{figure_section}

<hr>
<p class="meta" style="font-size:0.85em;">Generated by pycdr</p>
</body>
</html>
"""


def generate_html_report(adata, output_path, phenotype):
    """Generate a self-contained HTML report for CDR-g results.

    Parameters
    ----------
    adata : anndata.AnnData
    output_path : str or Path
    phenotype : str
    """
    n_cells, n_genes = adata.shape
    counts = adata.obs[phenotype].value_counts()
    n_groups = len(counts)
    group_detail = ", ".join(f"{k}={v}" for k, v in counts.items())

    fl = adata.uns.get("factor_loadings", {})
    n_factors = adata.uns.get("selected_loading", len(fl))
    factors_with_genes = sum(1 for v in fl.values() if len(v) > 0)

    # Factor summary table
    df = get_factor_summary(adata)
    df["factor"] = df["factor"].map(_pretty_factor_name)
    factor_table_html = _df_to_html_table(df)

    # Enrichment section and guide
    enrich_stats = adata.uns.get("enrichment_stats")
    has_enrichment = enrich_stats is not None and isinstance(enrich_stats, pd.DataFrame)
    if has_enrichment:
        enrich_display = enrich_stats.reset_index().rename(columns={"index": "factor"})
        enrich_display["factor"] = enrich_display["factor"].map(_pretty_factor_name)
        enrichment_section = "<h2>Enrichment Statistics</h2>\n" + _df_to_html_table(
            enrich_display
        )
        # Determine which method was used from column names
        if "stat" in enrich_stats.columns and "pval" in enrich_stats.columns:
            method_note = (
                "The Kruskal-Wallis test was used to assess whether each "
                "factor&#8217;s gene set is differentially enriched across "
                "conditions."
            )
            stat_label = "stat"
            stat_desc = (
                "The Kruskal-Wallis H-statistic. Larger values indicate "
                "stronger differences in enrichment scores between conditions."
            )
        else:
            method_note = (
                "A permutation-based enrichment test was used to assess "
                "whether each factor&#8217;s gene set shows "
                "condition-dependent activation."
            )
            stat_label = "statistic"
            stat_desc = (
                "The chi-square test statistic for differential activation "
                "across conditions."
            )

        enrichment_guide = (
            "<h3>Enrichment statistics</h3>\n"
            f"<p>{method_note}</p>\n"
            "<dl>\n"
            f"  <dt>{stat_label}</dt>\n"
            f"  <dd>{stat_desc}</dd>\n"
            "  <dt>pval</dt>\n"
            "  <dd>The raw p-value from the statistical test. Smaller values "
            "indicate more significant differences between conditions.</dd>\n"
            "  <dt>FDR</dt>\n"
            "  <dd>The Benjamini-Hochberg false discovery rate. Factors with "
            "FDR&nbsp;&lt;&nbsp;0.05 are considered significantly enriched "
            "across conditions. These are the most biologically meaningful "
            "factors to investigate further.</dd>\n"
            "  <dt>dominant_condition</dt>\n"
            "  <dd>The phenotype group that drives this factor most strongly. "
            "When enrichment is available, this is determined by the condition "
            "with the highest median ssGSEA score; otherwise it is derived "
            "from the condition with the highest mean absolute factor "
            "loading.</dd>\n"
            "</dl>"
        )
    else:
        enrichment_section = ""
        enrichment_guide = (
            '<p><em>No enrichment analysis was run. Use '
            '<code>--enrich</code> with <code>pycdr run</code> or '
            '<code>pycdr enrich</code> to test which factors show '
            'statistically significant differences between conditions.</em></p>'
        )

    # Figure section
    figure_section = _try_embed_figure(adata)

    # Run parameters section
    params_section = _build_params_section(adata)

    html = _HTML_TEMPLATE.format(
        n_cells=n_cells,
        n_genes=n_genes,
        phenotype=_html.escape(phenotype),
        n_groups=n_groups,
        group_detail=_html.escape(group_detail),
        n_factors=n_factors,
        factors_with_genes=factors_with_genes,
        params_section=params_section,
        factor_table_html=factor_table_html,
        enrichment_guide=enrichment_guide,
        enrichment_section=enrichment_section,
        figure_section=figure_section,
    )

    Path(output_path).write_text(html, encoding="utf-8")


def _build_params_section(adata):
    """Build an HTML fragment showing the run parameters from ``adata.uns["cdr_params"]``."""
    params = adata.uns.get("cdr_params")
    if params is None:
        return ""

    # Analysis parameters
    rows = []
    _add = lambda label, key, fmt=str: rows.append(
        (label, fmt(params[key])) if key in params else None
    )

    rows.append(("Phenotype", _html.escape(str(params.get("phenotype", "n/a")))))
    rows.append(("Variance threshold (capvar)", str(params.get("capvar", "n/a"))))
    rows.append(("Permutations (pernum)", str(params.get("pernum", "n/a"))))
    rows.append(("P-value threshold (thres)", str(params.get("thres", "n/a"))))

    # Subset — may be a list or numpy array after h5ad round-trip
    subset = params.get("subset", [])
    subset_list = list(subset) if hasattr(subset, '__iter__') and not isinstance(subset, str) else []
    if len(subset_list) > 0:
        rows.append(("Cell subset", _html.escape(", ".join(str(s) for s in subset_list))))

    # Filter parameters
    fm = str(params.get("filter_method", "none"))
    if fm and fm != "none":
        rows.append(("Filter method", _html.escape(fm)))
        if fm == "percent":
            rows.append(("  cell_fraction", str(params.get("cell_fraction", ""))))
            rows.append(("  median_count", str(params.get("median_count", ""))))
        elif fm == "numcells":
            rows.append(("  count_threshold", str(params.get("count_threshold", ""))))
            rows.append(("  min_cells", str(params.get("min_cells", ""))))

    # Enrichment parameters
    em = params.get("enrich_method")
    if em:
        rows.append(("Enrichment method", _html.escape(str(em))))
        if str(em) == "perm":
            rows.append(("  nperm", str(params.get("nperm", ""))))
            rows.append(("  enrich_thresh", str(params.get("enrich_thresh", ""))))
            gc = params.get("genecol")
            if gc:
                rows.append(("  genecol", _html.escape(str(gc))))
            rows.append(("  seed", str(params.get("seed", ""))))

    lines = ['<h2>Run Parameters</h2>', "<table>"]
    for label, value in rows:
        lines.append(f"<tr><td><strong>{label}</strong></td><td>{value}</td></tr>")
    lines.append("</table>")
    return "\n".join(lines)


def _df_to_html_table(df):
    """Convert a DataFrame to an HTML table string."""
    lines = ["<table>", "<thead><tr>"]
    for col in df.columns:
        lines.append(f"  <th>{_html.escape(str(col))}</th>")
    lines.append("</tr></thead>")
    lines.append("<tbody>")
    for _, row in df.iterrows():
        lines.append("<tr>")
        for col in df.columns:
            val = row[col]
            if isinstance(val, float):
                if pd.isna(val):
                    cell = ""
                elif abs(val) < 0.01 and val != 0:
                    cell = f"{val:.2e}"
                else:
                    cell = f"{val:.4f}"
            elif isinstance(val, list):
                cell = _html.escape(", ".join(str(x) for x in val))
            else:
                cell = _html.escape(str(val))
            lines.append(f"  <td>{cell}</td>")
        lines.append("</tr>")
    lines.append("</tbody></table>")
    return "\n".join(lines)


def _try_embed_figure(adata):
    """Try to generate and embed a summary figure as base64 PNG.

    Returns an HTML fragment or a note if matplotlib is not available.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        from .plotting import plot_summary
        import base64

        buf = io.BytesIO()
        plot_summary(adata, buf, dpi=100)
        buf.seek(0)
        b64 = base64.b64encode(buf.read()).decode("ascii")
        return (
            '<h2>Summary Figure</h2>\n'
            f'<img class="figure" src="data:image/png;base64,{b64}" '
            f'alt="CDR-g summary figure">'
        )
    except ImportError:
        return (
            '<p class="note">Install matplotlib for embedded plots: '
            '<code>pip install cdr-py[plot]</code></p>'
        )
    except Exception:
        return ""
