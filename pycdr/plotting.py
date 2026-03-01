"""Matplotlib-based summary figures for pycdr CDR-g results."""

import logging

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as _exc:
    raise ImportError(
        "matplotlib is required for pycdr plotting. "
        "Install with: pip install cdr-py[plot]"
    ) from _exc

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

DEFAULT_MAX_FACTORS = 30


def _rank_factors(fl, enrich_stats):
    """Rank factors by enrichment FDR (ascending), then gene count (descending).

    Returns an ordered list of factor names.
    """
    factor_names = list(fl.keys())
    has_fdr = (
        enrich_stats is not None
        and isinstance(enrich_stats, pd.DataFrame)
        and "fdr" in enrich_stats.columns
    )

    def sort_key(fname):
        fdr = enrich_stats.loc[fname, "fdr"] if has_fdr and fname in enrich_stats.index else 1.0
        gene_count = len(fl.get(fname, []))
        # Sort by FDR ascending, then gene count descending
        return (fdr, -gene_count)

    return sorted(factor_names, key=sort_key)


def plot_summary(adata, output_path, dpi=150, max_factors=DEFAULT_MAX_FACTORS):
    """Generate a multi-panel summary figure from CDR-g results.

    Panels:
    (a) Horizontal bar chart of gene counts per factor.
    (b) Heatmap of top genes x factors (z-scores).
    (c) Horizontal bar chart of -log10(FDR) per factor (if enrichment available).

    When there are more factors than *max_factors*, only the top factors
    are shown, ranked by enrichment FDR (if available) then gene count.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData after ``run_CDR_analysis`` (and optionally enrichment).
    output_path : str, Path, or file-like
        Where to save the figure (PNG format).
    dpi : int
        Figure resolution.
    max_factors : int or None
        Maximum number of factors to display. Set to ``None`` to show all.
        Defaults to 30.

    Returns
    -------
    int or None
        Number of factors omitted from the plot, or ``None`` if all are shown.
    """
    fl = adata.uns.get("factor_loadings", {})
    zscores = adata.uns.get("zscores")
    enrich_stats = adata.uns.get("enrichment_stats")
    gene_names = adata.var_names.tolist()

    # --- data prep ---
    all_factor_names = list(fl.keys())
    total_factors = len(all_factor_names)

    has_enrichment = (
        enrich_stats is not None
        and isinstance(enrich_stats, pd.DataFrame)
        and "fdr" in enrich_stats.columns
    )

    # Limit displayed factors
    n_omitted = None
    if max_factors is not None and total_factors > max_factors:
        ranked = _rank_factors(fl, enrich_stats)
        factor_names = ranked[:max_factors]
        n_omitted = total_factors - max_factors
        logger.info(
            "Showing top %d of %d factors (ranked by %s). "
            "%d factors omitted. Use --max-factors to adjust.",
            max_factors, total_factors,
            "FDR then gene count" if has_enrichment else "gene count",
            n_omitted,
        )
    else:
        factor_names = all_factor_names

    gene_counts = [len(fl[f]) for f in factor_names]

    n_panels = 3 if has_enrichment else 2

    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, max(4, len(factor_names) * 0.35)))
    if n_panels == 1:
        axes = [axes]

    # --- Panel (a): gene count bar chart ---
    ax = axes[0]
    y_pos = np.arange(len(factor_names))
    ax.barh(y_pos, gene_counts, color="#2563eb", edgecolor="white")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(factor_names, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Gene count")
    title_a = "(a) Genes per factor"
    if n_omitted:
        title_a += f" (top {len(factor_names)} of {total_factors})"
    ax.set_title(title_a)

    # --- Panel (b): z-score heatmap ---
    ax = axes[1]
    if zscores is not None:
        gene_to_idx = {g: i for i, g in enumerate(gene_names)}

        # Collect top genes across factors (up to 20 unique)
        top_gene_set = []
        seen = set()
        for fname in factor_names:
            genes = fl.get(fname, [])
            try:
                fi = int(fname.split(".")[-1])
            except (ValueError, IndexError):
                continue
            if fi >= zscores.shape[1]:
                continue
            idxs = [gene_to_idx[g] for g in genes if g in gene_to_idx]
            if not idxs:
                continue
            z_vals = zscores[idxs, fi]
            order = np.argsort(z_vals)[::-1]
            for j in order:
                g = genes[j] if j < len(genes) else None
                if g and g not in seen:
                    top_gene_set.append(g)
                    seen.add(g)
                if len(top_gene_set) >= 20:
                    break
            if len(top_gene_set) >= 20:
                break

        if top_gene_set:
            # Build heatmap matrix: genes x factors
            heatmap = np.zeros((len(top_gene_set), len(factor_names)))
            for j, fname in enumerate(factor_names):
                try:
                    fi = int(fname.split(".")[-1])
                except (ValueError, IndexError):
                    continue
                if fi >= zscores.shape[1]:
                    continue
                for i, g in enumerate(top_gene_set):
                    idx = gene_to_idx.get(g)
                    if idx is not None and idx < zscores.shape[0]:
                        heatmap[i, j] = zscores[idx, fi]

            vmax = max(abs(heatmap.max()), abs(heatmap.min()), 1)
            im = ax.imshow(heatmap, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
            ax.set_yticks(np.arange(len(top_gene_set)))
            ax.set_yticklabels(top_gene_set, fontsize=7)
            ax.set_xticks(np.arange(len(factor_names)))
            ax.set_xticklabels(factor_names, fontsize=7, rotation=45, ha="right")
            fig.colorbar(im, ax=ax, shrink=0.6, label="z-score")
        else:
            ax.text(0.5, 0.5, "No genes", ha="center", va="center", transform=ax.transAxes)
    else:
        ax.text(0.5, 0.5, "No z-scores", ha="center", va="center", transform=ax.transAxes)
    ax.set_title("(b) Top gene z-scores")

    # --- Panel (c): enrichment FDR ---
    if has_enrichment:
        ax = axes[2]
        fdr_vals = []
        for fname in factor_names:
            if fname in enrich_stats.index:
                fdr_vals.append(enrich_stats.loc[fname, "fdr"])
            else:
                fdr_vals.append(np.nan)
        neg_log_fdr = [-np.log10(f) if f > 0 else 0 for f in fdr_vals]

        colors = ["#dc2626" if f < 0.05 else "#94a3b8" for f in fdr_vals]
        ax.barh(y_pos, neg_log_fdr, color=colors, edgecolor="white")
        ax.axvline(-np.log10(0.05), color="black", linestyle="--", linewidth=0.8, label="FDR=0.05")
        ax.set_yticks(y_pos)
        ax.set_yticklabels(factor_names, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("-log10(FDR)")
        ax.set_title("(c) Enrichment significance")
        ax.legend(fontsize=7, loc="lower right")

    plt.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight", format="png")
    plt.close(fig)

    return n_omitted
