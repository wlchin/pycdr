"""Generate 2x3-panel permutation benchmark figure.

Panels:
  A: perm_time vs nperm (lines per n_genes) — confirms linearity
  B: per_perm_ms vs n_genes (lines per nperm) — per-iteration cost scaling
  C: perm_frac bar chart at nperm=2000 — the decision plot
  D: perm_time heatmap (n_genes x nperm) — absolute time overview
  E: perm_frac heatmap (n_genes x nperm) — fraction overview
  F: SVD vs Perm stacked bars at nperm=2000 — direct visual comparison
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


PROD_NPERM = 2000  # Production default for decision panels


def plot_permutation(summary_path, output_path):
    df = pd.read_csv(summary_path, sep="\t")

    # Compute medians across repeats
    group_cols = ["n_genes", "nperm"]
    median = df.groupby(group_cols).median(numeric_only=True).reset_index()

    # Derive columns if not present
    if "perm_fraction" not in median.columns:
        median["perm_fraction"] = median["perm_time_s"] / median["total_time_s"]
    if "per_perm_ms" not in median.columns:
        median["per_perm_ms"] = (median["perm_time_s"] / median["nperm"]) * 1000

    # Rename for consistency (Snakemake data uses _s suffix, standalone doesn't)
    col_map = {
        "perm_time_s": "perm_time",
        "svd_time_s": "svd_time",
        "total_time_s": "total_time",
    }
    for old, new in col_map.items():
        if old in median.columns and new not in median.columns:
            median[new] = median[old]

    genes_vals = sorted(median["n_genes"].unique())
    nperm_vals = sorted(median["nperm"].unique())

    fig, axes = plt.subplots(2, 3, figsize=(18, 11))

    # ==================== Panel A: perm_time vs nperm ====================
    ax = axes[0, 0]
    for ng in genes_vals:
        sub = median[median["n_genes"] == ng].sort_values("nperm")
        ax.plot(sub["nperm"], sub["perm_time"], marker="o", label=f"{ng} genes")
    ax.set_xlabel("nperm")
    ax.set_ylabel("Permutation time (s)")
    ax.set_title("A) Perm time vs nperm")
    ax.legend(title="n_genes", fontsize=8)
    ax.grid(True, alpha=0.3)

    # ==================== Panel B: per_perm_ms vs n_genes ====================
    ax = axes[0, 1]
    for np_val in nperm_vals:
        sub = median[median["nperm"] == np_val].sort_values("n_genes")
        ax.plot(sub["n_genes"], sub["per_perm_ms"], marker="s", label=f"nperm={np_val}")
    ax.set_xlabel("n_genes")
    ax.set_ylabel("Time per permutation (ms)")
    ax.set_title("B) Per-perm cost vs gene count")
    ax.legend(title="nperm", fontsize=8)
    ax.grid(True, alpha=0.3)

    # ==================== Panel C: perm_frac bars at nperm=2000 ====================
    ax = axes[0, 2]
    ref_nperm = PROD_NPERM
    # Fall back to closest available nperm if exact value missing
    if ref_nperm not in median["nperm"].values:
        ref_nperm = min(nperm_vals, key=lambda x: abs(x - PROD_NPERM))
    sub = median[median["nperm"] == ref_nperm].sort_values("n_genes")
    x = np.arange(len(sub))
    bars = ax.bar(x, sub["perm_fraction"] * 100, color="#DD8452", edgecolor="black", linewidth=0.5)
    ax.axhline(y=10, color="red", linestyle="--", alpha=0.7, label="10% threshold")
    ax.set_xticks(x)
    ax.set_xticklabels(sub["n_genes"].astype(int))
    ax.set_xlabel("n_genes")
    ax.set_ylabel("Permutation fraction (%)")
    ax.set_title(f"C) Perm fraction at nperm={ref_nperm}")
    ax.legend()
    ax.grid(True, alpha=0.3, axis="y")
    # Annotate bars
    for bar, val in zip(bars, sub["perm_fraction"]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                f"{val:.1%}", ha="center", va="bottom", fontsize=8)

    # ==================== Panel D: perm_time heatmap ====================
    ax = axes[1, 0]
    pivot_time = median.pivot(index="n_genes", columns="nperm", values="perm_time")
    pivot_time = pivot_time.reindex(index=genes_vals, columns=nperm_vals)

    vmin = max(pivot_time.values[~np.isnan(pivot_time.values)].min(), 0.01)
    vmax = pivot_time.values[~np.isnan(pivot_time.values)].max()
    im = ax.imshow(pivot_time.values, norm=LogNorm(vmin=vmin, vmax=vmax),
                   cmap="YlOrRd", aspect="auto")
    ax.set_xticks(range(len(nperm_vals)))
    ax.set_xticklabels(nperm_vals)
    ax.set_yticks(range(len(genes_vals)))
    ax.set_yticklabels(genes_vals)
    ax.set_xlabel("nperm")
    ax.set_ylabel("n_genes")
    ax.set_title("D) Perm time heatmap (s)")
    for i in range(len(genes_vals)):
        for j in range(len(nperm_vals)):
            val = pivot_time.values[i, j]
            if not np.isnan(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if val > vmax * 0.3 else "black")
    fig.colorbar(im, ax=ax, label="seconds")

    # ==================== Panel E: perm_frac heatmap ====================
    ax = axes[1, 1]
    pivot_frac = median.pivot(index="n_genes", columns="nperm", values="perm_fraction")
    pivot_frac = pivot_frac.reindex(index=genes_vals, columns=nperm_vals)

    im2 = ax.imshow(pivot_frac.values * 100, cmap="RdYlGn_r", aspect="auto",
                    vmin=0, vmax=max(pivot_frac.values[~np.isnan(pivot_frac.values)].max() * 100, 1))
    ax.set_xticks(range(len(nperm_vals)))
    ax.set_xticklabels(nperm_vals)
    ax.set_yticks(range(len(genes_vals)))
    ax.set_yticklabels(genes_vals)
    ax.set_xlabel("nperm")
    ax.set_ylabel("n_genes")
    ax.set_title("E) Perm fraction heatmap (%)")
    for i in range(len(genes_vals)):
        for j in range(len(nperm_vals)):
            val = pivot_frac.values[i, j]
            if not np.isnan(val):
                ax.text(j, i, f"{val:.1%}", ha="center", va="center", fontsize=7)
    fig.colorbar(im2, ax=ax, label="%")

    # ==================== Panel F: SVD vs Perm stacked bars at nperm=2000 ====================
    ax = axes[1, 2]
    sub = median[median["nperm"] == ref_nperm].sort_values("n_genes")
    x = np.arange(len(sub))
    width = 0.6
    ax.bar(x, sub["svd_time"], width, label="SVD", color="#4C72B0")
    ax.bar(x, sub["perm_time"], width, bottom=sub["svd_time"],
           label="Permutation", color="#DD8452")
    ax.set_xticks(x)
    ax.set_xticklabels(sub["n_genes"].astype(int))
    ax.set_xlabel("n_genes")
    ax.set_ylabel("Time (s)")
    ax.set_title(f"F) SVD vs Perm (nperm={ref_nperm})")
    ax.legend()
    ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    # Snakemake script interface
    plot_permutation(
        summary_path=str(snakemake.input[0]),
        output_path=str(snakemake.output[0]),
    )
