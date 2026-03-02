"""Generate 2x3-panel benchmark scaling figure (timing + memory)."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def plot_benchmark(summary_path, output_path):
    df = pd.read_csv(summary_path, sep="\t")

    # Filter to reference nperm so panels retain original semantics when TSV
    # contains multiple nperm values from the permutation benchmark grid.
    REF_NPERM = 500
    if "nperm" in df.columns:
        df = df[df["nperm"] == REF_NPERM]

    # Compute medians across repeats
    median = df.groupby(["n_genes", "n_cells"]).median(numeric_only=True).reset_index()

    genes_vals = sorted(median["n_genes"].unique())
    cells_vals = sorted(median["n_cells"].unique())

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    # ===================== TOP ROW: TIMING =====================

    # --- Panel A: Heatmap of median total time ---
    ax = axes[0, 0]
    pivot = median.pivot(index="n_genes", columns="n_cells", values="total_time_s")
    pivot = pivot.reindex(index=genes_vals, columns=cells_vals)

    vmin = max(pivot.values.min(), 0.01)
    vmax = pivot.values.max()
    im = ax.imshow(
        pivot.values,
        norm=LogNorm(vmin=vmin, vmax=vmax),
        cmap="YlOrRd",
        aspect="auto",
    )
    ax.set_xticks(range(len(cells_vals)))
    ax.set_xticklabels(cells_vals)
    ax.set_yticks(range(len(genes_vals)))
    ax.set_yticklabels(genes_vals)
    ax.set_xlabel("n_cells")
    ax.set_ylabel("n_genes")
    ax.set_title("A) Median total time (s)")

    # Annotate cells
    for i in range(len(genes_vals)):
        for j in range(len(cells_vals)):
            val = pivot.values[i, j]
            ax.text(j, i, f"{val:.1f}", ha="center", va="center", fontsize=8,
                    color="white" if val > vmax * 0.3 else "black")

    fig.colorbar(im, ax=ax, label="seconds")

    # --- Panel B: Line plot of median time vs n_genes ---
    ax = axes[0, 1]
    for nc in cells_vals:
        sub = median[median["n_cells"] == nc].sort_values("n_genes")
        ax.plot(sub["n_genes"], sub["total_time_s"], marker="o", label=f"{nc} cells")

    ax.set_yscale("log")
    ax.set_xlabel("n_genes")
    ax.set_ylabel("Median total time (s)")
    ax.set_title("B) Scaling with gene count")
    ax.legend(title="n_cells")
    ax.grid(True, alpha=0.3)

    # --- Panel C: Stacked bars SVD vs permutation at largest n_cells ---
    ax = axes[0, 2]
    max_cells = max(cells_vals)
    sub = median[median["n_cells"] == max_cells].sort_values("n_genes")

    x = np.arange(len(sub))
    width = 0.6
    ax.bar(x, sub["svd_time_s"], width, label="SVD", color="#4C72B0")
    ax.bar(x, sub["perm_time_s"], width, bottom=sub["svd_time_s"],
           label="Permutation", color="#DD8452")

    ax.set_xticks(x)
    ax.set_xticklabels(sub["n_genes"].astype(int))
    ax.set_xlabel("n_genes")
    ax.set_ylabel("Median time (s)")
    ax.set_title(f"C) Phase breakdown (n_cells={max_cells})")
    ax.legend()
    ax.grid(True, alpha=0.3, axis="y")

    # ===================== BOTTOM ROW: MEMORY =====================

    # --- Panel D: Memory heatmap (peak RSS) ---
    ax = axes[1, 0]
    pivot_mem = median.pivot(index="n_genes", columns="n_cells", values="peak_rss_mb")
    pivot_mem = pivot_mem.reindex(index=genes_vals, columns=cells_vals)

    mem_vmin = max(pivot_mem.values.min(), 1)
    mem_vmax = pivot_mem.values.max()
    im_mem = ax.imshow(
        pivot_mem.values,
        norm=LogNorm(vmin=mem_vmin, vmax=mem_vmax),
        cmap="YlGnBu",
        aspect="auto",
    )
    ax.set_xticks(range(len(cells_vals)))
    ax.set_xticklabels(cells_vals)
    ax.set_yticks(range(len(genes_vals)))
    ax.set_yticklabels(genes_vals)
    ax.set_xlabel("n_cells")
    ax.set_ylabel("n_genes")
    ax.set_title("D) Median peak RSS (MB)")

    for i in range(len(genes_vals)):
        for j in range(len(cells_vals)):
            val = pivot_mem.values[i, j]
            ax.text(j, i, f"{val:.0f}", ha="center", va="center", fontsize=8,
                    color="white" if val > mem_vmax * 0.3 else "black")

    fig.colorbar(im_mem, ax=ax, label="MB")

    # --- Panel E: Memory scaling lines (peak RSS vs n_genes) ---
    ax = axes[1, 1]
    for nc in cells_vals:
        sub = median[median["n_cells"] == nc].sort_values("n_genes")
        ax.plot(sub["n_genes"], sub["peak_rss_mb"], marker="o", label=f"{nc} cells")

    ax.set_yscale("log")
    ax.set_xlabel("n_genes")
    ax.set_ylabel("Median peak RSS (MB)")
    ax.set_title("E) Memory scaling with gene count")
    ax.legend(title="n_cells")
    ax.grid(True, alpha=0.3)

    # --- Panel F: Phase memory breakdown (stacked bars at largest n_cells) ---
    ax = axes[1, 2]
    sub = median[median["n_cells"] == max_cells].sort_values("n_genes")

    x = np.arange(len(sub))
    width = 0.6
    ax.bar(x, sub["svd_rss_delta_mb"], width, label="SVD", color="#4C72B0")
    ax.bar(x, sub["perm_rss_delta_mb"], width, bottom=sub["svd_rss_delta_mb"],
           label="Permutation", color="#DD8452")

    ax.set_xticks(x)
    ax.set_xticklabels(sub["n_genes"].astype(int))
    ax.set_xlabel("n_genes")
    ax.set_ylabel("Median RSS delta (MB)")
    ax.set_title(f"F) Phase memory breakdown (n_cells={max_cells})")
    ax.legend()
    ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    plot_benchmark(
        summary_path=str(snakemake.input[0]),
        output_path=str(snakemake.output[0]),
    )
