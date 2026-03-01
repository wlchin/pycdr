"""Generate 3-panel benchmark scaling figure."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def plot_benchmark(summary_path, output_path):
    df = pd.read_csv(summary_path, sep="\t")

    # Compute medians across repeats
    median = df.groupby(["n_genes", "n_cells"]).median(numeric_only=True).reset_index()

    genes_vals = sorted(median["n_genes"].unique())
    cells_vals = sorted(median["n_cells"].unique())

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # --- Panel A: Heatmap of median total time ---
    ax = axes[0]
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
    ax = axes[1]
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
    ax = axes[2]
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

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    plot_benchmark(
        summary_path=str(snakemake.input[0]),
        output_path=str(snakemake.output[0]),
    )
