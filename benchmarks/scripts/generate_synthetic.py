"""Generate synthetic h5ad files for benchmarking."""

import anndata as ad
import numpy as np
import pandas as pd


def generate_synthetic(n_genes, n_cells, n_phenotypes, seed, output_path):
    rng = np.random.default_rng(seed)

    # Gene-specific Poisson base rates (log-space) for realistic variance structure
    gene_log_means = np.log(rng.uniform(1, 20, size=n_genes))

    # Latent gene modules: co-expression structure like real scRNA-seq.
    # 30 latent factors keep SVD factor count at ~25-30 (matching typical HVG data)
    # and bound permutation-phase memory to ~1.3 GB for the largest grid point.
    n_latent = 30
    cell_factors = rng.normal(0, 1, size=(n_cells, n_latent))
    gene_loadings = rng.normal(0, 1, size=(n_latent, n_genes))

    # Balanced phenotype groups
    labels = np.repeat(
        [f"condition_{i}" for i in range(n_phenotypes)],
        n_cells // n_phenotypes,
    )
    if len(labels) < n_cells:
        labels = np.append(labels, ["condition_0"] * (n_cells - len(labels)))

    # Condition-specific perturbation of latent factors
    for c in range(n_phenotypes):
        mask = labels == f"condition_{c}"
        shift = rng.normal(0, 0.3, size=n_latent)
        cell_factors[mask] += shift

    # log(rate) = gene_mean + latent_signal * 0.3
    # Moderate signal strength gives realistic gene-gene correlations
    log_rates = gene_log_means + (cell_factors @ gene_loadings) * 0.3
    rates = np.exp(np.clip(log_rates, -2, 6))

    # Poisson count data + log1p-transform (CDR expects log-transformed data)
    counts = rng.poisson(lam=rates).astype(np.float32)
    X = np.log1p(counts)

    obs = pd.DataFrame({"condition": labels}, index=[f"cell_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.write_h5ad(output_path)


if __name__ == "__main__":
    # Snakemake script interface
    generate_synthetic(
        n_genes=int(snakemake.params.n_genes),
        n_cells=int(snakemake.params.n_cells),
        n_phenotypes=int(snakemake.params.n_phenotypes),
        seed=int(snakemake.params.seed),
        output_path=str(snakemake.output[0]),
    )
