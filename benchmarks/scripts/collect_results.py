"""Aggregate per-run timing JSONs and snakemake benchmark TSVs into a summary."""

import json
import os

import pandas as pd


def collect_results(json_files, benchmark_dir, output_path):
    records = []
    for jf in json_files:
        with open(jf) as f:
            rec = json.load(f)

        # Try to find matching snakemake benchmark TSV for max_rss
        n_genes = rec["n_genes"]
        n_cells = rec["n_cells"]
        nperm = rec.get("nperm", "")
        rep = rec["rep"]
        bm_path = os.path.join(
            benchmark_dir,
            f"genes{n_genes}_cells{n_cells}_nperm{nperm}_rep{rep}.tsv",
        )
        if os.path.exists(bm_path):
            bm = pd.read_csv(bm_path, sep="\t")
            if "max_rss" in bm.columns and len(bm) > 0:
                rec["max_rss_mb"] = float(bm["max_rss"].iloc[0])

        records.append(rec)

    df = pd.DataFrame(records)
    sort_cols = ["n_genes", "n_cells", "nperm", "rep"]
    sort_cols = [c for c in sort_cols if c in df.columns]
    df.sort_values(sort_cols, inplace=True)
    df.to_csv(output_path, sep="\t", index=False)


if __name__ == "__main__":
    collect_results(
        json_files=[str(f) for f in snakemake.input.jsons],
        benchmark_dir=str(snakemake.params.benchmark_dir),
        output_path=str(snakemake.output[0]),
    )
