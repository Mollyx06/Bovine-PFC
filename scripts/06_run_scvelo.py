#!/usr/bin/env python
"""Construct a scVelo object and infer velocity pseudotime."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
from scipy.io import mmread


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--velocity-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--mode", default="stochastic", choices=["stochastic", "deterministic", "dynamical"])
    return parser.parse_args()


def main() -> None:
    import scvelo as scv

    args = parse_args()
    velocity_dir = Path(args.velocity_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    spliced = mmread(velocity_dir / "velocity_spliced.mtx").T.tocsr()
    unspliced = mmread(velocity_dir / "velocity_unspliced.mtx").T.tocsr()
    genes = pd.read_csv(velocity_dir / "velocity_genes.txt", header=None)[0].astype(str).to_list()
    barcodes = pd.read_csv(velocity_dir / "velocity_barcodes.txt", header=None)[0].astype(str).to_list()
    umap = pd.read_csv(velocity_dir / "velocity_umap.csv", index_col=0)
    cell_types = pd.read_csv(velocity_dir / "velocity_cell_types.csv", index_col=0)

    adata = ad.AnnData(spliced)
    adata.layers["spliced"] = spliced
    adata.layers["unspliced"] = unspliced
    adata.obs_names = barcodes
    adata.var_names = genes
    adata.obs["cell_type"] = cell_types.loc[adata.obs_names, "cell_type"].astype("category")
    adata.obsm["X_umap"] = umap.loc[adata.obs_names].to_numpy()

    scv.pp.filter_and_normalize(adata)
    scv.pp.moments(adata)
    scv.tl.velocity(adata, mode=args.mode)
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_pseudotime(adata)

    adata.write_h5ad(out_dir / "velocity_reference.h5ad")


if __name__ == "__main__":
    main()
