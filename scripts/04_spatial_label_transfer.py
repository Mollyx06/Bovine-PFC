#!/usr/bin/env python
"""Transfer snRNA-seq labels to spatial bins with TACCO."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import tacco as tc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-h5ad", required=True)
    parser.add_argument("--section-sheet", required=True, help="CSV/TSV with section_id and h5ad_path.")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--label-key", default="anno_sub")
    parser.add_argument("--result-key", default="pred_celltype")
    parser.add_argument("--use-reference-prior", action="store_true")
    return parser.parse_args()


def read_sheet(path: Path) -> pd.DataFrame:
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    return pd.read_csv(path, sep=sep)


def clean_adata(adata: sc.AnnData) -> sc.AnnData:
    if adata.raw is not None:
        adata.raw = None
    for key in list(adata.varm.keys()):
        del adata.varm[key]
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.filter_cells(adata, min_counts=1)
    adata.X = adata.X.astype(np.float32)
    return adata


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    reference = sc.read_h5ad(args.reference_h5ad)
    if args.label_key not in reference.obs:
        raise KeyError(f"{args.label_key} not found in reference.obs")
    prior = reference.obs[args.label_key].value_counts(normalize=True)

    sheet = read_sheet(Path(args.section_sheet))
    for row in sheet.itertuples(index=False):
        section_id = str(row.section_id)
        adata = sc.read_h5ad(row.h5ad_path)

        common = reference.var_names.intersection(adata.var_names)
        ref = clean_adata(reference[:, common].copy())
        sp = clean_adata(adata[:, common].copy())

        kwargs = {"annotation_prior": prior} if args.use_reference_prior else {}
        tc.tl.annotate(sp, ref, args.label_key, result_key=args.result_key, verbose=False, **kwargs)
        prob = sp.obsm[args.result_key]
        sp.obs[args.result_key] = prob.idxmax(axis=1).astype("category")

        for key in ["spatial_domain", "spatial_leiden"]:
            if key in adata.obs:
                sp.obs[key] = adata.obs.loc[sp.obs_names, key]

        out_path = out_dir / f"{section_id}.label_transfer.h5ad"
        sp.write_h5ad(out_path)
        prob.to_csv(out_dir / f"{section_id}.label_transfer_probability.csv")


if __name__ == "__main__":
    main()
