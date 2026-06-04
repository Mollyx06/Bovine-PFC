#!/usr/bin/env python
"""Project velocity-derived pseudotime and target-label probability to spatial bins."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import tacco as tc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--velocity-h5ad", required=True)
    parser.add_argument("--section-sheet", required=True, help="CSV/TSV with section_id and h5ad_path.")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--celltype-key", default="cell_type")
    parser.add_argument("--pseudotime-key", default="velocity_pseudotime")
    parser.add_argument("--target-label", required=True)
    return parser.parse_args()


def read_sheet(path: Path) -> pd.DataFrame:
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    return pd.read_csv(path, sep=sep)


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ref0 = sc.read_h5ad(args.velocity_h5ad)
    if "spliced" in ref0.layers:
        ref0.X = ref0.layers["spliced"].copy()

    pt_by_type = ref0.obs.groupby(args.celltype_key, observed=True)[args.pseudotime_key].median().dropna()

    for row in read_sheet(Path(args.section_sheet)).itertuples(index=False):
        section_id = str(row.section_id)
        sp0 = sc.read_h5ad(row.h5ad_path)
        common = ref0.var_names.intersection(sp0.var_names)
        ref = ref0[:, common].copy()
        sp = sp0[:, common].copy()

        tc.tl.annotate(sp, ref, args.celltype_key, result_key="velocity_label_probability", verbose=False)
        prob = sp.obsm["velocity_label_probability"]
        common_types = [label for label in prob.columns if label in pt_by_type.index]
        sp.obs["spatial_velocity_pseudotime"] = prob[common_types].dot(pt_by_type.loc[common_types])
        sp.obs["target_label_probability"] = prob[args.target_label] if args.target_label in prob.columns else np.nan

        for key in ["spatial_domain", "spatial_leiden"]:
            if key in sp0.obs:
                sp.obs[key] = sp0.obs.loc[sp.obs_names, key]

        sp.write_h5ad(out_dir / f"{section_id}.spatial_velocity.h5ad")
        sp.obs[["spatial_velocity_pseudotime", "target_label_probability"]].to_csv(
            out_dir / f"{section_id}.spatial_velocity_scores.csv"
        )

    records = []
    for path in out_dir.glob("*.spatial_velocity.h5ad"):
        adata = sc.read_h5ad(path)
        df = adata.obs[["spatial_velocity_pseudotime", "target_label_probability"]].copy()
        df["section_id"] = path.name.split(".")[0]
        if "spatial_domain" in adata.obs:
            df["spatial_domain"] = adata.obs["spatial_domain"].astype(str)
        records.append(df)
    if records:
        stat = pd.concat(records).dropna()
        if stat.shape[0] > 2:
            rho = stat["spatial_velocity_pseudotime"].corr(stat["target_label_probability"], method="spearman")
            pd.DataFrame({"spearman_rho": [rho], "n_bins": [stat.shape[0]]}).to_csv(
                out_dir / "spatial_velocity_target_correlation.csv",
                index=False,
            )


if __name__ == "__main__":
    main()
