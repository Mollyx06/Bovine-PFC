#!/usr/bin/env python
"""Spatial clustering, optional domain annotation, markers, and Moran's I."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import issparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--section-sheet", required=True, help="CSV/TSV with section_id and h5ad_path.")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--domain-map", default=None, help="Optional CSV/TSV with cluster and annotation columns.")
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--cluster-resolution", type=float, default=1.0)
    return parser.parse_args()


def read_sheet(path: Path) -> pd.DataFrame:
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    return pd.read_csv(path, sep=sep)


def load_domain_map(path: str | None) -> dict[str, str]:
    if path is None:
        return {}
    table = read_sheet(Path(path))
    required = {"cluster", "annotation"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"Missing domain-map columns: {sorted(missing)}")
    return dict(zip(table["cluster"].astype(str), table["annotation"].astype(str)))


def main() -> None:
    import scanpy as sc

    try:
        import stereo as st
    except ImportError as exc:
        raise ImportError("stereopy is required for this analysis.") from exc

    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sheet = read_sheet(Path(args.section_sheet))
    domain_map = load_domain_map(args.domain_map)

    all_markers = []
    moran_records = []

    for row in sheet.itertuples(index=False):
        section_id = str(row.section_id)
        section_dir = out_dir / section_id
        section_dir.mkdir(parents=True, exist_ok=True)

        data = st.io.read_h5ad(file_path=row.h5ad_path, flavor="scanpy", spatial_key="spatial")
        data.tl.cal_qc()
        data.tl.normalize_total()
        data.tl.log1p()
        data.tl.highly_variable_genes(
            min_mean=0.0125,
            max_mean=3,
            min_disp=0.5,
            n_top_genes=args.n_top_genes,
            res_key="highly_variable_genes",
        )
        data.tl.pca(use_highly_genes=True, n_pcs=args.n_pcs, res_key="pca")
        data.tl.neighbors(pca_res_key="pca", n_pcs=args.n_pcs, res_key="neighbors")
        data.tl.spatial_neighbors(res_key="spatial_neighbors")
        data.tl.umap(pca_res_key="pca", neighbors_res_key="neighbors", res_key="umap")
        data.tl.leiden(neighbors_res_key="neighbors", res_key="leiden", resolution=args.cluster_resolution)
        data.tl.leiden(neighbors_res_key="spatial_neighbors", res_key="spatial_leiden", resolution=args.cluster_resolution)

        adata = st.io.stereo_to_anndata(
            data,
            flavor="seurat",
            output=str(section_dir / f"{section_id}.spatial_domain.h5ad"),
        )
        if "spatial_leiden" in adata.obs:
            adata.obs["spatial_domain"] = adata.obs["spatial_leiden"].astype(str)
            if domain_map:
                adata.obs["spatial_domain"] = adata.obs["spatial_domain"].map(domain_map).fillna("Unassigned")
            adata.write_h5ad(section_dir / f"{section_id}.spatial_annotated.h5ad")

        data.tl.find_marker_genes(
            cluster_res_key="spatial_leiden",
            method="t_test",
            use_highly_genes=False,
            use_raw=False,
            res_key="marker_genes",
        )
        for cluster_id, df in data.tl.result["marker_genes"].items():
            marker_df = df.copy()
            marker_df["cluster"] = cluster_id
            marker_df["section_id"] = section_id
            all_markers.append(marker_df)

        ad_moran = sc.read_h5ad(section_dir / f"{section_id}.spatial_domain.h5ad")
        if "spatial_connectivities" in ad_moran.obsp:
            try:
                from esda.moran import Moran
                from libpysal.weights import WSP

                weights = WSP(ad_moran.obsp["spatial_connectivities"]).to_W()
                hv = ad_moran.var.get("highly_variable", pd.Series(True, index=ad_moran.var_names))
                genes = ad_moran.var_names[hv].to_list()
                records = []
                for gene in genes:
                    vec = ad_moran[:, gene].X
                    if issparse(vec):
                        vec = vec.toarray()
                    mi = Moran(np.ravel(vec), weights)
                    records.append({"section_id": section_id, "gene": gene, "moran_i": mi.I, "p_norm": mi.p_norm})
                moran_records.append(pd.DataFrame(records))
            except ImportError:
                pass

    if all_markers:
        pd.concat(all_markers, ignore_index=True).to_csv(out_dir / "spatial_marker_table.csv", index=False)
    if moran_records:
        pd.concat(moran_records, ignore_index=True).to_csv(out_dir / "spatial_moran_i.csv", index=False)


if __name__ == "__main__":
    main()
