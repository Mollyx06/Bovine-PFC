#!/usr/bin/env python
"""Target-cell program scoring and spatial neighborhood summaries."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial import cKDTree


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-h5ad", required=True)
    parser.add_argument("--predicted-section-sheet", required=True,
                        help="CSV/TSV with section_id and h5ad_path for label-transfer spatial objects.")
    parser.add_argument("--annotated-section-sheet", default=None,
                        help="Optional CSV/TSV with section_id and h5ad_path containing spatial_domain.")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--target-label", required=True)
    parser.add_argument("--celltype-key", default="cell_type")
    parser.add_argument("--spatial-label-key", default="pred_celltype")
    parser.add_argument("--lineage-map", required=True,
                        help="CSV/TSV with label and lineage columns.")
    parser.add_argument("--fallback-gene-list", default=None,
                        help="Optional one-gene-per-line list used if DE yields too few target genes.")
    parser.add_argument("--neighbor-k", type=int, default=12)
    return parser.parse_args()


def read_sheet(path: Path | None) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame()
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    return pd.read_csv(path, sep=sep)


def read_gene_list(path: str | None) -> list[str]:
    if path is None:
        return []
    with open(path, "rt", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip() and not line.startswith("#")]


def load_lineage_map(path: str) -> dict[str, str]:
    table = read_sheet(Path(path))
    required = {"label", "lineage"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"Missing lineage-map columns: {sorted(missing)}")
    return dict(zip(table["label"].astype(str), table["lineage"].astype(str)))


def define_target_core(
    adata: sc.AnnData,
    out_dir: Path,
    celltype_key: str,
    target_label: str,
    lineage_map: dict[str, str],
    fallback_genes: list[str],
) -> list[str]:
    adata.obs["lineage_group"] = adata.obs[celltype_key].astype(str).map(lineage_map).fillna("Other").astype("category")
    target_lineage = lineage_map.get(target_label, target_label)
    compare_groups = [x for x in adata.obs["lineage_group"].cat.categories if x not in {target_lineage, "Other"}]
    adata_deg = adata[adata.obs["lineage_group"].isin(compare_groups + [target_lineage])].copy()
    all_deg = []
    for ref in compare_groups:
        key = f"target_vs_{ref}"
        sc.tl.rank_genes_groups(
            adata_deg,
            groupby="lineage_group",
            groups=[target_lineage],
            reference=ref,
            method="wilcoxon",
            key_added=key,
        )
        df = sc.get.rank_genes_groups_df(adata_deg, group=target_lineage, key=key).copy()
        df["comparison"] = key
        all_deg.append(df)

    if all_deg:
        deg_df = pd.concat(all_deg, ignore_index=True)
        deg_df.to_csv(out_dir / "target_pairwise_deg.csv", index=False)
        sig = deg_df[(deg_df["pvals_adj"] < 0.05) & (deg_df["logfoldchanges"] > 0.25)].copy()
        shared = set(sig["names"])
        for ref in compare_groups:
            shared &= set(sig.loc[sig["comparison"] == f"target_vs_{ref}", "names"])
        shared_df = (
            sig[sig["names"].isin(shared)]
            .groupby("names", as_index=False)["logfoldchanges"]
            .mean()
            .sort_values("logfoldchanges", ascending=False)
        )
        core = shared_df["names"].head(20).to_list()
    else:
        core = []

    if len(core) < 8:
        for gene in fallback_genes:
            if gene in adata.var_names and gene not in core:
                core.append(gene)
            if len(core) >= 20:
                break

    pd.DataFrame({"core_genes": core}).to_csv(out_dir / "target_core_genes.csv", index=False)
    return core


def run_enrichment(core_genes: list[str], out_dir: Path) -> None:
    try:
        import gseapy as gp
    except ImportError:
        return
    if len(core_genes) < 3:
        return
    for name, gene_set in [
        ("GO_BP", "GO_Biological_Process_2023"),
        ("Hallmark", "MSigDB_Hallmark_2020"),
    ]:
        res = gp.enrichr(gene_list=core_genes, gene_sets=gene_set, organism="Human", outdir=None)
        res.results.sort_values("Adjusted P-value").to_csv(out_dir / f"target_{name}_enrichment.csv", index=False)


def collapse_label(label: str, lineage_map: dict[str, str]) -> str:
    return lineage_map.get(str(label), "Other")


def summarize_neighbors(
    adata: sc.AnnData,
    label_key: str,
    target_label: str,
    lineage_map: dict[str, str],
    k: int,
) -> pd.DataFrame:
    coords = adata.obsm["spatial"]
    labels = adata.obs[label_key].astype(str).to_numpy()
    major = np.array([collapse_label(label, lineage_map) for label in labels])
    target_mask = labels == target_label
    if target_mask.sum() == 0:
        return pd.DataFrame()

    tree = cKDTree(coords)
    _, idx = tree.query(coords, k=k + 1)
    idx = idx[:, 1:]
    bg = pd.Series(major).value_counts(normalize=True)
    records = []
    for i in np.where(target_mask)[0]:
        records.extend(major[idx[i]].tolist())
    neigh = pd.Series(records).value_counts(normalize=True)
    out = pd.DataFrame(
        {
            "neighbor_group": neigh.index,
            "neighbor_freq": neigh.values,
            "bg_freq": [bg.get(x, 1e-9) for x in neigh.index],
        }
    )
    out["log2_enrichment"] = np.log2((out["neighbor_freq"] + 1e-9) / (out["bg_freq"] + 1e-9))
    return out.sort_values("log2_enrichment", ascending=False)


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    lineage_map = load_lineage_map(args.lineage_map)
    fallback_genes = read_gene_list(args.fallback_gene_list)

    adata = sc.read_h5ad(args.reference_h5ad)
    if args.celltype_key not in adata.obs:
        raise KeyError(f"{args.celltype_key} not found in reference object")

    core = define_target_core(
        adata,
        out_dir,
        args.celltype_key,
        args.target_label,
        lineage_map,
        fallback_genes,
    )
    run_enrichment(core, out_dir)

    score_genes = [gene for gene in core if gene in adata.var_names]
    if score_genes:
        sc.tl.score_genes(adata, gene_list=score_genes, score_name="target_core_score", use_raw=False)
        adata.obs[[args.celltype_key, "lineage_group", "velocity_pseudotime", "target_core_score"]].to_csv(
            out_dir / "target_core_score_by_cell.csv"
        )

    pred_sheet = read_sheet(Path(args.predicted_section_sheet))
    anno_sheet = read_sheet(Path(args.annotated_section_sheet)) if args.annotated_section_sheet else pd.DataFrame()
    anno_map = dict(zip(anno_sheet.get("section_id", []), anno_sheet.get("h5ad_path", [])))

    for row in pred_sheet.itertuples(index=False):
        section_id = str(row.section_id)
        sp = sc.read_h5ad(row.h5ad_path)
        genes = [gene for gene in score_genes if gene in sp.var_names]
        if genes:
            sc.tl.score_genes(sp, gene_list=genes, score_name="target_core_score", use_raw=False)
        if section_id in anno_map:
            anno = sc.read_h5ad(anno_map[section_id])
            if "spatial_domain" in anno.obs:
                common = sp.obs_names.intersection(anno.obs_names)
                sp.obs.loc[common, "spatial_domain"] = anno.obs.loc[common, "spatial_domain"].astype(str)

        sp.write_h5ad(out_dir / f"{section_id}.target_spatial_score.h5ad")
        if "target_core_score" in sp.obs:
            columns = ["target_core_score"] + (["spatial_domain"] if "spatial_domain" in sp.obs else [])
            sp.obs[columns].to_csv(out_dir / f"{section_id}.target_spatial_score.csv")
        neigh = summarize_neighbors(sp, args.spatial_label_key, args.target_label, lineage_map, args.neighbor_k)
        neigh.to_csv(out_dir / f"{section_id}.target_neighbor_enrichment.csv", index=False)


if __name__ == "__main__":
    main()
