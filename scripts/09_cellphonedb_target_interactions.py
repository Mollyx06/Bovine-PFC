#!/usr/bin/env python
"""Summarize CellPhoneDB interactions between a target label and partner labels."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from scipy import sparse

mpl.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "figure.dpi": 300, "savefig.dpi": 300})


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--significant-means", required=True)
    parser.add_argument("--adata", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--target-label", required=True)
    parser.add_argument("--partner-labels", required=True,
                        help="Text file with one partner label per line.")
    parser.add_argument("--celltype-key", default="anno_sub")
    parser.add_argument("--gene-list", default=None,
                        help="Optional one-gene-per-line list to use for expression violin plots.")
    parser.add_argument("--curated-pairs", default=None,
                        help="Optional one-interacting_pair-per-line order for selected pairs.")
    parser.add_argument("--top-n", type=int, default=10)
    return parser.parse_args()


def read_list(path: str | None) -> list[str]:
    if path is None:
        return []
    with open(path, "rt", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip() and not line.startswith("#")]


def cpdb_pair_columns(columns: pd.Index, target_label: str, partner_labels: list[str]) -> list[str]:
    keep = []
    partners = set(partner_labels)
    for col in columns:
        if "|" not in col:
            continue
        sender, receiver = col.split("|", 1)
        if (sender == target_label and receiver in partners) or (sender in partners and receiver == target_label):
            keep.append(col)
    return keep


def select_top(df: pd.DataFrame, pair_cols: list[str], top_n: int, curated_pairs: list[str]) -> pd.DataFrame:
    required = ["id_cp_interaction", "interacting_pair", "partner_a", "partner_b"]
    missing = [x for x in required if x not in df.columns]
    if missing:
        raise ValueError(f"Missing CellPhoneDB columns: {missing}")
    work = df[required + pair_cols].copy()
    work[pair_cols] = work[pair_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    work["total_significant"] = (work[pair_cols] > 0).sum(axis=1)
    work["mean_strength"] = work[pair_cols].replace(0, np.nan).mean(axis=1).fillna(0)
    if curated_pairs:
        present = [pair for pair in curated_pairs if pair in set(work["interacting_pair"])]
        work = work[work["interacting_pair"].isin(present)].copy()
        work["_order"] = pd.Categorical(work["interacting_pair"], categories=present, ordered=True)
        return work.sort_values("_order").drop(columns="_order")
    return work.sort_values(["total_significant", "mean_strength"], ascending=False).head(top_n)


def pretty_pair(col: str, target_label: str) -> str:
    sender, receiver = col.split("|", 1)
    sender = "Target" if sender == target_label else sender
    receiver = "Target" if receiver == target_label else receiver
    return f"{sender} -> {receiver}"


def plot_bubble(top: pd.DataFrame, pair_cols: list[str], target_label: str, out_pdf: Path) -> None:
    data = top.set_index("interacting_pair")[pair_cols].iloc[::-1]
    values = data.to_numpy(dtype=float)
    vmax = max(float(np.nanmax(values)), 1.0)
    fig, ax = plt.subplots(figsize=(max(5.2, 0.46 * len(pair_cols) + 2.4), max(4.0, 0.34 * data.shape[0] + 1.8)))
    for i, interaction in enumerate(data.index):
        for j, pair in enumerate(data.columns):
            val = data.loc[interaction, pair]
            if val > 0:
                ax.scatter(
                    j,
                    i,
                    s=25 + val * 80,
                    c=val,
                    cmap="RdBu_r",
                    vmin=0,
                    vmax=vmax,
                    edgecolor="#333333",
                    linewidth=0.35,
                )
    ax.set_xticks(range(len(pair_cols)))
    ax.set_xticklabels([pretty_pair(x, target_label) for x in pair_cols], rotation=45, ha="right")
    ax.set_yticks(range(data.shape[0]))
    ax.set_yticklabels(data.index)
    ax.grid(color="#E6E6E6", linewidth=0.5)
    fig.colorbar(
        mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=vmax), cmap="RdBu_r"),
        ax=ax,
        label="Interaction strength",
    )
    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def expression_vector(adata: ad.AnnData, gene: str) -> np.ndarray | None:
    if gene not in adata.var_names:
        return None
    x = adata[:, gene].X
    if sparse.issparse(x):
        x = x.toarray()
    return np.asarray(x).ravel()


def plot_violins(adata_path: Path, celltype_key: str, labels: list[str], genes: list[str], out_pdf: Path) -> None:
    if not genes:
        return
    data = ad.read_h5ad(adata_path)
    present = [label for label in labels if label in set(data.obs[celltype_key].astype(str))]
    sub = data[data.obs[celltype_key].isin(present)].copy()
    ncols = min(4, max(1, len(genes)))
    nrows = int(np.ceil(len(genes) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.35 * ncols, 2.25 * nrows), squeeze=False)
    for ax, gene in zip(axes.ravel(), genes):
        expr = expression_vector(sub, gene)
        if expr is None:
            ax.axis("off")
            continue
        groups = [expr[np.asarray(sub.obs[celltype_key] == label)] for label in present]
        ax.violinplot(groups, positions=np.arange(len(groups)), showmedians=True, showextrema=False, widths=0.85)
        ax.set_title(gene)
        ax.set_xticks(np.arange(len(present)))
        ax.set_xticklabels(["Target" if i == 0 else label for i, label in enumerate(present)], rotation=45, ha="right")
    for ax in axes.ravel()[len(genes):]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def plot_network(top: pd.DataFrame, out_pdf: Path) -> None:
    edges = [tuple(x.split("_", 1)) if "_" in x else (x, x) for x in top["interacting_pair"]]
    graph = nx.DiGraph()
    graph.add_edges_from(edges)
    pos = nx.spring_layout(graph, seed=7, k=1.15)
    fig, ax = plt.subplots(figsize=(5.8, 4.6))
    nx.draw_networkx_edges(graph, pos, ax=ax, arrowstyle="-|>", arrowsize=12, edge_color="#999999")
    nx.draw_networkx_nodes(graph, pos, ax=ax, node_color="#FDB462", node_size=900, edgecolors="white")
    nx.draw_networkx_labels(graph, pos, ax=ax, font_size=7, font_weight="bold")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    partners = read_list(args.partner_labels)
    genes = read_list(args.gene_list)
    curated_pairs = read_list(args.curated_pairs)

    significant = pd.read_csv(args.significant_means, sep="\t")
    pair_cols = cpdb_pair_columns(significant.columns, args.target_label, partners)
    if not pair_cols:
        raise ValueError("No target-partner CellPhoneDB pair columns were found.")
    top = select_top(significant, pair_cols, args.top_n, curated_pairs)
    top.to_csv(out_dir / "cellphonedb_target_top_interactions.csv", index=False)
    plot_bubble(top, pair_cols, args.target_label, out_dir / "cellphonedb_target_bubble.pdf")
    plot_network(top, out_dir / "cellphonedb_target_network.pdf")
    plot_violins(
        Path(args.adata),
        args.celltype_key,
        [args.target_label] + partners,
        genes,
        out_dir / "cellphonedb_target_expression_violin.pdf",
    )


if __name__ == "__main__":
    main()
