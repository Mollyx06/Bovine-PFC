#!/usr/bin/env python
"""Import filtered snRNA-seq matrices into per-library AnnData files."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import anndata as ad
import pandas as pd
from scipy.io import mmread


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sample-sheet",
        required=True,
        help="CSV/TSV with sample_id, library_id, and matrix_dir columns.",
    )
    parser.add_argument("--out-dir", required=True)
    parser.add_argument(
        "--delimiter",
        default=None,
        help="Optional sample-sheet delimiter. Defaults by file suffix.",
    )
    return parser.parse_args()


def read_table(path: Path, delimiter: str | None = None) -> pd.DataFrame:
    sep = delimiter if delimiter is not None else ("," if path.suffix.lower() == ".csv" else "\t")
    return pd.read_csv(path, sep=sep)


def first_existing(directory: Path, names: list[str]) -> Path:
    for name in names:
        path = directory / name
        if path.exists():
            return path
    raise FileNotFoundError(f"None of {names} found in {directory}")


def read_lines(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return [line.rstrip("\n").split("\t")[0] for line in handle]


def import_matrix(matrix_dir: Path, sample_id: str, library_id: str) -> ad.AnnData:
    matrix_path = first_existing(matrix_dir, ["matrix.mtx", "matrix.mtx.gz"])
    barcode_path = first_existing(matrix_dir, ["barcodes.tsv", "barcodes.tsv.gz"])
    feature_path = first_existing(
        matrix_dir,
        ["features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"],
    )

    counts = mmread(str(matrix_path)).T.tocsr()
    barcodes = [f"{sample_id}_{barcode}_{library_id}" for barcode in read_lines(barcode_path)]
    genes = read_lines(feature_path)

    adata = ad.AnnData(counts)
    adata.obs_names = barcodes
    adata.var_names = genes
    adata.var_names_make_unique()
    adata.obs["sample_id"] = sample_id
    adata.obs["library_id"] = library_id
    return adata


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sheet = read_table(Path(args.sample_sheet), args.delimiter)

    required = {"sample_id", "library_id", "matrix_dir"}
    missing = required - set(sheet.columns)
    if missing:
        raise ValueError(f"Missing sample-sheet columns: {sorted(missing)}")

    records = []
    for row in sheet.itertuples(index=False):
        sample_id = str(row.sample_id)
        library_id = str(row.library_id)
        adata = import_matrix(Path(row.matrix_dir), sample_id, library_id)
        out_path = out_dir / f"{sample_id}_{library_id}.snrna_raw.h5ad"
        adata.write_h5ad(out_path)
        records.append(
            {
                "sample_id": sample_id,
                "library_id": library_id,
                "n_cells": adata.n_obs,
                "n_genes": adata.n_vars,
                "file": str(out_path),
            }
        )

    pd.DataFrame(records).to_csv(out_dir / "snrna_raw_manifest.csv", index=False)


if __name__ == "__main__":
    main()
