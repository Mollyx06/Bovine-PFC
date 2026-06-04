#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(Matrix)
})

option_list <- list(
  make_option("--sample-sheet", type = "character",
              help = "CSV/TSV with sample_id, library_id, and matrix_dir or rds_path."),
  make_option("--major-annotation-map", type = "character",
              help = "CSV/TSV with sample_id, cluster, annotation."),
  make_option("--recluster-map", type = "character", default = NA,
              help = "CSV/TSV with subset_id, include_annotations, and optional resolution."),
  make_option("--subtype-annotation-map", type = "character", default = NA,
              help = "CSV/TSV with subset_id, cluster, subtype."),
  make_option("--out-dir", type = "character"),
  make_option("--dims-harmony", type = "integer", default = 30),
  make_option("--major-resolution", type = "double", default = 0.2),
  make_option("--sub-resolution", type = "double", default = 0.3),
  make_option("--min-features", type = "integer", default = 1000),
  make_option("--max-percent-mt", type = "double", default = 5),
  make_option("--mt-pattern", type = "character", default = "^MT-"),
  make_option("--batch-vars", type = "character", default = "library_id"),
  make_option("--merge-batch-vars", type = "character", default = "sample_id,library_id"),
  make_option("--exclude-final-labels", type = "character", default = "Other,Unclassified"),
  make_option("--export-h5ad", action = "store_true", default = FALSE)
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

read_sheet <- function(path) {
  sep <- ifelse(grepl("\\.csv$", path, ignore.case = TRUE), ",", "\t")
  read.table(path, header = TRUE, sep = sep, stringsAsFactors = FALSE, check.names = FALSE)
}

split_values <- function(x, sep = ",") {
  if (is.na(x) || !nzchar(x)) {
    return(character())
  }
  trimws(unlist(strsplit(x, sep, fixed = TRUE)))
}

split_semicolon <- function(x) {
  trimws(unlist(strsplit(x, "\\s*;\\s*")))
}

first_existing <- function(directory, names) {
  for (name in names) {
    path <- file.path(directory, name)
    if (file.exists(path)) {
      return(path)
    }
  }
  stop("None of the expected files were found in: ", directory)
}

ratio_filter_markers <- function(markers, ratio_cutoff = 0.5) {
  if (!all(c("pct.1", "pct.2") %in% colnames(markers))) {
    return(markers)
  }
  markers$ratio <- (markers$pct.1 - markers$pct.2) / pmax(markers$pct.1, .Machine$double.eps)
  markers[markers$ratio > ratio_cutoff, , drop = FALSE]
}

use_harmony <- function(obj, vars, dims) {
  vars <- intersect(vars, colnames(obj@meta.data))
  vars <- vars[vapply(vars, function(v) length(unique(obj@meta.data[[v]])) > 1, logical(1))]
  if (length(vars) == 0) {
    return(list(object = obj, reduction = "pca"))
  }
  obj <- RunHarmony(obj, group.by.vars = vars)
  list(object = obj, reduction = "harmony")
}

load_or_merge_sample <- function(sample_id, sheet) {
  rows <- sheet[sheet$sample_id == sample_id, , drop = FALSE]
  objs <- lapply(seq_len(nrow(rows)), function(i) {
    lib_id <- as.character(rows[i, "library_id"])
    has_rds <- "rds_path" %in% colnames(rows) &&
      !is.na(rows[i, "rds_path"]) && nzchar(rows[i, "rds_path"])
    has_matrix <- "matrix_dir" %in% colnames(rows) &&
      !is.na(rows[i, "matrix_dir"]) && nzchar(rows[i, "matrix_dir"])

    if (has_rds) {
      obj <- readRDS(rows[i, "rds_path"])
    } else if (has_matrix) {
      matrix_dir <- rows[i, "matrix_dir"]
      counts <- ReadMtx(
        mtx = first_existing(matrix_dir, c("matrix.mtx.gz", "matrix.mtx")),
        cells = first_existing(matrix_dir, c("barcodes.tsv.gz", "barcodes.tsv")),
        features = first_existing(matrix_dir, c("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv")),
        feature.column = 1
      )
      obj <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = 3, min.features = 200)
    } else {
      stop("Rows must contain either rds_path or matrix_dir.")
    }
    obj$sample_id <- sample_id
    obj$library_id <- lib_id
    obj
  })

  if (length(objs) == 1) {
    merged <- objs[[1]]
  } else {
    merged <- merge(objs[[1]], y = objs[-1], add.cell.ids = paste0(sample_id, "_", seq_along(objs)))
  }
  merged$sample_id <- sample_id
  merged
}

process_seurat <- function(obj, resolution, dims, batch_vars) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = opt$mt_pattern)
  obj <- subset(obj, subset = nFeature_RNA > opt$min_features & percent.mt < opt$max_percent_mt)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, features = VariableFeatures(obj))
  hm <- use_harmony(obj, batch_vars, dims)
  obj <- hm$object
  obj <- FindNeighbors(obj, reduction = hm$reduction, dims = 1:dims)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj, reduction = hm$reduction, dims = 1:dims)
  obj
}

annotate_major <- function(obj, sample_id, major_map) {
  map_rows <- major_map[major_map$sample_id %in% c(sample_id, "*"), , drop = FALSE]
  if (nrow(map_rows) == 0) {
    stop("No major annotation rows found for sample_id: ", sample_id)
  }
  cluster_map <- setNames(as.character(map_rows$annotation), as.character(map_rows$cluster))
  obj$annotation <- unname(cluster_map[as.character(obj$seurat_clusters)])
  missing <- sort(unique(as.character(obj$seurat_clusters)[is.na(obj$annotation)]))
  if (length(missing) > 0) {
    stop("Missing major annotations for sample ", sample_id, " clusters: ", paste(missing, collapse = ","))
  }
  obj
}

recluster_subset <- function(obj, subset_row, subtype_map) {
  subset_id <- as.character(subset_row[["subset_id"]])
  keep_annotations <- split_semicolon(subset_row[["include_annotations"]])
  resolution <- opt$sub_resolution
  if ("resolution" %in% names(subset_row) && !is.na(subset_row[["resolution"]])) {
    resolution <- as.numeric(subset_row[["resolution"]])
  }

  sub <- subset(obj, subset = annotation %in% keep_annotations)
  sub <- NormalizeData(sub)
  sub <- FindVariableFeatures(sub)
  sub <- ScaleData(sub)
  sub <- RunPCA(sub, npcs = 30)
  sub <- FindNeighbors(sub, dims = 1:15)
  sub <- FindClusters(sub, resolution = resolution)
  sub <- RunUMAP(sub, dims = 1:15)
  Idents(sub) <- sub$seurat_clusters

  markers <- FindAllMarkers(sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(markers, file.path(opt$out_dir, paste0(subset_id, "_markers.csv")), row.names = FALSE)
  write.csv(ratio_filter_markers(markers),
            file.path(opt$out_dir, paste0(subset_id, "_markers_ratio0.5.csv")),
            row.names = FALSE)

  map_rows <- subtype_map[subtype_map$subset_id == subset_id, , drop = FALSE]
  if (nrow(map_rows) > 0) {
    cluster_map <- setNames(as.character(map_rows$subtype), as.character(map_rows$cluster))
    sub$subtype <- unname(cluster_map[as.character(sub$seurat_clusters)])
    missing <- sort(unique(as.character(sub$seurat_clusters)[is.na(sub$subtype)]))
    if (length(missing) > 0) {
      stop("Missing subtype annotations for subset ", subset_id, " clusters: ", paste(missing, collapse = ","))
    }
  } else {
    sub$subtype <- paste0(subset_id, "_cluster_", sub$seurat_clusters)
  }
  sub
}

sheet <- read_sheet(opt$sample_sheet)
major_map <- read_sheet(opt$major_annotation_map)
required_major <- c("sample_id", "cluster", "annotation")
if (!all(required_major %in% colnames(major_map))) {
  stop("Major annotation map must contain: ", paste(required_major, collapse = ", "))
}

sample_objs <- list()
for (sample_id in unique(sheet$sample_id)) {
  obj <- load_or_merge_sample(sample_id, sheet)
  obj <- process_seurat(obj, opt$major_resolution, opt$dims_harmony, split_values(opt$batch_vars))
  obj <- annotate_major(obj, sample_id, major_map)
  saveRDS(obj, file.path(opt$out_dir, paste0(sample_id, ".snrna_harmony_major.rds")))
  sample_objs[[sample_id]] <- obj
}

if (length(sample_objs) == 1) {
  merged <- sample_objs[[1]]
} else {
  merged <- merge(sample_objs[[1]], y = sample_objs[-1], add.cell.ids = names(sample_objs))
}

merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
merged <- RunPCA(merged, features = VariableFeatures(merged))
hm <- use_harmony(merged, split_values(opt$merge_batch_vars), 20)
merged <- hm$object
merged <- RunUMAP(merged, reduction = hm$reduction, dims = 1:20)
saveRDS(merged, file.path(opt$out_dir, "snrna_harmony_major_merged.rds"))

merged$anno_sub <- as.character(merged$annotation)
if (!is.na(opt$recluster_map) && file.exists(opt$recluster_map)) {
  recluster_map <- read_sheet(opt$recluster_map)
  subtype_map <- data.frame(subset_id = character(), cluster = character(), subtype = character())
  if (!is.na(opt$subtype_annotation_map) && file.exists(opt$subtype_annotation_map)) {
    subtype_map <- read_sheet(opt$subtype_annotation_map)
  }
  for (i in seq_len(nrow(recluster_map))) {
    sub <- recluster_subset(merged, recluster_map[i, ], subtype_map)
    saveRDS(sub, file.path(opt$out_dir, paste0(recluster_map[i, "subset_id"], "_subtyped.rds")))
    subtype_values <- sub$subtype
    names(subtype_values) <- rownames(sub@meta.data)
    common <- intersect(rownames(merged@meta.data), names(subtype_values))
    merged$anno_sub[common] <- subtype_values[common]
  }
}

exclude <- split_values(opt$exclude_final_labels)
if (length(exclude) > 0) {
  merged <- subset(merged, subset = !anno_sub %in% exclude)
}

Idents(merged) <- merged$anno_sub
sub_markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(sub_markers, file.path(opt$out_dir, "snrna_subtype_markers.csv"), row.names = FALSE)
write.csv(ratio_filter_markers(sub_markers),
          file.path(opt$out_dir, "snrna_subtype_markers_ratio0.5.csv"),
          row.names = FALSE)
saveRDS(merged, file.path(opt$out_dir, "snrna_subtype_reference.rds"))

if (opt$export_h5ad) {
  if (!requireNamespace("sceasy", quietly = TRUE)) {
    warning("sceasy is not installed; skipping AnnData export.")
  } else {
    sceasy::convertFormat(
      file.path(opt$out_dir, "snrna_subtype_reference.rds"),
      from = "seurat",
      to = "anndata",
      outFile = file.path(opt$out_dir, "snrna_subtype_reference.h5ad")
    )
  }
}
