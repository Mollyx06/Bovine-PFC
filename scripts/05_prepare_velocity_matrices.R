#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
})

option_list <- list(
  make_option("--seurat-rds", type = "character",
              help = "Final Seurat object with cell-type annotations and UMAP."),
  make_option("--velocity-sheet", type = "character",
              help = "CSV/TSV with sample_id, library_id, spliced_mtx, unspliced_mtx, barcodes, genes."),
  make_option("--out-dir", type = "character"),
  make_option("--celltype-key", type = "character", default = "anno_sub"),
  make_option("--include-label-file", type = "character", default = NA,
              help = "Optional text file with one cell label per line to retain."),
  make_option("--include-labels", type = "character", default = NA,
              help = "Optional comma-separated labels to retain.")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

read_sheet <- function(path) {
  sep <- ifelse(grepl("\\.csv$", path, ignore.case = TRUE), ",", "\t")
  read.table(path, header = TRUE, sep = sep, stringsAsFactors = FALSE, check.names = FALSE)
}

read_labels <- function() {
  labels <- character()
  if (!is.na(opt$include_label_file) && file.exists(opt$include_label_file)) {
    labels <- c(labels, trimws(readLines(opt$include_label_file, warn = FALSE)))
  }
  if (!is.na(opt$include_labels) && nzchar(opt$include_labels)) {
    labels <- c(labels, trimws(unlist(strsplit(opt$include_labels, ",", fixed = TRUE))))
  }
  unique(labels[nzchar(labels)])
}

read_velocity_matrix <- function(path) {
  if (grepl("\\.gz$", path)) {
    as(Matrix::readMM(gzfile(path)), "dgCMatrix")
  } else {
    as(Matrix::readMM(path), "dgCMatrix")
  }
}

object <- readRDS(opt$seurat_rds)
if (!opt$celltype_key %in% colnames(object@meta.data)) {
  stop("celltype key not found in Seurat metadata: ", opt$celltype_key)
}

sheet <- read_sheet(opt$velocity_sheet)
required <- c("sample_id", "library_id", "spliced_mtx", "unspliced_mtx", "barcodes", "genes")
if (!all(required %in% colnames(sheet))) {
  stop("Velocity sheet must contain: ", paste(required, collapse = ", "))
}

emat_list <- list()
nmat_list <- list()
for (i in seq_len(nrow(sheet))) {
  row <- sheet[i, ]
  sample_id <- row[["sample_id"]]
  library_id <- row[["library_id"]]
  prefix <- paste0(sample_id, "_")
  suffix <- paste0("_", library_id)

  spliced <- read_velocity_matrix(row[["spliced_mtx"]])
  unspliced <- read_velocity_matrix(row[["unspliced_mtx"]])
  barcodes <- readLines(row[["barcodes"]])
  genes <- readLines(row[["genes"]])

  colnames(spliced) <- paste0(prefix, barcodes, suffix)
  colnames(unspliced) <- paste0(prefix, barcodes, suffix)
  rownames(spliced) <- genes
  rownames(unspliced) <- genes
  key <- paste(sample_id, library_id, sep = "_")
  emat_list[[key]] <- spliced
  nmat_list[[key]] <- unspliced
}

common_genes <- Reduce(intersect, lapply(emat_list, rownames))
emat <- do.call(cbind, lapply(emat_list, function(x) x[common_genes, , drop = FALSE]))
nmat <- do.call(cbind, lapply(nmat_list, function(x) x[common_genes, , drop = FALSE]))

common_cells <- intersect(colnames(emat), colnames(object))
emat <- emat[, common_cells, drop = FALSE]
nmat <- nmat[, common_cells, drop = FALSE]
object <- object[, common_cells]

include_labels <- read_labels()
if (length(include_labels) > 0) {
  keep_cells <- colnames(object)[object@meta.data[[opt$celltype_key]] %in% include_labels]
  emat <- emat[, keep_cells, drop = FALSE]
  nmat <- nmat[, keep_cells, drop = FALSE]
  object <- object[, keep_cells]
}

writeMM(emat, file.path(opt$out_dir, "velocity_spliced.mtx"))
writeMM(nmat, file.path(opt$out_dir, "velocity_unspliced.mtx"))
writeLines(rownames(emat), file.path(opt$out_dir, "velocity_genes.txt"))
writeLines(colnames(emat), file.path(opt$out_dir, "velocity_barcodes.txt"))
write.csv(Embeddings(object, "umap"), file.path(opt$out_dir, "velocity_umap.csv"))
metadata <- object@meta.data[, opt$celltype_key, drop = FALSE]
colnames(metadata) <- "cell_type"
write.csv(metadata, file.path(opt$out_dir, "velocity_cell_types.csv"))
saveRDS(object, file.path(opt$out_dir, "velocity_reference.rds"))
