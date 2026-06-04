#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

option_list <- list(
  make_option("--marker-csv", type = "character",
              help = "Marker table with cluster and gene columns."),
  make_option("--out-dir", type = "character"),
  make_option("--cluster-column", type = "character", default = "cluster"),
  make_option("--gene-column", type = "character", default = "gene"),
  make_option("--padj-column", type = "character", default = "p_val_adj"),
  make_option("--logfc-column", type = "character", default = "avg_log2FC"),
  make_option("--padj-cutoff", type = "double", default = 0.05),
  make_option("--logfc-cutoff", type = "double", default = 0.25),
  make_option("--keytype", type = "character", default = "SYMBOL"),
  make_option("--ontology", type = "character", default = "ALL")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

markers <- read.csv(opt$marker_csv, stringsAsFactors = FALSE, check.names = FALSE)
clusters <- unique(markers[[opt$cluster_column]])

run_go <- function(df, cluster_id) {
  genes <- df %>%
    filter(.data[[opt$cluster_column]] == cluster_id,
           .data[[opt$padj_column]] < opt$padj_cutoff,
           .data[[opt$logfc_column]] > opt$logfc_cutoff) %>%
    pull(.data[[opt$gene_column]]) %>%
    unique()
  if (length(genes) < 5) {
    return(NULL)
  }
  enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = opt$keytype,
    ont = opt$ontology,
    pAdjustMethod = "fdr",
    qvalueCutoff = 0.1,
    readable = TRUE
  )
}

for (cluster_id in clusters) {
  res <- run_go(markers, cluster_id)
  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    next
  }
  out <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(cluster_id))
  write.csv(
    as.data.frame(res),
    file.path(opt$out_dir, paste0("GO_enrichment_", out, ".csv")),
    row.names = FALSE
  )
}
