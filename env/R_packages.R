required_cran <- c(
  "dplyr",
  "harmony",
  "Matrix",
  "optparse",
  "Seurat"
)

required_bioc <- c(
  "clusterProfiler",
  "org.Hs.eg.db"
)

missing_cran <- setdiff(required_cran, rownames(installed.packages()))
if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

missing_bioc <- setdiff(required_bioc, rownames(installed.packages()))
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

message("Optional package for Seurat-to-AnnData export: sceasy")
