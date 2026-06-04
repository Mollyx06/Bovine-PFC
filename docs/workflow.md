# Workflow

Run scripts from the repository root after activating the Python and R
environments. The commands below use example file names; copy the templates in
`config/` to untracked local files before filling in real paths and labels.

```bash
python scripts/01_import_snrna_counts.py \
  --sample-sheet config/sample_sheet_local.csv \
  --out-dir processed/01_snrna_import

Rscript scripts/02_snrna_harmony_recluster.R \
  --sample-sheet config/sample_sheet_local.csv \
  --major-annotation-map config/major_annotation_map_local.csv \
  --recluster-map config/recluster_map_local.csv \
  --subtype-annotation-map config/subtype_annotation_map_local.csv \
  --out-dir processed/02_snrna_reference \
  --export-h5ad

python scripts/03_spatial_domain_analysis.py \
  --section-sheet config/section_sheet_local.csv \
  --domain-map config/domain_map_local.csv \
  --out-dir processed/03_spatial_domain

python scripts/04_spatial_label_transfer.py \
  --reference-h5ad processed/02_snrna_reference/snrna_subtype_reference.h5ad \
  --section-sheet config/spatial_annotated_sheet_local.csv \
  --out-dir processed/04_label_transfer

Rscript scripts/05_prepare_velocity_matrices.R \
  --seurat-rds processed/02_snrna_reference/snrna_subtype_reference.rds \
  --velocity-sheet config/velocity_sheet_local.csv \
  --include-label-file config/velocity_include_labels_local.txt \
  --out-dir processed/05_velocity_input

python scripts/06_run_scvelo.py \
  --velocity-dir processed/05_velocity_input \
  --out-dir processed/06_scvelo

python scripts/07_project_velocity_to_spatial.py \
  --velocity-h5ad processed/06_scvelo/velocity_reference.h5ad \
  --section-sheet config/spatial_predicted_sheet_local.csv \
  --target-label "<private target label>" \
  --out-dir processed/07_spatial_velocity

python scripts/08_target_program_neighborhood.py \
  --reference-h5ad processed/06_scvelo/velocity_reference.h5ad \
  --predicted-section-sheet config/spatial_predicted_sheet_local.csv \
  --target-label "<private target label>" \
  --lineage-map config/lineage_map_local.csv \
  --fallback-gene-list config/target_gene_list_local.txt \
  --out-dir processed/08_target_program

python scripts/09_cellphonedb_target_interactions.py \
  --significant-means data/cellphonedb/statistical_analysis_significant_means.txt \
  --adata processed/02_snrna_reference/snrna_subtype_reference.h5ad \
  --target-label "<private target label>" \
  --partner-labels config/partner_labels_local.txt \
  --out-dir processed/09_cellphonedb

Rscript scripts/10_marker_go_enrichment.R \
  --marker-csv processed/02_snrna_reference/snrna_subtype_markers.csv \
  --out-dir processed/10_go
```
