# Input Schema

Keep filled-in versions of these files outside Git. The tracked files in
`config/` are examples only.

## sample sheet

Columns:

- `sample_id`: private local sample identifier.
- `library_id`: private local library identifier.
- `matrix_dir`: optional 10x-style filtered matrix directory.
- `rds_path`: optional prebuilt Seurat object.

Each row must provide either `matrix_dir` or `rds_path`.

## major annotation map

Columns:

- `sample_id`
- `cluster`
- `annotation`

This maps first-pass Seurat clusters to broad cell classes. It is private
because it can reveal project-specific interpretation.

## recluster map

Columns:

- `subset_id`
- `include_annotations`: semicolon-separated broad annotations to recluster.
- `resolution`: optional subset-specific clustering resolution.

## subtype annotation map

Columns:

- `subset_id`
- `cluster`
- `subtype`

This maps subset-level clusters to final labels.

## section sheet

Columns:

- `section_id`
- `h5ad_path`

Use generic local section IDs if the real names are sensitive.

## velocity sheet

Columns:

- `sample_id`
- `library_id`
- `spliced_mtx`
- `unspliced_mtx`
- `barcodes`
- `genes`

Matrix files are expected in Matrix Market format.

## lineage map

Columns:

- `label`
- `lineage`

This maps final cell labels to broader groups for target-program and
neighborhood summaries.
