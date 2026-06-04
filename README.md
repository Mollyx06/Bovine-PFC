# Bovine-PFC

Reproducible analysis workflow for the manuscript-associated bovine prefrontal
cortex study.

This repository is intentionally limited to portable code, configuration
templates, and environment notes. It does not track raw data, processed data,
figures, manuscript files, reference PDFs, notebooks, or generated results.

## Repository layout

- `scripts/`: one ordered, manuscript-facing analysis workflow.
- `config/`: example input sheets and private-configuration templates.
- `env/`: Python and R environment specifications.
- `docs/`: input schema and workflow notes.

## Quick start

1. Create the Python environment:

   ```bash
   conda env create -f env/environment.yml
   conda activate bovine-pfc
   ```

2. Install R dependencies listed in `env/R_packages.R`.

3. Copy the example files in `config/` to private local files outside version
   control, then replace placeholder paths and labels with manuscript-specific
   values.

4. Run the scripts in numerical order, using untracked directories such as
   `data/`, `processed/`, and `results/` for all inputs and outputs.

## Data policy

Large or sensitive files must remain outside Git. The `.gitignore` file blocks
common single-cell, spatial-transcriptomics, image, document, and result-file
formats by default. If a reviewer needs data access, provide it through the
approved manuscript data repository or controlled-access mechanism, then point
the local config files to those downloaded files.

## Workflow summary

The public workflow covers:

1. importing filtered snRNA-seq count matrices;
2. batch correction, clustering, and subtype annotation with external maps;
3. spatial-domain analysis from AnnData inputs;
4. label transfer from snRNA-seq reference to spatial bins;
5. RNA-velocity matrix preparation and scVelo analysis;
6. projection of velocity-derived scores to spatial data;
7. target-program scoring and neighborhood summaries;
8. CellPhoneDB interaction summaries;
9. optional marker GO enrichment.

All sample identifiers, annotation maps, target labels, gene lists, interaction
lists, and file paths are supplied through local config files rather than being
hard-coded in the scripts.
