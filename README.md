# ImmunoVerse

ImmunoVerse is a curated single-cell reference atlas of 11 T-cell states grouped into six broader families. This repository contains the R analysis scripts and derived tables and figures for the comparison across 10x library chemistries and the analysis of paired blood and cerebrospinal fluid samples.

## Source data

Download the underlying single-cell RNA-seq data from NCBI GEO: [GSE144744](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144744) and [GSE138266](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138266). The scripts write downloaded files and intermediate objects to `data/raw/` and `data/processed/`. Keep substantially more workspace available than the compressed downloads require.

## R workflow

Run scripts from the repository root with `Rscript scripts/<filename>.R`. The stage order below follows the inputs and outputs of the numbered scripts. The full two-cohort workflow has not been rerun from a clean environment for this release candidate.

| Stage | Scripts |
| --- | --- |
| Environment and GSE144744 import | `00`–`03`, then `05_tcell_subclustering_safe.R` to write the full T-cell input |
| Curated v3 reference | `15`–`20` |
| v3 transfer, concordance and differential abundance | `21`–`25b` |
| v2 transfer and cross-chemistry analyses | `26`–`28` |
| Independent GSE138266 analysis | `29`–`35` |

The archive also includes earlier exploratory scripts for provenance. The core stages use Seurat, Harmony, Matrix, and the tidyverse packages. The differential-abundance stages use speckle and limma, and the external QC stage uses DropletUtils. Install these packages with CRAN or Bioconductor as appropriate before running the relevant stage. The package-version CSV in `results/tables/curated_V3_differential_abundance/` records the versions logged by that analysis.

## Results

- `results/tables/` contains the atlas definitions, family mappings, transfer diagnostics, and differential-abundance outputs, together with earlier exploratory tables.
- `results/figures/` contains generated figures for the discovery, transfer, and abundance analyses.

Raw GEO downloads, RStudio workspaces, and serialized intermediate objects are generated locally and are excluded from version control.
