# Epigenetic Clock Analysis of Methotrexate Response in Rheumatoid Arthritis

## What this does

Tests whether **epigenetic clock acceleration** — a composite "biological age" signal derived from baseline (pre-treatment) DNA methylation — predicts whether a rheumatoid arthritis (RA) patient responds to Methotrexate (MTX).

This is the Module 2 component of a larger project benchmarking nine epigenetic clocks (Horvath1, Horvath2, Hannum, PhenoAge, Zhang2019, Lin, Knight, HRS-InCh PhenoAge, VidalBralo) on RA case-control discrimination (Module 1) and, here, on MTX treatment response.

**Why it's novel:** the paper behind the original validation dataset (E-MTAB-10159; Gosselt et al. 2021) ran a standard single-CpG methylome-wide association study (MWAS) and found no genome-wide significant predictor of MTX response. This script asks a different question of the same kind of data: does treating methylation as a *composite aging signal* (an epigenetic clock), rather than testing CpGs one at a time, detect a predictive signal that single-CpG MWAS missed?

## What it needs from your data

- Raw IDAT files (EPIC or 450K array) for each patient, baseline/pre-treatment timepoint
- A metadata file (SDRF-style tab-delimited, or adapt Section 2 of the script) with, per patient:
  - age, sex
  - MTX response category (responder/non-responder, however your cohort labels it)
  - ideally: baseline disease activity score (e.g. DAS28) and delta/change in disease activity
  - ideally: cell-type proportions (Houseman/EpiDISH-derived CD4T, CD8T, B cell, NK, monocyte proportions) — if you don't have these, EpiDISH can compute them from the same IDAT files, or the relevant covariates can be dropped from the model formulas in Sections 8, 10 and 12

If your metadata file's column names differ from the SDRF convention used here, edit the right-hand-side strings in **Section 2** only — everything else in the script refers to the standardized column names it creates.

## How to run

1. Install dependencies (Bioconductor + methylCIPHER from GitHub):
   ```r
   install.packages("BiocManager")
   BiocManager::install(c("minfi", "limma",
     "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
     "IlluminaHumanMethylationEPICmanifest"))
   install.packages(c("pROC", "ggplot2", "dplyr", "glmnet", "devtools"))
   devtools::install_github("MorganLevineLab/methylCIPHER")
   ```
2. Open `scripts/mtx_clock_analysis.R` and edit **Section 0** (`IDAT_DIR`, `METADATA_FILE`).
3. If your metadata column names differ from the defaults, edit **Section 2**.
4. Run the full script (RStudio: Source, or `Rscript scripts/mtx_clock_analysis.R`). Runtime is roughly 15-30 minutes depending on cohort size, most of it in IDAT loading, normalization, and the MWAS.

## What it produces

All outputs are written to `outputs/`:

| File | Contains | Safe to share externally? |
|---|---|---|
| `internal_full_results.RData` | Per-patient clock scores, QC'd methylation matrix | **No — keep local.** This is individual-level data. |
| `summary_results_SHARE.RData` | Aggregate statistics only (group means, p-values, AUCs + CIs, counts) | **Yes** |
| `SHARE_*.csv` (7 files) | Same aggregate statistics as above, in plain CSV | **Yes** |
| `PhenoAge_by_response.png`, `HeadToHead_AUC.png`, `MWAS_Volcano.png` | Group-level plots, no per-patient identifiers | **Yes** |

The script prints this split explicitly at the end of its run, and deliberately separates individual-level outputs from aggregate ones so that only the `SHARE_*` files and the PNG plots need to be sent back — nothing that could be considered patient-level data leaves the file produced by `internal_full_results.RData`.

## Analyses performed

1. **Clock acceleration vs response** — t-tests (FDR-corrected) comparing responders vs non-responders on each clock's age acceleration
2. **Logistic regression** — does adding clock acceleration to a clinical-only model (age, sex, baseline DAS28, cell proportions) improve AUC?
3. **Continuous correlation** — Spearman correlation between clock acceleration and deltaDAS28 (degree of response)
4. **MWAS** — standard single-CpG methylome-wide association test against response, for comparison/replication
5. **DMARD Response Score (DRS)** — an elastic-net composite built from the top MWAS CpGs, cross-validated
6. **Head-to-head comparison** — AUC (with 95% CI) of DRS vs each clock vs a clinical-only baseline model

## Contact

Questions about the script or methodology: Rithvik, [your email] — Master's student, Bioinformatics and Data Science.
