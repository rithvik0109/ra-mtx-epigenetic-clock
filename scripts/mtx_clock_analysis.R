# =============================================================================
# Epigenetic Clock Analysis of Methotrexate (MTX) Treatment Response in RA
# =============================================================================
# This script tests whether epigenetic clock acceleration, measured from
# baseline (pre-treatment) DNA methylation, predicts whether a rheumatoid
# arthritis patient will respond to Methotrexate.
#
# It was originally developed and validated on E-MTAB-10159 (69 treatment-
# naive RA patients, baseline PBMCs, EPIC array, MTX response at 3 months).
# It is written here to run standalone against ANY cohort with the same
# structure: an EPIC/450K IDAT dataset + a metadata file giving age, sex,
# MTX response status, and (ideally) baseline disease activity and cell-type
# proportions.
#
# NOVELTY / RATIONALE
# The original E-MTAB-10159 paper (Gosselt et al. 2021) found no genome-wide
# significant baseline methylation predictor of MTX response using standard
# single-CpG MWAS. This script tests whether the epigenetic clock framework
# (treating methylation as a composite "biological age" signal rather than
# testing CpGs one at a time) detects a predictive signal that single-CpG
# MWAS missed.
#
# DATA SHARING / OUTPUT NOTE (read before running)
# Section 14 deliberately writes TWO separate files:
#   - "internal_full_results.RData"   -> KEEP LOCAL. Contains per-sample
#        clock scores and the QC'd methylation matrix. Do not send this
#        file outside your institution.
#   - "summary_results_SHARE.RData" and "summary_results_SHARE.csv" ->
#        SAFE TO SHARE. Contains only aggregate statistics (group means,
#        p-values, AUCs with confidence intervals, counts of significant
#        CpGs). No per-patient or per-CpG raw values are included.
# =============================================================================


# =============================================================================
# SECTION 0: CONFIGURATION — EDIT THESE PATHS BEFORE RUNNING
# =============================================================================

# Folder containing the raw .idat files
IDAT_DIR <- "/path/to/your/idat_folder"

# Metadata file (SDRF-style tab-delimited file, or adapt Section 2 to your
# own metadata format if it differs — see README for the column names
# this script expects)
METADATA_FILE <- file.path(IDAT_DIR, "metadata.sdrf.txt")

# Where outputs (plots, RData, csv) get written
OUTPUT_DIR <- "outputs"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("IDAT directory:", IDAT_DIR, "\n")
cat("Metadata file:", METADATA_FILE, "\n")
cat("Output directory:", OUTPUT_DIR, "\n")


# =============================================================================
# SECTION 1: LOAD LIBRARIES
# =============================================================================

required_pkgs <- c("minfi", "limma", "methylCIPHER", "pROC", "ggplot2",
                    "dplyr", "glmnet",
                    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
                    "IlluminaHumanMethylationEPICmanifest")

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(paste0(
    "Missing required package(s): ", paste(missing_pkgs, collapse = ", "),
    "\nInstall Bioconductor packages via BiocManager::install(), and\n",
    "methylCIPHER via devtools::install_github('MorganLevineLab/methylCIPHER')."
  ))
}

library(minfi)
library(limma)
library(methylCIPHER)
library(pROC)
library(ggplot2)
library(dplyr)
library(glmnet)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)


# =============================================================================
# SECTION 2: LOAD AND CLEAN METADATA
# =============================================================================
# EXPECTED COLUMNS (SDRF-style). If your metadata file uses different column
# names, edit the right-hand-side strings below to match your file's header
# — everything downstream refers only to the standardized sdrf_grn$... names.
# =============================================================================

cat("\nLoading metadata...\n")

sdrf <- read.delim(METADATA_FILE, stringsAsFactors = FALSE, check.names = FALSE)
cat("Raw rows:", nrow(sdrf), "\n")

# Deduplicate — keep one row per patient (SDRF files typically list one row
# per channel/file; keep the Grn IDAT rows only)
sdrf_grn <- sdrf[grep("_Grn.idat$", sdrf[["Array Data File"]]), ]
cat("After deduplication (Grn rows only):", nrow(sdrf_grn), "patients\n")

# --- Edit the column names on the right of each line below to match your file ---
sdrf_grn$age            <- as.numeric(sdrf_grn[["Characteristics[age]"]])
sdrf_grn$sex            <- as.character(sdrf_grn[["Characteristics[sex]"]])
sdrf_grn$response       <- as.character(sdrf_grn[["Characteristics[response to treament]"]])
sdrf_grn$delta_das28    <- as.numeric(sdrf_grn[["Characteristics[deltaDAS28]"]])
sdrf_grn$baseline_das28 <- as.numeric(sdrf_grn[["Characteristics[baseline DAS28]"]])
sdrf_grn$individual     <- as.character(sdrf_grn[["Characteristics[individual]"]])

# Cell-type proportions — if your metadata doesn't already have these
# (e.g. from Houseman deconvolution), either compute them with EpiDISH
# first and merge in, or comment these four lines out and remove
# CD4T/CD8T/Bcell/NK from the model formulas in Sections 8, 10, and 12.
sdrf_grn$CD4T  <- as.numeric(sdrf_grn[["Characteristics[cd4_tlymphocytes_houseman]"]])
sdrf_grn$CD8T  <- as.numeric(sdrf_grn[["Characteristics[cd8_tlymphocytes_houseman]"]])
sdrf_grn$Bcell <- as.numeric(sdrf_grn[["Characteristics[b_lymphocytes_houseman]"]])
sdrf_grn$NK    <- as.numeric(sdrf_grn[["Characteristics[natural_killer_cells_houseman]"]])

# Binary response: responder (good + moderate) = 1, non-responder = 0
# Edit the strings in %in% c(...) to match your response-category labels.
sdrf_grn$is_responder <- ifelse(
  sdrf_grn$response %in% c("good responder", "moderate responder"), 1, 0
)

# IDAT base names (strip _Grn.idat suffix)
sdrf_grn$idat_base <- gsub("_Grn\\.idat$", "", sdrf_grn[["Array Data File"]])

cat("\n=== Sample summary ===\n")
print(table(sdrf_grn$response))
cat("Binary — Responders:", sum(sdrf_grn$is_responder == 1), "\n")
cat("Binary — Non-responders:", sum(sdrf_grn$is_responder == 0), "\n")
cat("Age range:", min(sdrf_grn$age, na.rm = TRUE),
    "to", max(sdrf_grn$age, na.rm = TRUE), "\n")
cat("Sex distribution:\n")
print(table(sdrf_grn$sex))


# =============================================================================
# SECTION 3: LOAD RAW IDAT FILES
# =============================================================================

cat("\nLoading IDAT files from:", IDAT_DIR, "\n")
cat("This takes several minutes depending on cohort size...\n")

idat_files <- list.files(IDAT_DIR, pattern = "\\.idat$")
cat("IDAT files found:", length(idat_files), "\n")

if (length(idat_files) == 0) {
  stop(paste0(
    "\nNo IDAT files found in: ", IDAT_DIR,
    "\nCheck IDAT_DIR in Section 0."
  ))
}

targets <- data.frame(
  Sample_Name    = sdrf_grn$individual,
  Basename       = file.path(IDAT_DIR, sdrf_grn$idat_base),
  Response       = sdrf_grn$response,
  is_responder   = sdrf_grn$is_responder,
  age            = sdrf_grn$age,
  sex            = sdrf_grn$sex,
  delta_das28    = sdrf_grn$delta_das28,
  baseline_das28 = sdrf_grn$baseline_das28,
  CD4T           = sdrf_grn$CD4T,
  CD8T           = sdrf_grn$CD8T,
  Bcell          = sdrf_grn$Bcell,
  NK             = sdrf_grn$NK,
  stringsAsFactors = FALSE
)

basenames_exist <- file.exists(paste0(targets$Basename, "_Grn.idat"))
cat("Basenames with matching IDAT files:", sum(basenames_exist), "/",
    nrow(targets), "\n")

if (sum(basenames_exist) == 0) {
  cat("\nIDAT files not matching expected pattern.\n")
  cat("First expected path:", paste0(targets$Basename[1], "_Grn.idat"), "\n")
  cat("First actual file:", file.path(IDAT_DIR, idat_files[1]), "\n")
  stop("Fix IDAT_DIR path or file naming before continuing.")
}

targets <- targets[basenames_exist, ]

rgset <- read.metharray(basenames = targets$Basename, verbose = TRUE, force = TRUE)
cat("RGChannelSet loaded:", dim(rgset), "\n")


# =============================================================================
# SECTION 4: QUALITY CONTROL AND NORMALIZATION
# =============================================================================

cat("\nRunning QC and normalization...\n")

mset_norm <- preprocessFunnorm(rgset)
cat("After normalization:", dim(mset_norm), "\n")

beta_mtx <- getBeta(mset_norm)
colnames(beta_mtx) <- targets$Sample_Name
cat("Beta matrix dimensions:", dim(beta_mtx), "\n")

na_count <- rowSums(is.na(beta_mtx))
keep_probes <- na_count <= (0.05 * ncol(beta_mtx))
beta_mtx_qc <- beta_mtx[keep_probes, ]
cat("CpGs after QC:", nrow(beta_mtx_qc), "\n")

na_sample <- colSums(is.na(beta_mtx_qc))
keep_samples <- na_sample <= (0.05 * nrow(beta_mtx_qc))
beta_mtx_qc <- beta_mtx_qc[, keep_samples]
targets_qc  <- targets[keep_samples, ]
cat("Samples after QC:", ncol(beta_mtx_qc), "\n")
cat("Responders:", sum(targets_qc$is_responder == 1), "\n")
cat("Non-responders:", sum(targets_qc$is_responder == 0), "\n")


# =============================================================================
# SECTION 5: COMPUTE EPIGENETIC CLOCKS
# =============================================================================

cat("\nComputing epigenetic clocks...\n")

beta_mtx_t <- as.data.frame(t(beta_mtx_qc))

run_clock <- function(fn_name, ...) {
  cat("Computing", fn_name, "... ")
  tryCatch({
    fn  <- get(fn_name, envir = asNamespace("methylCIPHER"))
    res <- fn(beta_mtx_t, ...)
    cat("done\n")
    return(res)
  }, error = function(e) {
    cat("FAILED:", conditionMessage(e), "\n")
    return(NULL)
  })
}

extract_score <- function(res, col_hints) {
  if (is.null(res)) return(rep(NA_real_, nrow(beta_mtx_t)))
  if (is.numeric(res) && length(res) == nrow(beta_mtx_t)) return(as.numeric(res))
  df <- as.data.frame(res)
  for (nm in col_hints) {
    if (nm %in% colnames(df)) return(as.numeric(df[[nm]]))
  }
  num_cols <- which(sapply(df, is.numeric))
  if (length(num_cols) > 0) return(as.numeric(df[[num_cols[1]]]))
  return(rep(NA_real_, nrow(beta_mtx_t)))
}

cat("\n--- First Generation ---\n")
horvath1 <- run_clock("calcHorvath1")
horvath2 <- run_clock("calcHorvath2")
hannum   <- run_clock("calcHannum")

cat("\n--- Second Generation ---\n")
pheno  <- run_clock("calcPhenoAge")
zhang  <- run_clock("calcZhang2019")
lin    <- run_clock("calcLin")
knight <- run_clock("calcKnight")
hrs    <- run_clock("calcHRSInChPhenoAge")
vidal  <- run_clock("calcVidalBralo")

h1s  <- extract_score(horvath1, c("Horvath1", "Horvath", "DNAmAge"))
h2s  <- extract_score(horvath2, c("Horvath2", "DNAmAge"))
hns  <- extract_score(hannum,   c("Hannum", "DNAmAge"))
phs  <- extract_score(pheno,    c("PhenoAge", "DNAmPhenoAge"))
zhs  <- extract_score(zhang,    c("Zhang2019", "Zhang", "DNAmAge"))
lns  <- extract_score(lin,      c("Lin", "DNAmAge"))
kns  <- extract_score(knight,   c("Knight", "DNAmAge"))
hrs_s <- extract_score(hrs,     c("HRSInChPhenoAge", "HRS", "DNAmAge"))
vds  <- extract_score(vidal,    c("VidalBralo", "Vidal", "DNAmAge"))


# =============================================================================
# SECTION 6: BUILD RESULTS TABLE AND AGE ACCELERATION
# =============================================================================

cat("\nBuilding per-patient results table (kept local only)...\n")

results_mtx <- data.frame(
  sample         = colnames(beta_mtx_qc),
  is_responder   = targets_qc$is_responder,
  response       = targets_qc$Response,
  age            = targets_qc$age,
  sex            = targets_qc$sex,
  delta_das28    = targets_qc$delta_das28,
  baseline_das28 = targets_qc$baseline_das28,
  CD4T           = targets_qc$CD4T,
  CD8T           = targets_qc$CD8T,
  Bcell          = targets_qc$Bcell,
  NK             = targets_qc$NK,
  Horvath1       = h1s, Horvath2 = h2s, Hannum = hns,
  PhenoAge       = phs, Zhang2019 = zhs, Lin = lns,
  Knight         = kns, HRSPhenoAge = hrs_s, VidalBralo = vds,
  stringsAsFactors = FALSE
)

calc_accel <- function(clock_vals, age_vals) {
  if (sum(!is.na(clock_vals)) < 5) return(rep(NA_real_, length(clock_vals)))
  resid(lm(clock_vals ~ age_vals, na.action = na.exclude))
}

for (nm in c("Horvath1","Horvath2","Hannum","PhenoAge","Zhang2019",
             "Lin","Knight","HRSPhenoAge","VidalBralo")) {
  results_mtx[[paste0(nm, "_accel")]] <- calc_accel(results_mtx[[nm]], results_mtx$age)
}

clock_accels <- paste0(c("Horvath1","Horvath2","Hannum","PhenoAge","Zhang2019",
                          "Lin","Knight","HRSPhenoAge","VidalBralo"), "_accel")

cat("Results table:", nrow(results_mtx), "patients\n")


# =============================================================================
# SECTION 7: PRIMARY ANALYSIS — DO CLOCKS PREDICT MTX RESPONSE?
# =============================================================================

cat("\n=== PRIMARY ANALYSIS: Do clocks predict MTX response? ===\n")
cat("CRITICAL NOTE: Non-responder n =", sum(results_mtx$is_responder == 0),
    "— interpret with caution if small\n\n")

stat_results <- data.frame()
for (clock in clock_accels) {
  if (!clock %in% colnames(results_mtx)) next
  resp_vals <- results_mtx[[clock]][results_mtx$is_responder == 1]
  nonr_vals <- results_mtx[[clock]][results_mtx$is_responder == 0]
  resp_vals <- resp_vals[!is.na(resp_vals)]
  nonr_vals <- nonr_vals[!is.na(nonr_vals)]
  if (length(resp_vals) < 3 || length(nonr_vals) < 3) next
  tt <- t.test(nonr_vals, resp_vals)
  stat_results <- rbind(stat_results, data.frame(
    Clock = clock,
    Responder_n = length(resp_vals), NonResponder_n = length(nonr_vals),
    Responder_mean = round(mean(resp_vals), 3),
    NonResp_mean = round(mean(nonr_vals), 3),
    Difference = round(mean(nonr_vals) - mean(resp_vals), 3),
    P_value = signif(tt$p.value, 3),
    stringsAsFactors = FALSE
  ))
}
stat_results$FDR_p <- round(p.adjust(stat_results$P_value, "BH"), 5)
stat_results$Significant <- ifelse(stat_results$FDR_p < 0.05, "YES ***",
                             ifelse(stat_results$P_value < 0.05, "Nominal *", "No"))
rownames(stat_results) <- NULL
print(stat_results)


# =============================================================================
# SECTION 8: LOGISTIC REGRESSION — CLOCK + CLINICAL COVARIATES
# =============================================================================

cat("\n=== LOGISTIC REGRESSION: Clock + Clinical Predictors ===\n")

logistic_results <- data.frame()
for (clock in clock_accels) {
  if (!clock %in% colnames(results_mtx)) next
  df_model <- results_mtx[, c("is_responder", clock, "age", "sex",
                               "baseline_das28", "CD4T", "CD8T", "Bcell", "NK")]
  df_model <- df_model[complete.cases(df_model), ]
  df_model$sex_bin <- as.integer(grepl("female", df_model$sex, ignore.case = TRUE))
  if (nrow(df_model) < 15) next

  tryCatch({
    mod_clin <- glm(is_responder ~ age + sex_bin + baseline_das28 +
                       CD4T + CD8T + Bcell + NK, data = df_model, family = binomial)
    roc_clin <- roc(df_model$is_responder, fitted(mod_clin), quiet = TRUE)

    formula_str <- paste("is_responder ~", clock,
                          "+ age + sex_bin + baseline_das28 + CD4T + CD8T + Bcell + NK")
    mod_full <- glm(as.formula(formula_str), data = df_model, family = binomial)
    roc_full <- roc(df_model$is_responder, fitted(mod_full), quiet = TRUE)

    logistic_results <- rbind(logistic_results, data.frame(
      Clock = clock, N = nrow(df_model),
      AUC_clinical = round(as.numeric(auc(roc_clin)), 3),
      AUC_with_clock = round(as.numeric(auc(roc_full)), 3),
      AUC_improvement = round(as.numeric(auc(roc_full)) - as.numeric(auc(roc_clin)), 3),
      stringsAsFactors = FALSE
    ))
  }, error = function(e) cat("Logistic model failed for", clock, ":", conditionMessage(e), "\n"))
}
if (nrow(logistic_results) > 0) {
  logistic_results <- logistic_results[order(-logistic_results$AUC_with_clock), ]
  rownames(logistic_results) <- NULL
  print(logistic_results)
}


# =============================================================================
# SECTION 9: CONTINUOUS ANALYSIS — CLOCK vs DELTADAS28
# =============================================================================

cat("\n=== CONTINUOUS ANALYSIS: Clock vs DeltaDAS28 ===\n")

corr_results <- data.frame()
for (clock in clock_accels) {
  if (!clock %in% colnames(results_mtx)) next
  df_corr <- results_mtx[complete.cases(results_mtx[[clock]], results_mtx$delta_das28), ]
  if (nrow(df_corr) < 10) next
  ct <- cor.test(df_corr[[clock]], df_corr$delta_das28, method = "spearman")
  corr_results <- rbind(corr_results, data.frame(
    Clock = clock, N = nrow(df_corr),
    Spearman_r = round(ct$estimate, 3), P_value = signif(ct$p.value, 3),
    stringsAsFactors = FALSE
  ))
}
corr_results$FDR_p <- round(p.adjust(corr_results$P_value, "BH"), 5)
corr_results$Direction <- ifelse(corr_results$Spearman_r < 0,
  "Higher acceleration = better response", "Higher acceleration = worse response")
corr_results <- corr_results[order(corr_results$P_value), ]
rownames(corr_results) <- NULL
print(corr_results)


# =============================================================================
# SECTION 10: MWAS — BASELINE METHYLATION vs MTX RESPONSE
# =============================================================================

cat("\n=== MWAS: Baseline Methylation vs MTX Response ===\n")

sex_bin_vec <- as.integer(grepl("female", targets_qc$sex, ignore.case = TRUE))

design_mtx <- model.matrix(
  ~ targets_qc$is_responder + targets_qc$age + sex_bin_vec +
    targets_qc$baseline_das28 + targets_qc$CD4T + targets_qc$CD8T +
    targets_qc$Bcell + targets_qc$NK
)

cat("Running limma MWAS on", nrow(beta_mtx_qc), "CpGs (several minutes)...\n")

fit_mtx    <- lmFit(beta_mtx_qc, design_mtx)
fit_mtx_eb <- eBayes(fit_mtx)
mwas_mtx <- topTable(fit_mtx_eb, coef = 2, number = Inf,
                      adjust.method = "BH", sort.by = "p")

cat("Total CpGs tested:", nrow(mwas_mtx), "\n")
cat("FDR significant (q<0.05):", sum(mwas_mtx$adj.P.Val < 0.05), "\n")
cat("Nominal significant (p<0.05):", sum(mwas_mtx$P.Value < 0.05), "\n")


# =============================================================================
# SECTION 11: DMARD RESPONSE SCORE (DRS)
# =============================================================================

cat("\nBuilding DMARD Response Score (DRS)...\n")

n_features <- min(500, nrow(mwas_mtx))
feature_cpgs <- rownames(mwas_mtx)[1:n_features]

X_mtx <- t(beta_mtx_qc[feature_cpgs, ])
y_mtx <- targets_qc$is_responder
complete_mtx <- complete.cases(X_mtx) & !is.na(y_mtx)
X_mtx <- X_mtx[complete_mtx, ]
y_mtx <- y_mtx[complete_mtx]
cpg_ok <- apply(X_mtx, 2, function(x) sum(is.na(x)) == 0)
X_mtx <- X_mtx[, cpg_ok]

cat("Samples:", nrow(X_mtx), "| Features:", ncol(X_mtx), "\n")

set.seed(42)
cv_mtx <- cv.glmnet(X_mtx, y_mtx, family = "binomial", alpha = 0.5,
                     nfolds = 3, type.measure = "deviance")
final_mtx <- glmnet(X_mtx, y_mtx, family = "binomial", alpha = 0.5,
                     lambda = cv_mtx$lambda.min)

set.seed(42)
fold_ids <- sample(rep(1:3, length.out = nrow(X_mtx)))
cv_preds <- numeric(nrow(X_mtx))
for (fold in 1:3) {
  test_idx <- fold_ids == fold
  train_idx <- !test_idx
  if (sum(train_idx) < 5 || sum(test_idx) < 3) next
  fold_mod <- glmnet(X_mtx[train_idx, ], y_mtx[train_idx],
                      family = "binomial", alpha = 0.5, lambda = cv_mtx$lambda.min)
  cv_preds[test_idx] <- as.numeric(predict(fold_mod, newx = X_mtx[test_idx, ], type = "response"))
}

roc_drs <- roc(y_mtx, cv_preds, quiet = TRUE)
auc_drs <- auc(roc_drs)
ci_drs  <- ci.auc(roc_drs)
coef_drs <- coef(final_mtx)
drs_cpgs <- names(coef_drs[coef_drs[, 1] != 0, ])[-1]

cat("DRS CpGs selected:", length(drs_cpgs), "\n")
cat("DRS cross-validated AUC:", round(auc_drs, 3),
    "(95% CI:", round(ci_drs[1], 3), "-", round(ci_drs[3], 3), ")\n")


# =============================================================================
# SECTION 12: HEAD-TO-HEAD — CLOCKS vs DRS vs CLINICAL
# =============================================================================

cat("\n=== HEAD-TO-HEAD: DRS vs Clocks vs Clinical ===\n")

comparison_mtx <- data.frame(
  Model = "DRS (MWAS composite)", AUC = round(as.numeric(auc_drs), 3),
  CI_low = round(ci_drs[1], 3), CI_high = round(ci_drs[3], 3),
  Type = "MWAS-derived Score"
)

for (clock in clock_accels) {
  if (!clock %in% colnames(results_mtx)) next
  idx <- complete.cases(results_mtx[[clock]], results_mtx$is_responder)
  if (sum(idx) < 10) next
  roc_obj <- roc(results_mtx$is_responder[idx], results_mtx[[clock]][idx], quiet = TRUE)
  ci_obj  <- ci.auc(roc_obj)
  comparison_mtx <- rbind(comparison_mtx, data.frame(
    Model = clock, AUC = round(as.numeric(auc(roc_obj)), 3),
    CI_low = round(ci_obj[1], 3), CI_high = round(ci_obj[3], 3),
    Type = "Epigenetic Clock"
  ))
}

df_clin <- results_mtx[complete.cases(results_mtx$is_responder, results_mtx$baseline_das28,
                                       results_mtx$age, results_mtx$sex), ]
df_clin$sex_bin <- as.integer(grepl("female", df_clin$sex, ignore.case = TRUE))
tryCatch({
  mod_clin_only <- glm(is_responder ~ baseline_das28 + age + sex_bin,
                        data = df_clin, family = binomial)
  roc_clin_only <- roc(df_clin$is_responder, fitted(mod_clin_only), quiet = TRUE)
  ci_clin_only  <- ci.auc(roc_clin_only)
  comparison_mtx <- rbind(comparison_mtx, data.frame(
    Model = "Clinical Only (DAS28 + Age + Sex)",
    AUC = round(as.numeric(auc(roc_clin_only)), 3),
    CI_low = round(ci_clin_only[1], 3), CI_high = round(ci_clin_only[3], 3),
    Type = "Clinical Baseline"
  ))
}, error = function(e) cat("Clinical model failed:", conditionMessage(e), "\n"))

comparison_mtx <- comparison_mtx[order(-comparison_mtx$AUC), ]
rownames(comparison_mtx) <- NULL
print(comparison_mtx)


# =============================================================================
# SECTION 13: PLOTS (aggregate/statistical plots only — safe to share)
# =============================================================================

cat("\nGenerating plots...\n")

resp_colors <- c("good responder" = "#1F4E79", "moderate responder" = "#2E75B6",
                  "non responder" = "#C00000")

plot_1 <- ggplot(results_mtx[!is.na(results_mtx$PhenoAge_accel), ],
                  aes(x = response, y = PhenoAge_accel, fill = response)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = resp_colors) +
  labs(title = "PhenoAge Acceleration by MTX Response Group",
       x = "MTX Response at 3 Months", y = "PhenoAge Acceleration (years)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 15, hjust = 1))
ggsave(file.path(OUTPUT_DIR, "PhenoAge_by_response.png"), plot_1, width = 7, height = 5, dpi = 300)

plot_2 <- ggplot(comparison_mtx, aes(x = reorder(Model, AUC), y = AUC, fill = Type)) +
  geom_col(alpha = 0.85) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.3, color = "gray30") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "darkgreen") +
  coord_flip() +
  scale_fill_manual(values = c("MWAS-derived Score" = "#C00000",
                                "Epigenetic Clock" = "#1F4E79",
                                "Clinical Baseline" = "#7F7F7F")) +
  labs(title = "Head-to-Head: DRS vs Clocks vs Clinical — MTX Response",
       x = "Model", y = "AUC (95% CI)") +
  theme_classic(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(OUTPUT_DIR, "HeadToHead_AUC.png"), plot_2, width = 10, height = 7, dpi = 300)

mwas_plot_df <- data.frame(
  logFC = mwas_mtx$logFC, negLogP = -log10(mwas_mtx$P.Value), FDR = mwas_mtx$adj.P.Val
)
mwas_plot_df$sig <- ifelse(mwas_plot_df$FDR < 0.05, "FDR<0.05",
                     ifelse(mwas_plot_df$negLogP > -log10(0.05), "Nominal p<0.05", "NS"))
plot_3 <- ggplot(mwas_plot_df, aes(x = logFC, y = negLogP, color = sig)) +
  geom_point(alpha = 0.3, size = 0.5) +
  scale_color_manual(values = c("FDR<0.05" = "#C00000", "Nominal p<0.05" = "#FFC000", "NS" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "orange") +
  labs(title = "MWAS Volcano: Baseline Methylation vs MTX Response",
       x = "log Fold Change (Non-responder vs Responder)", y = expression(-log[10](p-value))) +
  theme_classic(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(OUTPUT_DIR, "MWAS_Volcano.png"), plot_3, width = 7, height = 6, dpi = 300)

cat("Plots saved to:", OUTPUT_DIR, "\n")
cat("These plots show group-level/aggregate patterns only — no plot contains\n")
cat("per-patient identifiers, so they are safe to include when sharing results.\n")


# =============================================================================
# SECTION 14: SAVE RESULTS — SPLIT INTO LOCAL-ONLY vs SAFE-TO-SHARE
# =============================================================================

cat("\nSaving results...\n")

# --- KEEP LOCAL: contains per-patient values and the methylation matrix ---
save(results_mtx, beta_mtx_qc, targets_qc,
     file = file.path(OUTPUT_DIR, "internal_full_results.RData"))
cat("Saved (DO NOT SHARE):", file.path(OUTPUT_DIR, "internal_full_results.RData"), "\n")

# --- SAFE TO SHARE: aggregate statistics only, no per-patient/per-CpG values ---
mwas_summary <- data.frame(
  total_CpGs_tested = nrow(mwas_mtx),
  FDR_significant_CpGs = sum(mwas_mtx$adj.P.Val < 0.05),
  nominal_significant_CpGs = sum(mwas_mtx$P.Value < 0.05),
  genome_wide_significant_p_lt_1e7 = sum(mwas_mtx$P.Value < 1e-7)
)

drs_summary <- data.frame(
  n_CpGs_selected = length(drs_cpgs),
  cv_AUC = round(as.numeric(auc_drs), 3),
  CI_low = round(ci_drs[1], 3),
  CI_high = round(ci_drs[3], 3)
)

cohort_summary <- data.frame(
  n_total = nrow(results_mtx),
  n_responders = sum(results_mtx$is_responder == 1),
  n_non_responders = sum(results_mtx$is_responder == 0),
  age_min = round(min(results_mtx$age, na.rm = TRUE), 1),
  age_max = round(max(results_mtx$age, na.rm = TRUE), 1)
)

save(cohort_summary, stat_results, logistic_results, corr_results,
     comparison_mtx, mwas_summary, drs_summary,
     file = file.path(OUTPUT_DIR, "summary_results_SHARE.RData"))

# Also write plain CSVs for easy inspection without opening R
write.csv(cohort_summary, file.path(OUTPUT_DIR, "SHARE_cohort_summary.csv"), row.names = FALSE)
write.csv(stat_results, file.path(OUTPUT_DIR, "SHARE_clock_vs_response_ttests.csv"), row.names = FALSE)
write.csv(logistic_results, file.path(OUTPUT_DIR, "SHARE_logistic_clock_plus_clinical.csv"), row.names = FALSE)
write.csv(corr_results, file.path(OUTPUT_DIR, "SHARE_clock_vs_deltaDAS28_correlation.csv"), row.names = FALSE)
write.csv(comparison_mtx, file.path(OUTPUT_DIR, "SHARE_headtohead_AUC_comparison.csv"), row.names = FALSE)
write.csv(mwas_summary, file.path(OUTPUT_DIR, "SHARE_mwas_summary.csv"), row.names = FALSE)
write.csv(drs_summary, file.path(OUTPUT_DIR, "SHARE_drs_summary.csv"), row.names = FALSE)

cat("\n=== SAFE TO SHARE ===\n")
cat("The following files contain ONLY aggregate statistics and are safe to\n")
cat("send back (plus the three PNG plots in", OUTPUT_DIR, "):\n")
cat(" -", file.path(OUTPUT_DIR, "summary_results_SHARE.RData"), "\n")
cat(" - SHARE_*.csv files in", OUTPUT_DIR, "\n")
cat("\nDo NOT share: internal_full_results.RData (contains per-patient data)\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
