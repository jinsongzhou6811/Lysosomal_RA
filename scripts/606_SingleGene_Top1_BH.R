# ========================================================
# 606_SingleGene_Top1_Steiger_BH.R  (v3)
# Single gene: select 1 strongest outcome per gene
#
# BH file:   {single_dir}/{gene}/{gene}_BH_corrected.csv
# Result folder (after 605): {single_dir}/{gene}/{gene}_{outcome}/
#
# QC logic (in order of strictness):
#   HARD filters (exclusion):
#     1. BH p-value < 0.05
#     2. F_pass (mean_F > 10, from BH file)
#     3. mean_F > 20 (additional threshold)
#     4. het_pass  (IVW Q-test p > 0.05, from BH file)
#     5. nsnp >= 3
#   SOFT flags (recorded but NOT used for exclusion):
#     6. egger_pass  → recorded in output as egger_flag
#     7. MR-PRESSO:  if significant (p<=0.05), record flag;
#                    use PRESSO-corrected OR if available
#     8. Direction consistency: record % agreement, do NOT exclude
#   Steiger: REMOVED (units="unit" causes universal failure)
#
# Output: outputs/06_UVMR_Analysis/04_UVMR_Data_filtering/SingleGene_Top1/
# ========================================================

library(rprojroot)
library(fs)
library(data.table)
library(dplyr)
library(stringr)

project_root <- find_rstudio_root_file()
single_dir   <- file.path(project_root, "outputs", "06_UVMR_Analysis", "02_UVMR_Single_Gene")

output_dir <- file.path(project_root, "outputs", "06_UVMR_Analysis",
                        "04_UVMR_Data_filtering", "SingleGene_Top1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(single_dir)) stop("Single gene directory not found: ", single_dir)

all_single_folders <- fs::dir_ls(single_dir, type = "directory")

results_list   <- list()
skipped_no_bh  <- 0L
skipped_no_sig <- 0L
skipped_qc     <- 0L
added          <- 0L

for (folder in all_single_folders) {
  
  exposure <- basename(folder)
  
  # ── (1) BH file ──────────────────────────────────────────────
  bh_file <- file.path(folder, paste0(exposure, "_BH_corrected.csv"))
  if (!file.exists(bh_file)) {
    cat("  [SKIP] No BH file:", exposure, "\n")
    skipped_no_bh <- skipped_no_bh + 1L
    next
  }
  bh_dt <- fread(bh_file)
  
  # ── (2) BH significance filter ───────────────────────────────
  sig_rows <- bh_dt[!is.na(ivw_pval_BH) & ivw_pval_BH < 0.05]
  if (nrow(sig_rows) == 0) {
    cat("  [SKIP] No BH-sig outcomes:", exposure, "\n")
    skipped_no_sig <- skipped_no_sig + 1L
    next
  }
  
  cat("-> Processing:", exposure, "| BH-sig:", nrow(sig_rows), "\n")
  
  for (i in seq_len(nrow(sig_rows))) {
    row     <- sig_rows[i]
    outcome <- row$outcome_id
    
    # ── (3) HARD QC: F statistic ─────────────────────────────
    f_ok <- isTRUE(row$F_pass) && !is.na(row$mean_F) && row$mean_F > 20
    if (!f_ok) {
      cat("    [HARD FAIL - F] outcome:", outcome,
          "| F_pass:", row$F_pass, "mean_F:", round(row$mean_F, 1), "\n")
      skipped_qc <- skipped_qc + 1L
      next
    }
    
    # ── (4) HARD QC: Heterogeneity ───────────────────────────
    het_ok <- isTRUE(row$het_pass)
    if (!het_ok) {
      cat("    [HARD FAIL - het] outcome:", outcome,
          "| het_Q_pval:", round(row$het_Q_pval, 4), "\n")
      skipped_qc <- skipped_qc + 1L
      next
    }
    
    # ── (5) Locate result folder ──────────────────────────────
    out_folder <- file.path(folder, paste0(exposure, "_", outcome))
    if (!dir.exists(out_folder)) {
      cat("    [SKIP] Folder missing:", basename(out_folder), "\n")
      next
    }
    
    # ── (6) Load MR results ───────────────────────────────────
    mr_res <- tryCatch(
      fread(file.path(out_folder, "mr_results.csv")),
      error = function(e) { cat("    [SKIP] mr_results error\n"); NULL }
    )
    if (is.null(mr_res)) next
    
    filtered <- tryCatch(
      fread(file.path(out_folder, "filtered_data.csv")),
      error = function(e) { cat("    [SKIP] filtered_data error\n"); NULL }
    )
    if (is.null(filtered)) next
    
    # ── (7) HARD QC: nsnp ────────────────────────────────────
    if (nrow(filtered) < 3) {
      cat("    [HARD FAIL - nsnp <3]", basename(out_folder), "\n")
      skipped_qc <- skipped_qc + 1L
      next
    }
    
    # ── (8) IVW row ───────────────────────────────────────────
    ivw_row <- mr_res[method == "Inverse variance weighted"]
    if (nrow(ivw_row) == 0) {
      cat("    [SKIP] No IVW result\n")
      next
    }
    
    # ── SOFT flags (record only, do NOT exclude) ──────────────
    
    # Egger intercept flag
    egger_flag <- !isTRUE(row$egger_pass)
    
    # MR-PRESSO flag + corrected OR if available
    presso_sig    <- FALSE
    presso_b      <- NA_real_
    presso_file   <- file.path(out_folder, "mr_presso_main.csv")
    if (file.exists(presso_file)) {
      presso_res <- tryCatch(fread(presso_file), error = function(e) NULL)
      if (!is.null(presso_res) && nrow(presso_res) > 0) {
        pval_val <- presso_res$`P-value`[1]
        presso_sig <- !is.na(pval_val) && pval_val <= 0.05
        # Try to get corrected estimate (Outlier-corrected row)
        if (presso_sig && nrow(presso_res) >= 2) {
          corr_row <- presso_res[grepl("Outlier", presso_res$`MR Analysis`,
                                       ignore.case = TRUE)]
          if (nrow(corr_row) > 0)
            presso_b <- as.numeric(corr_row$`Causal Estimate`[1])
        }
      }
    }
    
    # Direction consistency (soft, just record ratio)
    non_ivw_b      <- mr_res$b[mr_res$method != "Inverse variance weighted" &
                                 !is.na(mr_res$b)]
    n_agree        <- sum(sign(non_ivw_b) == sign(ivw_row$b))
    n_methods      <- length(non_ivw_b)
    direction_ratio <- if (n_methods > 0) n_agree / n_methods else NA_real_
    
    # ── Composite score ───────────────────────────────────────
    mean_F <- row$mean_F
    score  <- -log10(row$ivw_pval_BH + 1e-300) *
      abs(ivw_row$b) / ivw_row$se *
      mean_F * 1.2
    
    results_list[[paste(exposure, outcome, sep = "_")]] <- data.table(
      exposure             = exposure,
      outcome              = outcome,
      b                    = ivw_row$b,
      se                   = ivw_row$se,
      pval_raw             = ivw_row$pval,
      pval_BH              = row$ivw_pval_BH,
      OR                   = row$OR,
      OR_lower95           = row$OR_lower95,
      OR_upper95           = row$OR_upper95,
      nsnp                 = nrow(filtered),
      mean_F               = mean_F,
      het_Q_pval           = row$het_Q_pval,
      egger_intercept_pval = row$egger_intercept_pval,
      egger_flag           = egger_flag,         # SOFT: TRUE = concern
      presso_sig           = presso_sig,          # SOFT: TRUE = outlier detected
      presso_corrected_b   = presso_b,            # SOFT: corrected estimate if available
      direction_ratio      = direction_ratio,      # SOFT: proportion of methods agreeing
      all_QC_pass          = row$all_QC_pass,
      score                = score
    )
    added <- added + 1L
    cat("    [OK]", outcome,
        "| egger_flag:", egger_flag,
        "| presso_sig:", presso_sig,
        "| dir:", round(direction_ratio, 2),
        "| score:", round(score, 2), "\n")
  }
}

# ── Summary ──────────────────────────────────────────────────────
cat("\n=== Summary ===\n")
cat("Skipped (no BH file):        ", skipped_no_bh,  "\n")
cat("Skipped (no BH-sig outcome): ", skipped_no_sig, "\n")
cat("Skipped (hard QC fail):      ", skipped_qc,     "\n")
cat("Passed hard filters:         ", added,           "\n")

if (length(results_list) == 0) {
  cat("No results passed hard filters.\n")
} else {
  all_top <- rbindlist(results_list, fill = TRUE)
  
  # Save full passing table
  all_file <- file.path(output_dir, "SingleGene_AllPassing.csv")
  fwrite(all_top %>% arrange(desc(score)), all_file)
  cat("All passing ->", all_file, "\n")
  
  # Top1 per gene (highest score)
  one_per_gene <- all_top %>%
    group_by(exposure) %>%
    slice_max(score, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(score))
  
  out_file <- file.path(output_dir, "SingleGene_Top1_All.csv")
  fwrite(one_per_gene, out_file)
  cat("Top1 per gene (", nrow(one_per_gene), "genes) ->", out_file, "\n")
  
  # Print clean summary table
  cat("\n=== Top results (sorted by score) ===\n")
  print(all_top %>%
          as.data.frame() %>%
          arrange(desc(score)) %>%
          dplyr::select(exposure, outcome, OR, OR_lower95, OR_upper95,
                        pval_raw, pval_BH, mean_F, egger_flag, presso_sig,
                        direction_ratio, score) %>%
          head(20))
}

cat("606 Single-gene Top1 completed!\n")