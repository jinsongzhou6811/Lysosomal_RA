# ========================================================
# 607_ComboWeights_Top10_Per_Module_Steiger_BH.R  (v3)
# Combo-weights: 1 strongest outcome per combo -> Top10 per module
#
# BH file:   {weights_dir}/{combo}/{combo}_weights_BH_corrected.csv
# Result folder (after 605): {weights_dir}/{combo}/{combo}_weights_{outcome}/
#
# QC logic  — same as 606 v3:
#   HARD filters (exclusion):
#     1. BH p-value < 0.05
#     2. F_pass && mean_F > 20
#     3. het_pass
#     4. nsnp >= 3
#   SOFT flags (recorded, NOT used for exclusion):
#     5. egger_pass
#     6. MR-PRESSO significance + corrected estimate
#     7. Direction consistency ratio
#   Steiger: REMOVED
#
# Output: outputs/06_UVMR_Analysis/04_UVMR_Data_filtering/ComboWeights_Top10/
# ========================================================

library(rprojroot)
library(fs)
library(data.table)
library(dplyr)
library(stringr)

project_root <- find_rstudio_root_file()
weights_dir  <- file.path(project_root, "outputs", "06_UVMR_Analysis", "03_UVMR_Combo_Weights")

output_dir <- file.path(project_root, "outputs", "06_UVMR_Analysis",
                        "04_UVMR_Data_filtering", "ComboWeights_Top10")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(weights_dir)) stop("Weights directory not found: ", weights_dir)

get_module <- function(name) {
  if (str_detect(name, "GZMB")) return("GZMB")
  if (str_detect(name, "NKG7")) return("NKG7")
  if (str_detect(name, "LAMP3")) return("LAMP3")
  return("Unknown")
}

all_combo_folders <- fs::dir_ls(weights_dir, type = "directory")

results_list   <- list()
skipped_no_bh  <- 0L
skipped_no_sig <- 0L
skipped_qc     <- 0L
added          <- 0L

for (top_folder in all_combo_folders) {
  
  combo_name  <- basename(top_folder)
  module_name <- get_module(combo_name)
  
  # ── (1) BH file ──────────────────────────────────────────────
  bh_file <- file.path(top_folder, paste0(combo_name, "_weights_BH_corrected.csv"))
  if (!file.exists(bh_file)) {
    cat("  [SKIP] No BH file:", combo_name, "\n")
    skipped_no_bh <- skipped_no_bh + 1L
    next
  }
  bh_dt <- fread(bh_file)
  
  # ── (2) BH significance filter ───────────────────────────────
  sig_rows <- bh_dt[!is.na(ivw_pval_BH) & ivw_pval_BH < 0.05]
  if (nrow(sig_rows) == 0) {
    cat("  [SKIP] No BH-sig outcomes:", combo_name, "\n")
    skipped_no_sig <- skipped_no_sig + 1L
    next
  }
  
  cat("-> Processing:", combo_name, "| module:", module_name,
      "| BH-sig:", nrow(sig_rows), "\n")
  
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
    out_folder <- file.path(top_folder, paste0(combo_name, "_weights_", outcome))
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
    
    # ── SOFT flags ────────────────────────────────────────────
    
    egger_flag <- !isTRUE(row$egger_pass)
    
    presso_sig  <- FALSE
    presso_b    <- NA_real_
    presso_file <- file.path(out_folder, "mr_presso_main.csv")
    if (file.exists(presso_file)) {
      presso_res <- tryCatch(fread(presso_file), error = function(e) NULL)
      if (!is.null(presso_res) && nrow(presso_res) > 0) {
        pval_val   <- presso_res$`P-value`[1]
        presso_sig <- !is.na(pval_val) && pval_val <= 0.05
        if (presso_sig && nrow(presso_res) >= 2) {
          corr_row <- presso_res[grepl("Outlier", presso_res$`MR Analysis`,
                                       ignore.case = TRUE)]
          if (nrow(corr_row) > 0)
            presso_b <- as.numeric(corr_row$`Causal Estimate`[1])
        }
      }
    }
    
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
    
    results_list[[paste(combo_name, outcome, sep = "_weights_")]] <- data.table(
      module               = module_name,
      exposure             = combo_name,
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
      egger_flag           = egger_flag,
      presso_sig           = presso_sig,
      presso_corrected_b   = presso_b,
      direction_ratio      = direction_ratio,
      all_QC_pass          = row$all_QC_pass,
      score                = score
    )
    added <- added + 1L
    cat("    [OK]", outcome,
        "| egger:", egger_flag,
        "| presso:", presso_sig,
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
  cat("No results passed hard filters. Nothing saved.\n")
  cat("607 Combo-weights Top10 completed!\n")
  stop("No results.", call. = FALSE)
}

all_top <- rbindlist(results_list, fill = TRUE)

# ── Per-exposure Top1 ─────────────────────────────────────────────
one_per_exposure <- all_top %>%
  group_by(exposure) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>%
  ungroup()

# ── Per-module Top10 ──────────────────────────────────────────────
all_modules <- unique(c("GZMB", "NKG7",
                        unique(one_per_exposure$module)))

for (mod in all_modules) {
  mod_data <- one_per_exposure %>% filter(module == mod)
  if (nrow(mod_data) == 0) next
  top10    <- mod_data %>% arrange(desc(score)) %>% slice_head(n = 10)
  out_file <- file.path(output_dir, paste0("Top10_", mod, ".csv"))
  fwrite(top10, out_file)
  cat("Module", mod, "(", nrow(top10), ") ->", out_file, "\n")
}

# ── Global summary ────────────────────────────────────────────────
global_top10 <- one_per_exposure %>%
  group_by(module) %>%
  arrange(desc(score)) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  arrange(module, desc(score))

fwrite(global_top10, file.path(output_dir, "All_Modules_Top10_Summary.csv"))
cat("Global summary ->", file.path(output_dir, "All_Modules_Top10_Summary.csv"), "\n")

# Full passing table
fwrite(all_top %>% arrange(module, desc(score)),
       file.path(output_dir, "ComboWeights_AllPassing.csv"))
cat("All passing ->", file.path(output_dir, "ComboWeights_AllPassing.csv"), "\n")

# ── Print top results ─────────────────────────────────────────────
cat("\n=== Top results per module (sorted by score) ===\n")
print(one_per_exposure %>%
        as.data.frame() %>%
        arrange(module, desc(score)) %>%
        dplyr::select(module, exposure, outcome, OR, OR_lower95, OR_upper95,
                      pval_BH, mean_F, egger_flag, presso_sig, direction_ratio, score))

cat("607 Combo-weights Top10 completed!\n")