# ========================================================
# Enhanced UVMR Results Organizer + Cleaner (v3.3.2) - Weights Version
# 修复：跳过 logs 文件夹 + 避免 Windows 长路径删除失败
# ========================================================
library(rprojroot)
library(fs)
library(stringr)
library(dplyr)

# ====================== Smart Path Management ======================
project_root <- find_rstudio_root_file()
single_gene_dir <- file.path(project_root, "outputs", "06_UVMR_Analysis", "02_UVMR_Single_Gene")
weights_dir <- file.path(project_root, "outputs", "06_UVMR_Analysis", "03_UVMR_Combo_Weights")

if (!dir_exists(single_gene_dir)) stop("Single gene directory does not exist: ", single_gene_dir)
if (!dir_exists(weights_dir)) stop("Weights directory does not exist: ", weights_dir)

# ====================== Helper: Skip logs ======================
skip_logs <- function(folder) {
  if (path_file(folder) == "logs") {
    cat(" → Skipped logs directory: ", path_file(folder), "\n")
    return(TRUE)
  }
  return(FALSE)
}

# ====================== Main Logic: Single Gene ======================
success_single <- 0; fail_single <- 0; cleaned_single <- 0
all_single_folders <- fs::dir_ls(single_gene_dir, type = "directory")
for (folder in all_single_folders) {
  folder_name <- fs::path_file(folder)
  if (!str_starts(folder_name, "exposure_")) next
  
  gene_name <- str_remove(folder_name, "^exposure_")
  
  if (!str_detect(gene_name, "_")) {
    cat("→ Processing single gene: ", gene_name, "\n")
  } else {
    cat("⚠ Skipped (invalid format): ", folder_name, "\n")
    next
  }
  
  target_dir <- fs::path(single_gene_dir, gene_name)
  if (fs::dir_exists(target_dir)) fs::dir_delete(target_dir)
  fs::file_move(folder, target_dir)
  
  if (fs::dir_exists(target_dir)) {
    success_single <- success_single + 1
    cat("✓ Rename successful: ", folder_name, "→", gene_name, "\n")
  } else {
    fail_single <- fail_single + 1
    cat("✗ Rename failed: ", folder_name, "\n")
    next
  }
  
  # Step 2: Remove "results_" prefix
  sub_folders <- fs::dir_ls(target_dir, type = "directory", recurse = FALSE)
  for (sub in sub_folders) {
    sub_name <- fs::path_file(sub)
    if (str_starts(sub_name, "results_")) {
      new_name <- str_remove(sub_name, "^results_")
      if (!str_starts(new_name, paste0(gene_name, "_"))) next
      
      new_path <- fs::path(fs::path_dir(sub), new_name)
      if (fs::dir_exists(new_path)) fs::dir_delete(new_path)
      fs::file_move(sub, new_path)
      cat("  Subfolder renamed: ", sub_name, "→", new_name, "\n")
    }
  }
  
  # Step 3: Clean incomplete folders (跳过 logs)
  outcome_folders <- fs::dir_ls(target_dir, type = "directory", recurse = FALSE)
  for (outcome_folder in outcome_folders) {
    if (skip_logs(outcome_folder)) next
    
    if (length(fs::dir_ls(outcome_folder, type = "directory")) == 0) {
      existing_files <- basename(fs::dir_ls(outcome_folder, regexp = "\\.csv$"))
      required_files <- c("mr_results.csv", "heterogeneity_results.csv",
                          "egger_intercept.csv", "exposure_SNPs.csv")
      
      if (!all(required_files %in% existing_files)) {
        fs::dir_delete(outcome_folder)
        cleaned_single <- cleaned_single + 1
        cat(" ❌ Deleted incomplete folder: ", fs::path_file(outcome_folder), "\n")
      }
    }
  }
}

# ====================== Main Logic: Weights ======================
success_weights <- 0; fail_weights <- 0; cleaned_weights <- 0
all_weights_folders <- fs::dir_ls(weights_dir, type = "directory")
for (folder in all_weights_folders) {
  folder_name <- fs::path_file(folder)
  if (!str_starts(folder_name, "exposure_")) next
  
  gene_combo_full <- str_remove(folder_name, "^exposure_")
  gene_combo <- str_remove(gene_combo_full, "_weights$")
  
  if (str_ends(gene_combo_full, "_weights")) {
    cat("→ Processing Weights combination: ", gene_combo, "\n")
  } else {
    cat("⚠ Skipped (invalid format): ", folder_name, "\n")
    next
  }
  
  target_dir <- fs::path(weights_dir, gene_combo)
  if (fs::dir_exists(target_dir)) fs::dir_delete(target_dir)
  fs::file_move(folder, target_dir)
  
  if (fs::dir_exists(target_dir)) {
    success_weights <- success_weights + 1
    cat("✓ Rename successful: ", folder_name, "→", gene_combo, "\n")
  } else {
    fail_weights <- fail_weights + 1
    cat("✗ Rename failed: ", folder_name, "\n")
    next
  }
  
  # Step 2: Remove "results_" prefix
  sub_folders <- fs::dir_ls(target_dir, type = "directory", recurse = FALSE)
  for (sub in sub_folders) {
    sub_name <- fs::path_file(sub)
    if (str_starts(sub_name, "results_")) {
      new_name <- str_remove(sub_name, "^results_")
      if (!str_starts(new_name, paste0(gene_combo, "_weights_"))) next
      
      new_path <- fs::path(fs::path_dir(sub), new_name)
      if (fs::dir_exists(new_path)) fs::dir_delete(new_path)
      fs::file_move(sub, new_path)
      cat("  Subfolder renamed: ", sub_name, "→", new_name, "\n")
    }
  }
  
  # Step 3: Clean incomplete folders (跳过 logs)
  outcome_folders <- fs::dir_ls(target_dir, type = "directory", recurse = FALSE)
  for (outcome_folder in outcome_folders) {
    if (skip_logs(outcome_folder)) next
    
    if (length(fs::dir_ls(outcome_folder, type = "directory")) == 0) {
      existing_files <- basename(fs::dir_ls(outcome_folder, regexp = "\\.csv$"))
      required_files <- c("mr_results.csv",
                          "heterogeneity_results.csv",
                          "egger_intercept.csv",
                          "exposure_SNPs.csv")
      
      if (!all(required_files %in% existing_files)) {
        fs::dir_delete(outcome_folder)
        cleaned_weights <- cleaned_weights + 1
        cat(" ❌ Deleted incomplete folder: ", fs::path_file(outcome_folder),
            " (missing required files)\n")
      }
    }
  }
}

# ====================== Final Statistics ======================
cat("\n=== Processing completed! ===\n")
cat("Single Gene - Successfully renamed: ", success_single, "\n")
cat("Single Gene - Rename failed: ", fail_single, "\n")
cat("Single Gene - Cleaned incomplete: ", cleaned_single, "\n")
cat("Weights - Successfully renamed: ", success_weights, "\n")
cat("Weights - Rename failed: ", fail_weights, "\n")
cat("Weights - Cleaned incomplete: ", cleaned_weights, "\n")
cat("Single gene directory: ", single_gene_dir, "\n")
cat("Weights directory: ", weights_dir, "\n")
