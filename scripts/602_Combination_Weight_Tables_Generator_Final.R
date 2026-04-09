# ============================================================
# 602_Combination_Weights_MR_FIXED_v2_EN.R
# Module combination weight table generator (English version)
# Two-stage strategy:
# Stage 1 (this script): Loose threshold (CLUMP_P = 1e-5) to avoid killing small combinations
# Stage 2 (604 script): Strict screening (5e-8 + full MR metrics)
# ============================================================

library(rprojroot)
library(TwoSampleMR)
library(data.table)
library(dplyr)

project_root <- find_rstudio_root_file()
data_dir <- file.path(project_root, "data")

# ============================================================
# Part 1: Generate combination lists per module
# ============================================================
cat("Scanning data/ folder for Target Genes files...\n")
txt_files <- list.files(data_dir,
                        pattern = "^Target Genes .*\\.txt$",
                        full.names = TRUE)

groups <- list()
for (f in txt_files) {
  lines <- readLines(f, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[lines != ""]
  
  if (length(lines) <= 1) {
    cat(" Skipping file (only 1 line or empty):", basename(f), "\n")
    next
  }
  
  core <- lines[1]                    # First line = core gene
  extended <- lines[-1]               # Remaining lines = extended genes
  
  # Extract module name from filename (more robust)
  module_name <- gsub("^Target Genes |\\.txt$", "", basename(f))
  
  groups[[length(groups) + 1]] <- list(
    core = core,
    extended = extended,
    module_name = module_name
  )
  
  cat(" ✓ Loaded module:", module_name, " | Core:", core,
      " | Extended genes:", length(extended), "\n")
}
cat("Total loaded", length(groups), "independent modules (GZMB / LAMP3 / NKG7 / SLC39A8)\n\n")

# ==================== Generate independent list files per module ====================
output_dir <- file.path(data_dir, "Exposure_IDs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (group in groups) {
  core <- group$core
  extended <- group$extended
  module_name <- group$module_name
  
  module_list_file <- file.path(output_dir, paste0("Combo_Weights_List_", module_name, ".txt"))
  file.create(module_list_file)
  
  m <- length(extended)
  comb_count <- 0
  
  for (n in 1:m) {                                 # Start from 1 (no core-only)
    combs <- combn(extended, n, simplify = FALSE)
    for (comb in combs) {
      line <- paste(c(core, comb), collapse = ",")
      cat(line, "\n", sep = "", file = module_list_file, append = TRUE)
      comb_count <- comb_count + 1
    }
  }
  
  cat(" ✓ Module", module_name, "→", comb_count,
      "combinations saved to", basename(module_list_file), "\n")
}

cat("\n=== Part 1 completed! Independent combination lists generated per module ===\n")
cat("(Modules are completely separate - GZMB will never mix with LAMP3)\n\n")

cat("\nMerging all module lists into one master file...\n")
master_list_file <- file.path(output_dir, "Combo_Weights_List.txt")
file.create(master_list_file)

total_combs <- 0
for (group in groups) {
  module_name <- group$module_name
  module_list_file <- file.path(output_dir, paste0("Combo_Weights_List_", module_name, ".txt"))
  
  if (!file.exists(module_list_file)) {
    cat(" Warning: missing", basename(module_list_file), "→ skipped\n")
    next
  }
  
  lines <- readLines(module_list_file)
  valid_lines <- lines[trimws(lines) != ""]
  
  if (length(valid_lines) == 0) {
    cat(" Warning:", basename(module_list_file), "is empty → skipped\n")
    next
  }
  
  cat(valid_lines, file = master_list_file, sep = "\n", append = TRUE)
  
  cat(" Merged:", basename(module_list_file), "→", length(valid_lines), "lines\n")
  total_combs <- total_combs + length(valid_lines)
}

cat("\n=== Merge completed ===\n")
cat("Total combinations in Combo_Weights_List.txt :", total_combs, "\n")
cat("Master file path:", master_list_file, "\n")
cat("(Each combination still comes from only one module - no cross-module mixing)\n\n")

# ============================================================
# Part 1.5: Generate UVMR_exposure_Combo_Weights_ids.txt
# Copy Combo_Weights_List.txt → data/UVMR_exposure_Combo_Weights_ids.txt
# Transform: "GZMB,TUBA1A,TUBA1B" → "GZMB_TUBA1A_TUBA1B_weights"
# ============================================================
combo_lines <- readLines(master_list_file)
combo_lines <- trimws(combo_lines)
combo_lines <- combo_lines[combo_lines != ""]

uvmr_ids <- gsub(",", "_", combo_lines)       # comma → underscore
uvmr_ids <- paste0(uvmr_ids, "_weights")      # append _weights

uvmr_file <- file.path(data_dir, "UVMR_exposure_Combo_Weights_ids.txt")
writeLines(uvmr_ids, uvmr_file)

cat("=== UVMR exposure ID file generated ===\n")
cat("Lines:", length(uvmr_ids), "\n")
cat("File:", uvmr_file, "\n")
cat("Preview:\n")
cat(paste("  ", head(uvmr_ids, 5)), sep = "\n")
cat("\n\n")

# ============================================================
# Part 2: Process combinations → Generate weight tables
# ============================================================
list_file <- file.path(data_dir, "Exposure_IDs", "Combo_Weights_List.txt")
input_dir  <- file.path(data_dir, "Exposure_IDs")
output_dir <- file.path(data_dir, "Exposure_IDs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==================== Custom clumping function (defensive) ====================
clump_snps <- function(df, clump_r2 = 0.001, clump_kb = 10000, clump_p = 1e-5) {
  df <- df[!is.na(df$chr.exposure) & !is.na(df$pos.exposure), ]
  if (nrow(df) == 0) {
    warning("No valid rows after removing NA in chr.exposure or pos.exposure.")
    return(df)
  }
  
  df$chr.exposure <- gsub("^chr", "", df$chr.exposure)
  
  if (!"pval.exposure" %in% colnames(df)) stop("Column 'pval.exposure' not found.")
  if (!is.numeric(df$pval.exposure)) {
    warning("Converting pval.exposure to numeric...")
    df$pval.exposure <- as.numeric(as.character(df$pval.exposure))
    df <- df[!is.na(df$pval.exposure), ]
  }
  
  tryCatch({
    clumped <- clump_data(df,
                          clump_r2 = clump_r2,
                          clump_kb = clump_kb,
                          clump_p1 = clump_p,
                          clump_p2 = clump_p,
                          pop = "EUR")
    return(clumped)
  }, error = function(e) {
    warning("LD clumping failed: ", conditionMessage(e), "\nReturning original data.")
    return(df)
  })
}

# ==================== Main processing loop ====================
combo_list <- readLines(list_file)
for (line in combo_list) {
  if (trimws(line) == "") next
  
  genes <- unlist(strsplit(line, ","))
  combo_name <- paste(genes, collapse = "_")
  output_file <- file.path(output_dir, paste0(combo_name, "_weights.csv"))
  
  df_list <- list()
  for (gene in genes) {
    input_file <- file.path(input_dir, paste0(gene, ".csv"))
    if (file.exists(input_file)) {
      df <- fread(input_file, header = TRUE, sep = ",")
      df$exposure <- combo_name
      
      required_cols <- c("SNP", "beta", "se", "pval", "effect_allele", "chr", "pos")
      missing_cols <- setdiff(required_cols, colnames(df))
      if (length(missing_cols) > 0) {
        warning(paste("Missing columns in", input_file, ":", paste(missing_cols, collapse = ", ")))
        next
      }
      
      cat("Loaded file:", input_file, "\n")
      df_list[[gene]] <- df
    } else {
      warning(paste("File does not exist:", input_file))
    }
  }
  
  if (length(df_list) == 0) {
    warning(paste("No data for combination:", combo_name))
    next
  }
  
  combined <- rbindlist(df_list, fill = TRUE)
  
  # === Keep SNP with smallest p-value (MR best practice) ===
  setDT(combined)
  combined <- combined[, .SD[which.min(pval)], by = SNP]
  combined <- as.data.frame(combined)
  
  if (nrow(combined) == 0) {
    warning(paste("No SNPs after deduplication for combination:", combo_name))
    next
  }
  
  # Allele cleaning
  combined$effect_allele[combined$effect_allele == ""] <- NA
  combined$other_allele[combined$other_allele == ""] <- NA
  valid_alleles <- c("A", "T", "C", "G", "AT", "TC", "GC", "CCT")
  combined$effect_allele[!combined$effect_allele %in% valid_alleles] <- NA
  
  cat("Processing combination:", combo_name, "\n")
  
  tryCatch({
    exp_dat <- format_data(
      combined,
      type = "exposure",
      snp_col = "SNP",
      beta_col = "beta",
      se_col = "se",
      eaf_col = "EAF",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "pval",
      chr_col = "chr",
      pos_col = "pos",
      phenotype_col = "exposure"
    )
    
    clumped <- clump_snps(exp_dat, clump_r2 = 0.001, clump_kb = 10000, clump_p = 1e-5)
    cat("Number of SNPs after clumping:", nrow(clumped), "\n")
    
    if (nrow(clumped) == 0) {
      warning(paste("No SNPs after clumping for combination:", combo_name))
      next
    }
    
    # Add F-statistic
    if ("beta.exposure" %in% colnames(clumped) && "se.exposure" %in% colnames(clumped)) {
      clumped$F <- (clumped$beta.exposure / clumped$se.exposure)^2
      clumped$F[is.na(clumped$F) | is.infinite(clumped$F)] <- NA
    }
    
    # Rename to clean format
    if (all(c("SNP", "beta.exposure", "se.exposure", "pval.exposure",
              "effect_allele.exposure", "other_allele.exposure",
              "eaf.exposure", "F") %in% colnames(clumped))) {
      clumped <- clumped %>%
        dplyr::rename(
          beta = beta.exposure,
          se = se.exposure,
          pval = pval.exposure,
          effect_allele = effect_allele.exposure,
          other_allele = other_allele.exposure,
          EAF = eaf.exposure
        )
    }
    
    fwrite(clumped, output_file, sep = ",", quote = TRUE)
    cat("Weight table saved to:", output_file, "| IVs:", nrow(clumped), "\n")
    
  }, error = function(e) {
    warning("Error processing combination ", combo_name, ": ", conditionMessage(e))
  })
}

cat("\n=== Stage 1 completed! All weight tables generated ===\n")
cat("Loose threshold used: CLUMP_P = 1e-5\n")
cat("Output folder:", output_dir, "\n")
cat("Next step: Run your 604 script with strict 5e-8 filtering!\n")