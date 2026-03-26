library(rprojroot)
project_root <- find_rstudio_root_file()
data_dir <- file.path(project_root, "data")
# Adjustable parameters (concentrated at the beginning of the code)
# Base working directory path
BASE_DIR <- file.path(data_dir, "GWAS")
# Gene name list file path
GENE_LIST_FILE <- file.path(data_dir, "Target Genes.txt")
# P-value threshold for filtering significant SNPs
PVAL_THRESHOLD <- 1e-5
# LD clumping distance threshold (unit: bp)
DIST_THRESHOLD <- 10000
# Load required packages
library(data.table)
library(dplyr)
library(uuid)
library(biomaRt) # for mapping rsID
# Check if gene list file exists
if (!file.exists(GENE_LIST_FILE)) {
  stop("Gene list file does not exist: ", GENE_LIST_FILE)
}
# Read gene name list and trim leading/trailing spaces
gene_list <- read.table(GENE_LIST_FILE, header = FALSE, stringsAsFactors = FALSE, encoding = "UTF-8")$V1
gene_list <- trimws(gene_list)
# Function: calculate standard error (SE) from confidence interval
calculate_se <- function(ci_lower, ci_upper) {
  if (is.na(ci_lower) || is.na(ci_upper)) return(NA)
  (ci_upper - ci_lower) / (2 * 1.96)
}
# Function: calculate F-statistic
calculate_f <- function(beta, se) {
  if (is.na(beta) || is.na(se) || se == 0) return(NA)
  (beta / se)^2
}
# Function: extract effect allele and other allele from Variant ID
extract_alleles <- function(variant_id, beta) {
  alleles <- tryCatch({
    parts <- unlist(strsplit(as.character(variant_id), "_"))
    if (length(parts) >= 4) parts[3:4] else c(NA, NA)
  }, error = function(e) c(NA, NA))
  if (is.na(beta) || length(alleles) < 2 || any(is.na(alleles))) return(c(NA, NA))
  if (beta >= 0) {
    effect_allele <- alleles[2]
    other_allele <- alleles[1]
  } else {
    effect_allele <- alleles[1]
    other_allele <- alleles[2]
  }
  return(c(effect_allele, other_allele))
}
# Function: map rsID to chr and pos using biomaRt
get_snp_pos <- function(rsids) {
  ensembl <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
  snp_info <- tryCatch({
    getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
          filters = "snp_filter", values = rsids[!is.na(rsids) & rsids != ""],
          mart = ensembl)
  }, error = function(e) {
    message("biomaRt query failed: ", e$message)
    return(data.frame(refsnp_id = character(), chr_name = character(), chrom_start = numeric()))
  })
  return(snp_info)
}
# LD clumping function: based on distance
clump_snps <- function(df, dist_threshold = DIST_THRESHOLD) {
  # Remove rows with NA in chr or pos
  df <- df[!is.na(df$chr) & !is.na(df$pos), ]
  if (nrow(df) == 0) return(df)
  # Sort by pval ascending
  df <- df[order(df$pval), ]
  kept <- c()
  for (i in 1:nrow(df)) {
    if (!any(kept)) {
      kept <- c(kept, i)
      next
    }
    is_indep <- TRUE
    for (j in kept) {
      if (df$chr[i] == df$chr[j] && abs(df$pos[i] - df$pos[j]) < dist_threshold) {
        is_indep <- FALSE
        break
      }
    }
    if (is_indep) kept <- c(kept, i)
  }
  return(df[kept, ])
}
# Traverse gene list
for (gene in gene_list) {
  cat("Processing gene: ", gene, "\n")
  gene_dir <- file.path(BASE_DIR, gene)
  if (!dir.exists(gene_dir)) {
    cat("Directory does not exist, skipping: ", gene_dir, "\n")
    next
  }
  setwd(gene_dir)
  
  # Initialize gene-specific result data frame, add chr and pos columns
  result <- data.frame(
    SNP = character(),
    beta = numeric(),
    se = numeric(),
    pval = numeric(),
    effect_allele = character(),
    other_allele = character(),
    exposure = character(),
    EAF = numeric(),
    F = numeric(),
    chr = character(),
    pos = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 1. Process gwas-association-ensemblMappedGenes.tsv
  gwas_catalog_file <- paste0(gene, " gwas-association-ensemblMappedGenes.tsv")
  if (file.exists(gwas_catalog_file)) {
    try({
      gwas_catalog <- fread(gwas_catalog_file, header = TRUE, stringsAsFactors = FALSE, encoding = "UTF-8")
      cat("Gene ", gene, " GWAS catalog data preview:\n")
      print(head(gwas_catalog))
      gwas_catalog$SNPS <- as.character(gwas_catalog$SNPS)
      gwas_catalog$`STRONGEST SNP-RISK ALLELE` <- as.character(gwas_catalog$`STRONGEST SNP-RISK ALLELE`)
      gwas_catalog$SNPS <- sub("-.*", "", gwas_catalog$SNPS)
      gwas_catalog$`STRONGEST SNP-RISK ALLELE` <- sub("-\\?$", "", gwas_catalog$`STRONGEST SNP-RISK ALLELE`)
      gwas_catalog$`P-VALUE` <- as.numeric(gwas_catalog$`P-VALUE`)
      gwas_catalog <- gwas_catalog[!is.na(gwas_catalog$`P-VALUE`) & gwas_catalog$`P-VALUE` < PVAL_THRESHOLD, ]
      for (i in 1:nrow(gwas_catalog)) {
        row <- gwas_catalog[i, ]
        snp <- as.character(row$SNPS)
        pval <- as.numeric(row$`P-VALUE`)
        beta <- tryCatch({
          as.numeric(row$`OR or BETA`)
        }, warning = function(w) NA, error = function(e) NA)
        ci_text <- as.character(row$`95% CI (TEXT)`)
        ci_values <- tryCatch({
          nums <- as.numeric(unlist(regmatches(ci_text, gregexpr("[0-9.]+", ci_text))))
          if (length(nums) >= 2) nums[1:2] else c(NA, NA)
        }, warning = function(w) c(NA, NA), error = function(e) c(NA, NA))
        se <- calculate_se(ci_values[1], ci_values[2])
        if (is.na(beta) || is.na(se)) {
          cat("Skip gene ", gene, " row ", i, ": invalid beta or se\n")
          next
        }
        f_stat <- calculate_f(beta, se)
        effect_allele <- tryCatch({
          sub(".*-", "", as.character(row$`STRONGEST SNP-RISK ALLELE`))
        }, error = function(e) NA)
        other_allele <- tryCatch({
          alleles <- unlist(strsplit(as.character(row$SNPS), "_"))[3:4]
          if (length(alleles) >= 2 && !is.na(effect_allele) && effect_allele != alleles[1]) {
            alleles[1]
          } else {
            NA
          }
        }, error = function(e) NA)
        eaf <- tryCatch({
          as.numeric(row$`RISK ALLELE FREQUENCY`)
        }, warning = function(w) NA, error = function(e) NA)
        chr <- as.character(row$CHR_ID)
        pos <- as.numeric(row$CHR_POS)
        result <- rbind(result, data.frame(
          SNP = as.character(snp),
          beta = as.numeric(beta),
          se = as.numeric(se),
          pval = as.numeric(pval),
          effect_allele = as.character(effect_allele),
          other_allele = as.character(other_allele),
          exposure = as.character(gene),
          EAF = as.numeric(eaf),
          F = as.numeric(f_stat),
          chr = as.character(chr),
          pos = as.numeric(pos),
          stringsAsFactors = FALSE
        ))
      }
    }, silent = TRUE)
  }
  
  # 2. Process eQTLs GTEx Portal.csv
  eqtl_file <- paste0(gene, " eQTLs GTEx Portal.csv")
  if (file.exists(eqtl_file)) {
    try({
      eqtl_data <- fread(eqtl_file, header = TRUE, stringsAsFactors = FALSE, encoding = "UTF-8")
      cat("Gene ", gene, " eQTL data preview:\n")
      print(head(eqtl_data))
      eqtl_data$`Variant Id` <- as.character(eqtl_data$`Variant Id`)
      eqtl_data$`SNP Id` <- as.character(eqtl_data$`SNP Id`)
      eqtl_data$`SNP Id`[is.na(eqtl_data$`SNP Id`)] <- paste0("NA_", UUIDgenerate())
      eqtl_data$`P-Value` <- as.numeric(eqtl_data$`P-Value`)
      eqtl_data$NES <- as.numeric(eqtl_data$NES)
      eqtl_data <- eqtl_data[!is.na(eqtl_data$`P-Value`) & eqtl_data$`P-Value` < PVAL_THRESHOLD & !is.na(eqtl_data$NES), ]
      for (i in 1:nrow(eqtl_data)) {
        row <- eqtl_data[i, ]
        snp <- as.character(row$`SNP Id`)
        pval <- as.numeric(row$`P-Value`)
        beta <- as.numeric(row$NES)
        z <- tryCatch({
          sqrt(2 * log(1 / pval))
        }, error = function(e) NA)
        se <- if (!is.na(z) && z != 0) abs(beta) / z else NA
        if (is.na(beta) || is.na(se)) {
          cat("Skip gene ", gene, " row ", i, ": invalid beta or se\n")
          next
        }
        f_stat <- calculate_f(beta, se)
        variant_id <- as.character(row$`Variant Id`)
        alleles <- tryCatch({
          parts <- unlist(strsplit(variant_id, "_"))
          if (length(parts) >= 4) parts[3:4] else c(NA, NA)
        }, error = function(e) {
          cat("Parse gene ", gene, " Variant Id error: ", variant_id, "\n")
          c(NA, NA)
        })
        eaf <- NA
        chr <- tryCatch({
          unlist(strsplit(variant_id, "_"))[1]
        }, error = function(e) NA)
        pos <- tryCatch({
          as.numeric(unlist(strsplit(variant_id, "_")))
        }, error = function(e) NA)
        result <- rbind(result, data.frame(
          SNP = as.character(snp),
          beta = as.numeric(beta),
          se = as.numeric(se),
          pval = as.numeric(pval),
          effect_allele = as.character(alleles[2]),
          other_allele = as.character(alleles[1]),
          exposure = as.character(gene),
          EAF = as.numeric(eaf),
          F = as.numeric(f_stat),
          chr = as.character(chr),
          pos = as.numeric(pos),
          stringsAsFactors = FALSE
        ))
      }
    }, silent = TRUE)
  }
  
  # 3. Process sQTLs GTEx Portal.csv
  sqtl_file <- paste0(gene, " sQTLs GTEx Portal.csv")
  if (file.exists(sqtl_file)) {
    try({
      sqtl_data <- fread(sqtl_file, header = TRUE, stringsAsFactors = FALSE, encoding = "UTF-8")
      cat("Gene ", gene, " sQTL data preview:\n")
      print(head(sqtl_data))
      sqtl_data$`Variant Id` <- as.character(sqtl_data$`Variant Id`)
      sqtl_data$`SNP Id` <- as.character(sqtl_data$`SNP Id`)
      sqtl_data$`SNP Id`[is.na(sqtl_data$`SNP Id`)] <- paste0("NA_", UUIDgenerate())
      sqtl_data$`P-Value` <- as.numeric(sqtl_data$`P-Value`)
      sqtl_data$NES <- as.numeric(sqtl_data$NES)
      sqtl_data <- sqtl_data[!is.na(sqtl_data$`P-Value`) & sqtl_data$`P-Value` < PVAL_THRESHOLD & !is.na(sqtl_data$NES), ]
      for (i in 1:nrow(sqtl_data)) {
        row <- sqtl_data[i, ]
        snp <- as.character(row$`SNP Id`)
        pval <- as.numeric(row$`P-Value`)
        beta <- as.numeric(row$NES)
        z <- tryCatch({
          sqrt(2 * log(1 / pval))
        }, error = function(e) NA)
        se <- if (!is.na(z) && z != 0) abs(beta) / z else NA
        if (is.na(beta) || is.na(se)) {
          cat("Skip gene ", gene, " row ", i, ": invalid beta or se\n")
          next
        }
        f_stat <- calculate_f(beta, se)
        variant_id <- as.character(row$`Variant Id`)
        alleles <- tryCatch({
          parts <- unlist(strsplit(variant_id, "_"))
          if (length(parts) >= 4) parts[3:4] else c(NA, NA)
        }, error = function(e) {
          cat("Parse gene ", gene, " Variant Id error: ", variant_id, "\n")
          c(NA, NA)
        })
        eaf <- NA
        chr <- tryCatch({
          unlist(strsplit(variant_id, "_"))[1]
        }, error = function(e) NA)
        pos <- tryCatch({
          as.numeric(unlist(strsplit(variant_id, "_"))[2]) # fix bracket
        }, error = function(e) NA)
        result <- rbind(result, data.frame(
          SNP = as.character(snp),
          beta = as.numeric(beta),
          se = as.numeric(se),
          pval = as.numeric(pval),
          effect_allele = as.character(alleles[2]),
          other_allele = as.character(alleles[1]),
          exposure = as.character(gene),
          EAF = as.numeric(eaf),
          F = as.numeric(f_stat),
          chr = as.character(chr),
          pos = as.numeric(pos),
          stringsAsFactors = FALSE
        ))
      }
    }, silent = TRUE)
  }
  
  # 4. Process top_associations.csv
  top_assoc_file <- paste0(gene, "_top_associations.csv")
  if (file.exists(top_assoc_file)) {
    try({
      top_assoc <- fread(top_assoc_file, header = TRUE, stringsAsFactors = FALSE, encoding = "UTF-8")
      cat("Gene ", gene, " top_assoc data preview:\n")
      print(head(top_assoc))
      top_assoc$rsids <- as.character(top_assoc$rsids)
      top_assoc$rsids <- sapply(strsplit(top_assoc$rsids, ","), function(x) x[1])
      top_assoc$rsids[is.na(top_assoc$rsids) | top_assoc$rsids == "NA"] <- paste0("NA_", UUIDgenerate())
      top_assoc$pval <- as.numeric(top_assoc$pval)
      top_assoc$beta <- as.numeric(top_assoc$beta)
      top_assoc$sebeta <- as.numeric(top_assoc$sebeta)
      top_assoc <- top_assoc[!is.na(top_assoc$pval) & top_assoc$pval < PVAL_THRESHOLD, ]
      
      # Initialize top_assoc specific result data frame
      top_assoc_result <- data.frame(
        SNP = character(),
        beta = numeric(),
        se = numeric(),
        pval = numeric(),
        effect_allele = character(),
        other_allele = character(),
        exposure = character(),
        EAF = numeric(),
        F = numeric(),
        chr = character(),
        pos = numeric(),
        stringsAsFactors = FALSE
      )
      
      # Use biomaRt to map rsID to chr and pos
      valid_rsids <- top_assoc$rsids[!grepl("^NA_", top_assoc$rsids)]
      if (length(valid_rsids) > 0) {
        snp_info <- get_snp_pos(valid_rsids)
        top_assoc$chr <- NA
        top_assoc$pos <- NA
        for (k in 1:nrow(top_assoc)) {
          if (!grepl("^NA_", top_assoc$rsids[k])) {
            match <- snp_info[snp_info$refsnp_id == top_assoc$rsids[k], ]
            if (nrow(match) > 0) {
              top_assoc$chr[k] <- match$chr_name[1]
              top_assoc$pos[k] <- match$chrom_start[1]
            }
          }
        }
      } else {
        top_assoc$chr <- NA
        top_assoc$pos <- NA
      }
      
      # Fill top_assoc_result and result data frame
      for (i in 1:nrow(top_assoc)) {
        row <- top_assoc[i, ]
        snp <- as.character(row$rsids)
        pval <- as.numeric(row$pval)
        beta <- as.numeric(row$beta)
        se <- as.numeric(row$sebeta)
        if (is.na(beta) || is.na(se)) {
          cat("Skip gene ", gene, " row ", i, ": invalid beta or se\n")
          next
        }
        f_stat <- calculate_f(beta, se)
        chr <- as.character(row$chr)
        pos <- as.numeric(row$pos)
        
        # Add to top_assoc_result
        top_assoc_result <- rbind(top_assoc_result, data.frame(
          SNP = as.character(snp),
          beta = as.numeric(beta),
          se = as.numeric(se),
          pval = as.numeric(pval),
          effect_allele = "",
          other_allele = "",
          exposure = as.character(gene),
          EAF = NA,
          F = as.numeric(f_stat),
          chr = as.character(chr),
          pos = as.numeric(pos),
          stringsAsFactors = FALSE
        ))
        
        # Add to result (original logic)
        result <- rbind(result, data.frame(
          SNP = as.character(snp),
          beta = as.numeric(beta),
          se = as.numeric(se),
          pval = as.numeric(pval),
          effect_allele = "",
          other_allele = "",
          exposure = as.character(gene),
          EAF = NA,
          F = as.numeric(f_stat),
          chr = as.character(chr),
          pos = as.numeric(pos),
          stringsAsFactors = FALSE
        ))
      }
      
      # Deduplicate and perform LD clumping on top_assoc_result
      if (nrow(top_assoc_result) > 0) {
        setDT(top_assoc_result)
        top_assoc_result <- top_assoc_result[!is.na(pval), .SD[pval == min(pval, na.rm = TRUE)], by = SNP]
        top_assoc_result$EAF <- as.numeric(top_assoc_result$EAF)
        
        # Apply LD clumping
        top_assoc_clumped <- clump_snps(top_assoc_result, dist_threshold = DIST_THRESHOLD)
        
        # Output to top_associations_LD.csv
        top_assoc_output_file <- paste0(gene, "_top_associations_LD.csv")
        write.csv(top_assoc_clumped, file = top_assoc_output_file, row.names = FALSE, quote = TRUE)
        cat("Output saved (top_assoc LD clumped): ", top_assoc_output_file, "\n")
      }
    }, silent = TRUE)
  }
  
  # Deduplicate, keep record with lowest p-value
  if (nrow(result) > 0) {
    cat("Gene ", gene, " result data structure:\n")
    str(result)
    setDT(result)
    result <- result[!is.na(pval), .SD[pval == min(pval, na.rm = TRUE)], by = SNP]
    result$EAF <- as.numeric(result$EAF)
    
    # Apply clumping
    result_clumped <- clump_snps(result, dist_threshold = DIST_THRESHOLD)
    
    # Output
    output_file <- paste0(gene, "_UVMR_SNPs_clumped.csv")
    write.csv(result_clumped, file = output_file, row.names = FALSE, quote = TRUE)
    cat("Output saved (LD clumped): ", output_file, "\n")
  } else {
    cat("No valid SNPs found for gene ", gene, "\n")
  }
}
cat("Analysis completed.\n")

# =====================================================================
# New script: Copy and rename *_UVMR_SNPs_clumped.csv files
# =====================================================================
cat("\n=== Start copying and renaming _UVMR_SNPs_clumped.csv files ===\n")

source_base <- file.path(data_dir, "GWAS")
target_dir  <- file.path(data_dir, "Exposure_IDs")

if (!dir.exists(target_dir)) {
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Target directory created: ", target_dir, "\n")
}

gene_folders <- list.dirs(source_base, full.names = TRUE, recursive = FALSE)

success <- 0
skipped <- 0

for (folder in gene_folders) {
  gene_name <- basename(folder)
  source_file <- file.path(folder, paste0(gene_name, "_UVMR_SNPs_clumped.csv"))
  
  if (!file.exists(source_file)) {
    cat("Skipped (file not found): ", gene_name, "\n")
    skipped <- skipped + 1
    next
  }
  
  final_target <- file.path(target_dir, paste0(gene_name, ".csv"))
  
  if (file.exists(final_target)) file.remove(final_target)
  
  file.copy(source_file, final_target, overwrite = TRUE)
  
  cat("Success: ", gene_name, " → ", basename(final_target), "\n")
  success <- success + 1
}

cat("\nCopy and rename completed! Success: ", success, " files, Skipped: ", skipped, " files\n")
cat("Target directory: ", target_dir, "\n")