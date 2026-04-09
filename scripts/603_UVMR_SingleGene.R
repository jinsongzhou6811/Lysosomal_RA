# ========================================================
# Standard Path Management Template (All Modules)
# ========================================================
# install.packages("rprojroot") # Only need to execute once
library(rprojroot)
# Automatically detect project root directory (effective on any computer, any drive)
project_root <- find_rstudio_root_file()
# Fixed paths (shared by all modules)
data_dir <- file.path(project_root, "data")
outputs_dir <- file.path(project_root, "outputs")
# ------------------- Only need to modify this line -------------------
module_name <- "06_UVMR_Analysis" # ←←← Change here!!
# ------------------------------------------------------
# Module actual directory (automatically created)
module_dir <- file.path(outputs_dir, module_name)
dir.create(module_dir, showWarnings = FALSE, recursive = TRUE)
# ========================================================
# Configuration parameters (centralized for easy modification)
config <- list(
  base_dir = data_dir,
  output_base_dir = file.path(module_dir, "02_UVMR_Single_Gene"),
  exposure_ids_file = file.path(data_dir, "UVMR_exposure_Single_gene_ids.txt"), # Use absolute path
  outcome_ids_file = file.path(data_dir, "UVMR_outcome_local_GWAS_ids.txt"),
  num_cores = 20,
  batch_size = 14,
  F_MIN = 10,
  mr_presso_params = list(
    NbDistribution = 1000,
    SignifThreshold = 0.1
  )
)

# Install and load necessary R packages
library(foreach)
library(doParallel)
library(TwoSampleMR)
library(data.table)
library(MRPRESSO)
library(ggplot2)
library(logger)
library(vcfR)
library(dplyr)

# Set up logging
log_dir <- file.path(config$output_base_dir, "logs")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
global_log_file <- file.path(log_dir, paste0("mr_analysis_global_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_appender(appender_file(global_log_file))
log_threshold(INFO)
log_formatter(formatter_sprintf) # Use sprintf formatter to avoid glue issues
log_info("Base directory set to: %s", config$base_dir)
log_info("Output base directory set to: %s", config$output_base_dir)
# Read exposure and outcome IDs
if (!file.exists(config$exposure_ids_file)) {
  log_error("Exposure ID file does not exist: %s", config$exposure_ids_file)
  stop("Exposure ID file does not exist: ", config$exposure_ids_file)
}
# Read exposure IDs and thoroughly clean
exposure_ids <- readLines(config$exposure_ids_file)
exposure_ids <- trimws(exposure_ids) # Remove leading and trailing spaces
exposure_ids <- exposure_ids[nzchar(exposure_ids)] # Remove empty lines
exposure_ids <- na.omit(exposure_ids) # Remove NA
exposure_ids <- as.character(exposure_ids) # Force convert to character
exposure_ids <- unique(exposure_ids) # Remove duplicates (for safety)
log_info("Number of exposure_ids after final cleanup: %d", length(exposure_ids))
log_info("Read %d exposure IDs from %s: %s", length(exposure_ids), config$exposure_ids_file, paste(exposure_ids, collapse=', '))
if (!file.exists(config$outcome_ids_file)) {
  log_error("Outcome ID file does not exist: %s", config$outcome_ids_file)
  stop("Outcome ID file does not exist: ", config$outcome_ids_file)
}
outcome_ids <- readLines(config$outcome_ids_file)
outcome_ids <- outcome_ids[nzchar(outcome_ids)]
log_info("Read %d outcome IDs from %s: %s", length(outcome_ids), config$outcome_ids_file, paste(outcome_ids, collapse=', '))
# Set up parallel computing with the number of cores
registerDoParallel(cores = config$num_cores)
log_info("Parallel computing enabled, using %d cores", config$num_cores)
# Function: Read and convert VCF file (unchanged)
read_and_convert_vcf <- function(vcf_file, outcome_id) {
  tryCatch({
    vcf <- read.vcfR(vcf_file, verbose = FALSE)
    fixed <- vcf@fix
    fixed_df <- data.frame(
      CHROM = fixed[, "CHROM"],
      POS = fixed[, "POS"],
      SNP = fixed[, "ID"],
      REF = fixed[, "REF"],
      ALT = fixed[, "ALT"],
      stringsAsFactors = FALSE
    )
    gt <- vcf@gt
    format_fields <- strsplit(gt[, "FORMAT"], ":")[[1]]
    log_info("FORMAT fields in VCF file %s: %s", vcf_file, paste(format_fields, collapse=', '))
    
    required_fields <- c("ES", "SE", "LP")
    if (!all(required_fields %in% format_fields)) {
      missing_fields <- paste(required_fields[!required_fields %in% format_fields], collapse = ', ')
      log_error("VCF file %s missing required FORMAT fields: %s", vcf_file, missing_fields)
      stop("VCF file missing required FORMAT fields: ", missing_fields)
    }
    
    gt_data <- gt[, 2]
    gt_split <- do.call(rbind, lapply(gt_data, function(x) {
      split_vals <- strsplit(x, ":")[[1]]
      if (length(split_vals) != length(format_fields)) {
        return(rep(NA, length(format_fields)))
      }
      split_vals
    }))
    if (ncol(gt_split) != length(format_fields)) {
      log_error("VCF file %s FORMAT fields do not match data column count", vcf_file)
      stop("VCF file FORMAT fields do not match data column count")
    }
    gt_df <- as.data.frame(gt_split, stringsAsFactors = FALSE)
    colnames(gt_df) <- format_fields
    
    gt_df$ES <- as.numeric(as.character(gt_df$ES))
    gt_df$SE <- as.numeric(as.character(gt_df$SE))
    gt_df$LP <- as.numeric(as.character(gt_df$LP))
    
    if ("AF" %in% format_fields) {
      gt_df$AF <- as.numeric(as.character(gt_df$AF))
    } else {
      log_info("VCF file %s missing AF field, default filling eaf as 0.5", vcf_file)
      gt_df$AF <- 0.5
    }
    
    gt_df$pval <- 10^(-gt_df$LP)
    
    invalid_rows <- is.na(gt_df$ES) | is.na(gt_df$SE) | is.na(gt_df$LP) | !is.finite(gt_df$ES) | !is.finite(gt_df$SE)
    if (any(invalid_rows)) {
      log_info("Removed %d invalid records (containing NA or non-finite values) from VCF file %s", sum(invalid_rows), vcf_file)
      fixed_df <- fixed_df[!invalid_rows, ]
      gt_df <- gt_df[!invalid_rows, ]
    }
    
    outcome_df <- cbind(fixed_df, gt_df)
    outcome_df <- outcome_df %>%
      dplyr::select(
        SNP = SNP,
        beta = ES,
        se = SE,
        pval = pval,
        effect_allele = ALT,
        other_allele = REF,
        eaf = AF
      ) %>%
      mutate(outcome = outcome_id)
    return(outcome_df)
  }, error = function(e) {
    log_error("Failed to read VCF file %s: %s", vcf_file, e$message)
    stop("Failed to read VCF file: ", e$message)
  })
}
# Function: Process exposure-outcome pair (optimized: uniform SNP threshold <3 as insufficient; ensure F.exposure retained)
process_exposure_outcome <- function(exposure_ID, outcome_ID, exposure_dir, outcome_raw = NULL) {
  tryCatch({
    log_info("Starting analysis of exposure %s and outcome %s", exposure_ID, outcome_ID)
    
    # Create output directory
    output_dir <- file.path(exposure_dir, paste0("results_", exposure_ID, "_", outcome_ID))
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    log_info("Output directory created: %s", output_dir)
    
    # Read exposure data
    exposure_file <- file.path(config$base_dir, "Exposure_IDs", paste0(exposure_ID, ".csv"))
    if (!file.exists(exposure_file)) {
      log_error("Exposure data file does not exist: %s", exposure_file)
      stop("Exposure data file does not exist: ", exposure_file)
    }
    exposure_dat <- fread(exposure_file)
    required_columns <- c("SNP", "beta", "se", "pval", "effect_allele", "other_allele", "EAF", "F")
    if (!all(required_columns %in% colnames(exposure_dat))) {
      log_error("Exposure data file %s missing required columns", exposure_file)
      stop("Exposure data file missing required columns: ", exposure_file)
    }
    
    initial_snp_count <- nrow(exposure_dat)
   exposure_dat <- exposure_dat[!is.na(F) & F > config$F_MIN & is.finite(F), ]
    filtered_snp_count <- nrow(exposure_dat)
    if (filtered_snp_count == 0) {
      log_error("No remaining SNPs after filtering invalid F values: %s", exposure_file)
      stop("No remaining SNPs after filtering invalid F values")
    }
    if (initial_snp_count > filtered_snp_count) {
      log_info("Removed %d SNPs with invalid F values (NA or <= %d) from exposure %s", 
         initial_snp_count - filtered_snp_count, config$F_MIN, exposure_ID)
    }
    
    valid_alleles <- c("A", "C", "T", "G", "D", "I")
    invalid_allele_snps <- exposure_dat[!effect_allele %in% valid_alleles | !other_allele %in% valid_alleles | is.na(effect_allele) | is.na(other_allele), SNP]
    if (length(invalid_allele_snps) > 0) {
      log_info("Removed %d SNPs due to invalid effect_allele or other_allele: %s", length(invalid_allele_snps), paste(invalid_allele_snps, collapse=', '))
      exposure_dat <- exposure_dat[!SNP %in% invalid_allele_snps, ]
    }
    if (nrow(exposure_dat) == 0) {
      log_error("No remaining SNPs after filtering invalid alleles")
      stop("No remaining SNPs after filtering invalid alleles")
    }
    
    na_eaf_count <- sum(is.na(exposure_dat$EAF))
    exposure_dat$EAF <- ifelse(is.na(exposure_dat$EAF), 0.5, exposure_dat$EAF)
    if (na_eaf_count > 0) {
      log_info("Filled default value 0.5 for EAF of %d SNPs", na_eaf_count)
    }
    
    write.csv(exposure_dat, file.path(output_dir, "exposure_data_list.csv"), row.names = FALSE)
    log_info("Exposure data list saved: %s", file.path(output_dir, 'exposure_data_list.csv'))
    
    exposure_dat$exposure <- exposure_ID
    exposure_dat <- as.data.frame(exposure_dat)
    exposure_SNPs <- format_data(
      exposure_dat,
      type = "exposure",
      snps = NULL,
      header = TRUE,
      phenotype_col = "exposure",
      snp_col = "SNP",
      beta_col = "beta",
      se_col = "se",
      pval_col = "pval",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      eaf_col = "EAF"
    )
    exposure_SNPs <- merge(exposure_SNPs, exposure_dat[, c("SNP", "F")], by = "SNP", all.x = TRUE)
    colnames(exposure_SNPs)[colnames(exposure_SNPs) == "F"] <- "F.exposure"
    if (!"F.exposure" %in% colnames(exposure_SNPs)) {
      log_error("Unable to retain F.exposure column in exposure_SNPs")
      stop("Unable to retain F.exposure column in exposure_SNPs")
    }
    log_info("Extracted %d SNPs for exposure %s, F.exposure column retained", nrow(exposure_SNPs), exposure_ID)
    write.csv(exposure_SNPs, file.path(output_dir, "exposure_SNPs.csv"), row.names = FALSE)
    
    # Use pre-read outcome data or read new outcome data
    if (is.null(outcome_raw)) {
      outcome_file <- file.path(config$base_dir, "GWAS", paste0(outcome_ID, ".vcf.gz"))
      if (!file.exists(outcome_file)) {
        log_error("Outcome data file does not exist: %s", outcome_file)
        stop("Outcome data file does not exist: ", outcome_file)
      }
      outcome_raw <- read_and_convert_vcf(outcome_file, outcome_ID)
    }
    missing_snps <- setdiff(exposure_SNPs$SNP, outcome_raw$SNP)
    if (length(missing_snps) > 0) {
      log_info("Missing %d exposure SNPs in VCF file: %s", length(missing_snps), paste(missing_snps, collapse=', '))
    }
    outcome_SNPs <- format_data(
      outcome_raw,
      type = "outcome",
      snps = exposure_SNPs$SNP,
      header = TRUE,
      phenotype_col = "outcome",
      snp_col = "SNP",
      beta_col = "beta",
      se_col = "se",
      pval_col = "pval",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      eaf_col = "eaf"
    )
    log_info("Extracted %d SNPs for outcome %s", nrow(outcome_SNPs), outcome_ID)
    if (nrow(outcome_SNPs) == 0) {
      log_error("No matching exposure SNPs in outcome %s, unable to continue analysis", outcome_ID)
      return(list(
        outcome_id = outcome_ID,
        id.exposure = exposure_ID,
        status = "failed",
        error = "No matching SNPs found in outcome data",
        snps_after_filter = 0,
        f_exposure_min = NA,
        f_exposure_max = NA,
        f_exposure_mean = NA,
        passed_screening = FALSE
      ))
    }
    
    write.csv(outcome_SNPs, file.path(output_dir, "outcome_data_list.csv"), row.names = FALSE)
    log_info("Outcome data list saved: %s", file.path(output_dir, 'outcome_data_list.csv'))
    
    # Handle NA values in other_allele
    na_other_allele <- is.na(exposure_SNPs$other_allele.exposure)
    if (any(na_other_allele)) {
      matching_snps <- exposure_SNPs$SNP[na_other_allele]
      outcome_matches <- outcome_SNPs[outcome_SNPs$SNP %in% matching_snps, c("SNP", "other_allele.outcome")]
      if (nrow(outcome_matches) > 0) {
        exposure_SNPs$other_allele.exposure[na_other_allele] <- outcome_matches$other_allele.outcome[
          match(exposure_SNPs$SNP[na_other_allele], outcome_matches$SNP)
        ]
        log_info("Replaced %d NA values in other_allele for exposure %s using alleles from outcome %s", sum(!is.na(exposure_SNPs$other_allele.exposure[na_other_allele])), exposure_ID, outcome_ID)
      } else {
        log_info("No matching SNPs found in outcome data to replace NA values in other_allele for exposure %s", exposure_ID)
      }
    } else {
      log_info("No NA values in other_allele found in exposure %s", exposure_ID)
    }
    
    write.csv(exposure_SNPs, file.path(output_dir, "exposure_allele_replaced_list.csv"), row.names = FALSE)
    log_info("Exposure data list after replacing other_allele saved: %s", file.path(output_dir, 'exposure_allele_replaced_list.csv'))
    
    # Data harmonization (optimized: explicitly retain F.exposure)
    exposure_SNPs$units.exposure <- "unit"
    outcome_SNPs$units.outcome <- "unit"
    dat <- harmonise_data(exposure_SNPs, outcome_SNPs, action=2)
    if (nrow(dat) == 0) {
      log_error("No common SNPs after harmonization for exposure %s and outcome %s", exposure_ID, outcome_ID)
      return(list(
        outcome_id = outcome_ID,
        id.exposure = exposure_ID,
        status = "failed",
        error = "No common SNPs after harmonisation",
        snps_after_filter = 0,
        f_exposure_min = NA,
        f_exposure_max = NA,
        f_exposure_mean = NA,
        passed_screening = FALSE
      ))
    }
    
    log_info("Remaining %d SNPs after harmonization", nrow(dat))
    write.csv(dat, file.path(output_dir, "harmonised_data_before_filter.csv"), row.names = FALSE)
    log_info("Harmonized data saved: %s", file.path(output_dir, 'harmonised_data_before_filter.csv'))
    
    if ("F.exposure" %in% colnames(dat)) {
        log_info("F.exposure distribution: min=%.2f, max=%.2f, mean=%.2f", 
                 min(dat$F.exposure, na.rm=TRUE), 
                 max(dat$F.exposure, na.rm=TRUE), 
                 mean(dat$F.exposure, na.rm=TRUE))
        
        palindromic_snps <- dat$SNP[dat$palindromic & !dat$mr_keep]
        if (length(palindromic_snps) > 0) {
            log_info("Removed %d palindromic SNPs: %s", length(palindromic_snps), 
                     paste(palindromic_snps, collapse=', '))
        }
    } else {
      log_error("Missing F.exposure column in harmonized data")
      return(list(
        outcome_id = outcome_ID,
        id.exposure = exposure_ID,
        status = "failed",
        error = "Missing F.exposure column in harmonised data",
        snps_after_filter = 0,
        f_exposure_min = NA,
        f_exposure_max = NA,
        f_exposure_mean = NA,
        passed_screening = FALSE
      ))
    }
    
    dat_filtered <- dat[!is.na(dat$F.exposure) & is.finite(dat$F.exposure) & dat$F.exposure > config$F_MIN & dat$mr_keep, ]
    if (nrow(dat_filtered) == 0) {
      log_error("No remaining SNPs after F statistic filtering (F>%d) for exposure %s and outcome %s", config$F_MIN, exposure_ID, outcome_ID)
      write.csv(dat, file.path(output_dir, "harmonised_data_before_filter.csv"), row.names = FALSE)
      log_info("Harmonized data saved for inspection: %s", file.path(output_dir, 'harmonised_data_before_filter.csv'))
      return(list(
        outcome_id = outcome_ID,
        id.exposure = exposure_ID,
        status = "failed",
        error = "No SNPs after F statistic filtering",
        snps_after_filter = 0,
        f_exposure_min = NA,
        f_exposure_max = NA,
        f_exposure_mean = NA,
        passed_screening = FALSE
      ))
    }
    
    log_info("Filtered data: Remaining %d SNPs after F statistic > %d", nrow(dat_filtered), config$F_MIN)
    
    invalid_rows <- is.na(dat_filtered$beta.outcome) | is.na(dat_filtered$se.outcome) |
      !is.finite(dat_filtered$beta.outcome) | !is.finite(dat_filtered$se.outcome)
    if (any(invalid_rows)) {
      log_info("Removed %d invalid records (containing NA or non-finite values)", sum(invalid_rows))
      dat_filtered <- dat_filtered[!invalid_rows, ]
    }
    
    write.csv(dat_filtered, file.path(output_dir, "filtered_data.csv"), row.names = FALSE)
    log_info("Filtered data saved: %s", file.path(output_dir, 'filtered_data.csv'))
    
    write.csv(dat_filtered, file.path(output_dir, "filtered_data_list.csv"), row.names = FALSE)
    log_info("Filtered data list saved: %s", file.path(output_dir, 'filtered_data_list.csv'))
    
    # Perform MR analysis (uniform threshold <3 as insufficient)
    if (nrow(dat_filtered) < 3) {
      log_error("Insufficient number of SNPs (%d < 3) for exposure %s and outcome %s, unable to perform MR analysis", nrow(dat_filtered), exposure_ID, outcome_ID)
      empty_res <- data.table(
        exposure = character(),
        outcome = character(),
        id.exposure = character(),
        id.outcome = character(),
        method = character(),
        nsnp = numeric(),
        b = numeric(),
        se = numeric(),
        pval = numeric()
      )
      fwrite(empty_res, file.path(output_dir, "mr_results.csv"))
      fwrite(empty_res, file.path(output_dir, "mr_or_results.csv"))
      log_info("Insufficient SNPs, saved empty MR results: %s", file.path(output_dir, 'mr_results.csv'))
      return(list(
        outcome_id = outcome_ID,
        id.exposure = exposure_ID,
        status = "failed",
        error = "Insufficient SNPs for MR analysis",
        snps_after_filter = nrow(dat_filtered),
        f_exposure_min = NA,
        f_exposure_max = NA,
        f_exposure_mean = NA,
        passed_screening = FALSE
      ))
    }
    res <- mr(dat_filtered, method_list = c("mr_egger_regression", "mr_weighted_median", "mr_ivw", "mr_simple_mode", "mr_weighted_mode"))
    OR <- generate_odds_ratios(res)
    fwrite(res, file.path(output_dir, "mr_results.csv"))
    fwrite(OR, file.path(output_dir, "mr_or_results.csv"))
    log_info("MR analysis completed, results saved")
    
    # Perform MR-PRESSO analysis
    tryCatch({
      if (nrow(dat_filtered) < 3 || !all(c("beta.exposure", "beta.outcome", "se.exposure", "se.outcome") %in% colnames(dat_filtered))) {
        log_error("dat_filtered insufficient data or missing required columns, unable to perform MR-PRESSO analysis")
        presso_results <- data.table(
          Exposure = "beta.exposure",
          `MR Analysis` = "Raw",
          `Causal Estimate` = NA,
          Sd = NA,
          `T-stat` = NA,
          `P-value` = NA
        )
        fwrite(presso_results, file.path(output_dir, "mr_presso_main.csv"))
        outlier_results <- data.table(
          SNP = character(),
          RSSobs = numeric(),
          Pvalue = numeric()
        )
        fwrite(outlier_results, file.path(output_dir, "mr_presso_outliers.csv"))
        log_info("MR-PRESSO failed, saved empty results: %s", file.path(output_dir, 'mr_presso_main.csv'))
        log_info("MR-PRESSO failed, saved empty Outlier Test results: %s", file.path(output_dir, 'mr_presso_outliers.csv'))
      } else {
        log_info("Class of dat_filtered: %s", class(dat_filtered))
        dat_filtered <- as.data.frame(dat_filtered)
        log_info("Class after converting dat_filtered to data.frame: %s", class(dat_filtered))
        write.csv(dat_filtered, file.path(output_dir, "mr_presso_input.csv"), row.names = FALSE)
        log_info("MR-PRESSO input data saved: %s", file.path(output_dir, 'mr_presso_input.csv'))
        log_info("Column names of dat_filtered: %s", paste(colnames(dat_filtered), collapse=', '))
        log_info("Number of rows in dat_filtered: %d", nrow(dat_filtered))
        
        presso <- mr_presso(
          BetaOutcome = "beta.outcome",
          BetaExposure = "beta.exposure",
          SdOutcome = "se.outcome",
          SdExposure = "se.exposure",
          data = dat_filtered,
          OUTLIERtest = TRUE,
          DISTORTIONtest = TRUE,
          NbDistribution = config$mr_presso_params$NbDistribution,
          SignifThreshold = config$mr_presso_params$SignifThreshold
        )
        log_info("MR-PRESSO result structure: %s", str(presso))
        
        if (is.null(presso$`Main MR results`) || nrow(presso$`Main MR results`) == 0) {
          log_info("MR-PRESSO Main MR results empty, saving default results")
          presso_results <- data.table(
            Exposure = "beta.exposure",
            `MR Analysis` = "Raw",
            `Causal Estimate` = NA,
            Sd = NA,
            `T-stat` = NA,
            `P-value` = NA
          )
          fwrite(presso_results, file.path(output_dir, "mr_presso_main.csv"))
        } else {
          fwrite(presso$`Main MR results`, file.path(output_dir, "mr_presso_main.csv"))
          log_info("MR-PRESSO Main MR results saved: %s", file.path(output_dir, 'mr_presso_main.csv'))
        }
        
        if (is.null(presso$`MR-PRESSO results`$`Outlier Test`) || nrow(presso$`MR-PRESSO results`$`Outlier Test`) == 0) {
          log_info("MR-PRESSO no outliers detected, saving empty Outlier Test results")
          outlier_results <- data.table(
            SNP = character(),
            RSSobs = numeric(),
            Pvalue = numeric()
          )
          fwrite(outlier_results, file.path(output_dir, "mr_presso_outliers.csv"))
        } else {
          fwrite(presso$`MR-PRESSO results`$`Outlier Test`, file.path(output_dir, "mr_presso_outliers.csv"))
          log_info("MR-PRESSO Outlier Test results saved: %s", file.path(output_dir, 'mr_presso_outliers.csv'))
        }
        
        log_info("MR-PRESSO analysis completed")
      }
    }, error = function(e) {
      log_error("MR-PRESSO analysis failed: %s", e$message)
      log_info("Class of dat_filtered: %s", class(dat_filtered))
      log_info("Skipping MR-PRESSO analysis, continuing with subsequent steps")
      presso_results <- data.table(
        Exposure = "beta.exposure",
        `MR Analysis` = "Raw",
        `Causal Estimate` = NA,
        Sd = NA,
        `T-stat` = NA,
        `P-value` = NA
      )
      fwrite(presso_results, file.path(output_dir, "mr_presso_main.csv"))
      log_info("MR-PRESSO failed, saved empty results: %s", file.path(output_dir, 'mr_presso_main.csv'))
      
      outlier_results <- data.table(
        SNP = character(),
        RSSobs = numeric(),
        Pvalue = numeric()
      )
      fwrite(outlier_results, file.path(output_dir, "mr_presso_outliers.csv"))
      log_info("MR-PRESSO failed, saved empty Outlier Test results: %s", file.path(output_dir, 'mr_presso_outliers.csv'))
    })
    
    # Generate visualization charts (uniform threshold <3 as insufficient)
    tryCatch({
      if (nrow(dat_filtered) < 3 || !all(c("id.exposure", "id.outcome") %in% colnames(dat_filtered))) {
        log_error("dat_filtered insufficient data or missing required columns, unable to generate visualizations")
        empty_plot <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "Visualization failed: Insufficient data or missing columns")
        ggsave(file.path(output_dir, "scatter_plot.pdf"), plot = empty_plot, width = 6, height = 6)
        ggsave(file.path(output_dir, "forest_plot.pdf"), plot = empty_plot, width = 8, height = 6)
        ggsave(file.path(output_dir, "funnel_plot.pdf"), plot = empty_plot, width = 6, height = 6)
        log_info("Saved empty visualization charts: %s", file.path(output_dir, 'scatter_plot.pdf'))
      } else {
        p <- mr_scatter_plot(res, dat_filtered)
        ggsave(file.path(output_dir, "scatter_plot.pdf"), plot = p[[1]], width = 6, height = 6)
        log_info("Scatter plot generated: %s", file.path(output_dir, 'scatter_plot.pdf'))
        
        log_info("dat_filtered contains %d SNPs", nrow(dat_filtered))
        dat_subset <- dat_filtered[order(dat_filtered$F.exposure, decreasing = TRUE)[1:min(20, nrow(dat_filtered))], ]
        log_info("dat_subset contains %d SNPs", nrow(dat_subset))
        
        if (nrow(dat_subset) < 3) {
          log_info("dat_subset SNP count less than 3 (%d), unable to generate res_single.csv", nrow(dat_subset))
          empty_res_single <- data.table(
            exposure = character(),
            outcome = character(),
            id.exposure = character(),
            id.outcome = character(),
            samplesize = numeric(),
            SNP = character(),
            b = numeric(),
            se = numeric(),
            p = numeric()
          )
          fwrite(empty_res_single, file.path(output_dir, "res_single.csv"))
          log_info("Insufficient SNPs, saved empty res_single.csv: %s", file.path(output_dir, 'res_single.csv'))
          
          log_info("Individual SNP count less than 3 (%d), not generating forest plot", nrow(dat_subset))
          empty_plot <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = paste("Insufficient SNPs for forest plot (<3):", nrow(dat_subset), "SNPs"))
          ggsave(file.path(output_dir, "forest_plot.pdf"), plot = empty_plot, width = 8, height = 6)
          log_info("Saved placeholder forest plot: %s", file.path(output_dir, 'forest_plot.pdf'))
          
          log_info("Individual SNP count less than 3 (%d), not generating funnel plot", nrow(dat_subset))
          empty_plot <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = paste("Insufficient SNPs for funnel plot (<3):", nrow(dat_subset), "SNPs"))
          ggsave(file.path(output_dir, "funnel_plot.pdf"), plot = empty_plot, width = 6, height = 6)
          log_info("Saved placeholder funnel plot: %s", file.path(output_dir, 'funnel_plot.pdf'))
          
          res_single <- empty_res_single
        } else {
          tryCatch({
            res_single <- mr_singlesnp(dat_subset)
            log_info("mr_singlesnp completed, res_single contains %d records", nrow(res_single))
            fwrite(res_single, file.path(output_dir, "res_single.csv"))
            log_info("Single SNP MR results saved, containing %d SNPs: %s", nrow(res_single), file.path(output_dir, 'res_single.csv'))
          }, error = function(e) {
            log_error("mr_singlesnp failed: %s", e$message)
            empty_res_single <- data.table(
              exposure = character(),
              outcome = character(),
              id.exposure = character(),
              id.outcome = character(),
              samplesize = numeric(),
              SNP = character(),
              b = numeric(),
              se = numeric(),
              p = numeric()
            )
            fwrite(empty_res_single, file.path(output_dir, "res_single.csv"))
            log_info("mr_singlesnp failed, saved empty res_single.csv: %s", file.path(output_dir, 'res_single.csv'))
            res_single <- empty_res_single
          })
          
          res_single <- fread(file.path(output_dir, "res_single.csv"))
          
          res_single[, SNP := trimws(SNP)]
          res_single[, b := as.numeric(b)]
          res_single[, se := as.numeric(se)]
          
          if ("samplesize" %in% colnames(res_single)) {
            res_single[is.na(samplesize), samplesize := 0]
            log_info("Filled %d NA values in samplesize column with 0", sum(is.na(res_single$samplesize)))
          }
          
          log_info("Loaded data from res_single.csv, containing %d records", nrow(res_single))
          log_info("res_single column names: %s", paste(colnames(res_single), collapse=', '))
          log_info("res_single SNP column content: %s", paste(res_single$SNP, collapse=', '))
          log_info("res_single b column data type: %s", class(res_single$b))
          log_info("res_single se column data type: %s", class(res_single$se))
          
          if (!all(c("SNP", "b", "se") %in% colnames(res_single))) {
            log_error("res_single missing required columns: %s", paste(setdiff(c('SNP', 'b', 'se'), colnames(res_single)), collapse=', '))
            stop("res_single missing required columns")
          }
          
          na_rows <- res_single[is.na(b) | is.na(se), ]
          if (nrow(na_rows) > 0) {
            log_info("Found %d records removed due to NA in b or se columns", nrow(na_rows))
            log_info("SNPs of NA records: %s", paste(na_rows$SNP, collapse=', '))
          }
          
          non_finite_rows <- res_single[!is.finite(b) | !is.finite(se), ]
          if (nrow(non_finite_rows) > 0) {
            log_info("Found %d records removed due to non-finite values in b or se columns", nrow(non_finite_rows))
            log_info("SNPs of non-finite value records: %s", paste(non_finite_rows$SNP, collapse=', '))
          }
          
          res_single <- res_single[!is.na(b) & !is.na(se), ]
          res_single <- res_single[is.finite(b) & is.finite(se), ]
          log_info("res_single after filtering contains %d records", nrow(res_single))
          
          individual_snps <- res_single[!grepl("^All", res_single$SNP), ]
          num_individual_snps <- nrow(individual_snps)
          log_info("Individual SNP count (excluding summary estimates): %d", num_individual_snps)
          log_info("Individual SNPs: %s", paste(individual_snps$SNP, collapse=', '))
          
          if (num_individual_snps < 3) {
            log_info("Individual SNP count less than 3 (%d), not generating forest plot", num_individual_snps)
            empty_plot <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = paste("Insufficient SNPs for forest plot (<3):", num_individual_snps, "SNPs"))
            ggsave(file.path(output_dir, "forest_plot.pdf"), plot = empty_plot, width = 8, height = 6)
            log_info("Saved placeholder forest plot: %s", file.path(output_dir, 'forest_plot.pdf'))
          } else {
            summary_rows <- res_single[grepl("^All", res_single$SNP), ]
            if (nrow(summary_rows) == 0) {
              log_info("res_single missing summary estimates, adding IVW estimate")
              ivw_res <- mr(dat_filtered, method_list = "mr_ivw")
              ivw_row <- data.table(
                exposure = exposure_ID,
                outcome = outcome_ID,
                id.exposure = res_single$id.exposure[1],
                id.outcome = res_single$id.outcome[1],
                samplesize = NA,
                SNP = "All - Inverse variance weighted",
                b = ivw_res$b,
                se = ivw_res$se,
                p = ivw_res$pval
              )
              res_single <- rbind(res_single, ivw_row, fill = TRUE)
            }
            
            log_info("Attempting to generate forest plot...")
            p_forest <- mr_forest_plot(res_single)
            
            if (length(p_forest) == 0 || is.null(p_forest[[1]])) {
              log_error("mr_forest_plot returned empty result")
              empty_plot <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "Forest plot generation failed")
              ggsave(file.path(output_dir, "forest_plot.pdf"), plot = empty_plot, width = 8, height = 6)
              log_info("Saved empty forest plot: %s", file.path(output_dir, 'forest_plot.pdf'))
            } else {
              plot_height <- max(6, nrow(res_single) * 0.25)
              tryCatch({
                ggsave(file.path(output_dir, "forest_plot.pdf"), plot = p_forest[[1]], width = 8, height = plot_height, device = "pdf")
                log_info("Forest plot generated: %s", file.path(output_dir, 'forest_plot.pdf'))
              }, error = function(e) {
                log_error("ggsave failed to save forest plot: %s", e$message)
                empty_plot <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "Forest plot saving failed")
                ggsave(file.path(output_dir, "forest_plot.pdf"), plot = empty_plot, width = 8, height = 6)
                log_info("Saved empty forest plot: %s", file.path(output_dir, 'forest_plot.pdf'))
              })
            }
          }
          
          tryCatch({
            log_info("Starting to generate funnel plot, res_single contains %d records", nrow(res_single))
            log_info("res_single column names: %s", paste(colnames(res_single), collapse=', '))
            
            funnel_data <- res_single[, .(SNP, b, se, p)]
            fwrite(funnel_data, file.path(output_dir, "funnel_plot_data.csv"))
            log_info("Funnel plot data saved: %s", file.path(output_dir, 'funnel_plot_data.csv'))
            
            if (num_individual_snps < 3) {
              log_info("Individual SNP count less than 3 (%d), not generating funnel plot", num_individual_snps)
              empty_plot <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = paste("Insufficient SNPs for funnel plot (<3):", num_individual_snps, "SNPs"))
              ggsave(file.path(output_dir, "funnel_plot.pdf"), plot = empty_plot, width = 6, height = 6)
              log_info("Saved placeholder funnel plot: %s", file.path(output_dir, 'funnel_plot.pdf'))
            } else {
              p_funnel <- mr_funnel_plot(res_single)
              if (length(p_funnel) == 0 || is.null(p_funnel[[1]])) {
                log_error("mr_funnel_plot returned empty result")
                empty_plot <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "Funnel plot generation failed")
                ggsave(file.path(output_dir, "funnel_plot.pdf"), plot = empty_plot, width = 6, height = 6)
                log_info("Saved empty funnel plot: %s", file.path(output_dir, 'funnel_plot.pdf'))
              } else {
                tryCatch({
                  ggsave(file.path(output_dir, "funnel_plot.pdf"), plot = p_funnel[[1]], width = 6, height = 6, dpi = 300, device = "pdf")
                  log_info("Funnel plot generated: %s", file.path(output_dir, 'funnel_plot.pdf'))
                }, error = function(e) {
                  log_error("ggsave failed to save funnel plot: %s", e$message)
                  empty_plot <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "Funnel plot saving failed")
                  ggsave(file.path(output_dir, "funnel_plot.pdf"), plot = empty_plot, width = 6, height = 6)
                  log_info("Saved empty funnel plot: %s", file.path(output_dir, 'funnel_plot.pdf'))
                })
              }
            }
          }, error = function(e) {
            log_error("Funnel plot generation failed: %s", e$message)
            empty_plot <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = paste("Funnel plot generation failed:", e$message))
            ggsave(file.path(output_dir, "funnel_plot.pdf"), plot = empty_plot, width = 6, height = 6)
            log_info("Saved empty funnel plot: %s", file.path(output_dir, 'funnel_plot.pdf'))
          })
        }
      }
      
      # Perform leave-one-out analysis
      tryCatch({
        log_info("Starting leave-one-out analysis, dat_filtered contains %d SNPs", nrow(dat_filtered))
        log_info("dat_filtered column names: %s", paste(colnames(dat_filtered), collapse=', '))
        
        required_cols <- c("beta.exposure", "beta.outcome", "se.exposure", "se.outcome")
        missing_cols <- setdiff(required_cols, colnames(dat_filtered))
        if (length(missing_cols) > 0) {
          log_error("dat_filtered missing required columns: %s", paste(missing_cols, collapse=', '))
          empty_res_loo <- data.table(
            SNP = character(),
            b = numeric(),
            se = numeric(),
            pval = numeric()
          )
          fwrite(empty_res_loo, file.path(output_dir, "leave_one_out_results.csv"))
          log_info("Missing required columns, saved empty leave-one-out results: %s", file.path(output_dir, 'leave_one_out_results.csv'))
          empty_plot <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "Leave-one-out plot failed: Missing required columns")
          ggsave(file.path(output_dir, "leave_one_out_plot.pdf"), plot = empty_plot, width = 8, height = 6)
          log_info("Saved empty leave-one-out plot: %s", file.path(output_dir, 'leave_one_out_plot.pdf'))
        } else if (nrow(dat_filtered) < 3) {
          log_error("dat_filtered SNP count insufficient (<3), unable to perform leave-one-out analysis")
          empty_res_loo <- data.table(
            SNP = character(),
            b = numeric(),
            se = numeric(),
            pval = numeric()
          )
          fwrite(empty_res_loo, file.path(output_dir, "leave_one_out_results.csv"))
          log_info("Insufficient SNPs, saved empty leave-one-out results: %s", file.path(output_dir, 'leave_one_out_results.csv'))
          empty_plot <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "Leave-one-out plot failed: Insufficient SNPs")
          ggsave(file.path(output_dir, "leave_one_out_plot.pdf"), plot = empty_plot, width = 8, height = 6)
          log_info("Saved empty leave-one-out plot: %s", file.path(output_dir, 'leave_one_out_plot.pdf'))
        } else {
          res_loo <- mr_leaveoneout(dat_filtered)
          log_info("Leave-one-out analysis completed, res_loo contains %d records", nrow(res_loo))
          fwrite(res_loo, file.path(output_dir, "leave_one_out_results.csv"))
          log_info("Leave-one-out results saved: %s", file.path(output_dir, 'leave_one_out_results.csv'))
          
          p_loo <- mr_leaveoneout_plot(res_loo)
          if (length(p_loo) == 0 || is.null(p_loo[[1]])) {
            log_error("mr_leaveoneout_plot returned empty result")
            empty_plot <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "Leave-one-out plot generation failed")
            ggsave(file.path(output_dir, "leave_one_out_plot.pdf"), plot = empty_plot, width = 8, height = 6)
            log_info("Saved empty leave-one-out plot: %s", file.path(output_dir, 'leave_one_out_plot.pdf'))
          } else {
            tryCatch({
              ggsave(file.path(output_dir, "leave_one_out_plot.pdf"), plot = p_loo[[1]], width = 8, height = max(8, nrow(res_loo) * 0.3), dpi = 300, device = "pdf")
              log_info("Leave-one-out plot generated: %s", file.path(output_dir, 'leave_one_out_plot.pdf'))
            }, error = function(e) {
              log_error("ggsave failed to save leave-one-out plot: %s", e$message)
              empty_plot <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "Leave-one-out plot saving failed")
              ggsave(file.path(output_dir, "leave_one_out_plot.pdf"), plot = empty_plot, width = 8, height = 6)
              log_info("Saved empty leave-one-out plot: %s", file.path(output_dir, 'leave_one_out_plot.pdf'))
            })
          }
        }
      }, error = function(e) {
        log_error("Leave-one-out plot generation failed: %s", e$message)
        empty_res_loo <- data.table(
          SNP = character(),
          b = numeric(),
          se = numeric(),
          pval = numeric()
        )
        fwrite(empty_res_loo, file.path(output_dir, "leave_one_out_results.csv"))
        log_info("Leave-one-out analysis failed, saved empty leave-one-out results: %s", file.path(output_dir, 'leave_one_out_results.csv'))
        empty_plot <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = paste("Leave-one-out plot generation failed:", e$message))
        ggsave(file.path(output_dir, "leave_one_out_plot.pdf"), plot = empty_plot, width = 8, height = 6)
        log_info("Saved empty leave-one-out plot: %s", file.path(output_dir, 'leave_one_out_plot.pdf'))
      })
      
      log_info("Visualization charts and leave-one-out results generated")
    }, error = function(e) {
      log_error("Visualization failed: %s", e$message)
      log_info("Continuing with subsequent analysis")
    })
    
    # Perform heterogeneity test
    tryCatch({
      het_result <- mr_heterogeneity(dat_filtered, method_list = c("mr_egger_regression", "mr_ivw"))
      fwrite(het_result, file.path(output_dir, "heterogeneity_results.csv"))
      log_info("Heterogeneity test completed")
    }, error = function(e) {
      log_error("Heterogeneity test failed: %s", e$message)
      log_info("Continuing with subsequent analysis")
    })
    
    # Perform pleiotropy test
    tryCatch({
      egger_intercept <- mr_pleiotropy_test(dat_filtered)
      fwrite(egger_intercept, file.path(output_dir, "egger_intercept.csv"))
      log_info("Pleiotropy test completed")
    }, error = function(e) {
      log_error("Pleiotropy test failed: %s", e$message)
      log_info("Continuing with subsequent analysis")
    })
    
    # Screen and copy significant results
    focus_dir <- file.path(exposure_dir, "focus")
    if (!dir.exists(focus_dir)) dir.create(focus_dir, recursive = TRUE)
    
    mr_results <- fread(file.path(output_dir, "mr_results.csv"))
    if (nrow(mr_results) == 0 || !("method" %in% colnames(mr_results))) {
      log_info("mr_results empty or missing method column, skipping significance screening")
      return(list(
        outcome_id = outcome_ID,
        id.exposure = exposure_ID,
        status = "failed",
        error = "Invalid or empty MR results",
        snps_after_filter = nrow(dat_filtered),
        f_exposure_min = NA,
        f_exposure_max = NA,
        f_exposure_mean = NA,
        passed_screening = FALSE
      ))
    }
    
    ivw_result <- mr_results[method == "Inverse variance weighted"]
    passed_screening <- FALSE
    
    if (nrow(ivw_result) > 0) {
      ivw_pval <- as.numeric(ivw_result$pval)
      if (!is.na(ivw_pval) && ivw_pval < 0.05) {
        new_focus_dir <- file.path(focus_dir, paste0("results_", exposure_ID, "_", outcome_ID))
        if (!dir.exists(new_focus_dir)) dir.create(new_focus_dir, recursive = TRUE)
        file.copy(list.files(output_dir, full.names = TRUE), new_focus_dir, recursive = TRUE)
        log_info("Significant results copied to focus directory (IVW p-value = %.4f): %s", ivw_pval, new_focus_dir)
        
        result_all <- list(
          mr_results = mr_results,
          egger_intercept = if (file.exists(file.path(output_dir, "egger_intercept.csv")))
            fread(file.path(output_dir, "egger_intercept.csv")) else NULL,
          heterogeneity_results = if (file.exists(file.path(output_dir, "heterogeneity_results.csv")))
            fread(file.path(output_dir, "heterogeneity_results.csv")) else NULL,
          mr_or_results = fread(file.path(output_dir, "mr_or_results.csv")),
          mr_presso_main = if (file.exists(file.path(output_dir, "mr_presso_main.csv")))
            fread(file.path(output_dir, "mr_presso_main.csv")) else NULL,
          mr_presso_outliers = if (file.exists(file.path(output_dir, "mr_presso_outliers.csv")))
            fread(file.path(output_dir, "mr_presso_outliers.csv")) else NULL,
          leave_one_out_results = fread(file.path(output_dir, "leave_one_out_results.csv"))
        )
        
        result_file <- file.path(new_focus_dir, "Result_All.txt")
        file_conn <- file(result_file, open = "wt")
        for (name in names(result_all)) {
          if (!is.null(result_all[[name]])) {
            cat(paste(name, ":"), file = file_conn)
            cat("\n", file = file_conn)
            write.table(result_all[[name]], file = file_conn, append = TRUE, row.names = FALSE, col.names = TRUE, sep = "\t")
            cat("\n", file = file_conn)
          }
        }
        close(file_conn)
        log_info("Combined results saved to %s", result_file)
        passed_screening <- TRUE
      } else {
        log_info("Screening criteria not met: IVW p-value = %.4f >= 0.05", ivw_pval)
      }
    } else {
      log_info("No result for Inverse variance weighted method found in mr_results")
    }
    
    log_info("Analysis of exposure %s and outcome %s successfully completed", exposure_ID, outcome_ID)
    cat("Analysis successfully completed.\n")
    
    return(list(
      outcome_id = outcome_ID,
      id.exposure = exposure_ID,
      status = "success",
      snps_after_filter = nrow(dat_filtered),
      f_exposure_min = if ("F.exposure" %in% colnames(dat)) min(dat$F.exposure, na.rm=TRUE) else NA,
      f_exposure_max = if ("F.exposure" %in% colnames(dat)) max(dat$F.exposure, na.rm=TRUE) else NA,
      f_exposure_mean = if ("F.exposure" %in% colnames(dat)) mean(dat$F.exposure, na.rm=TRUE) else NA,
      passed_screening = passed_screening
    ))
  }, error = function(e) {
    log_error("Error processing exposure %s and outcome %s: %s", exposure_ID, outcome_ID, e$message)
    return(list(
      outcome_id = outcome_ID,
      id.exposure = exposure_ID,
      status = "failed",
      error = e$message,
      snps_after_filter = NA,
      f_exposure_min = NA,
      f_exposure_max = NA,
      f_exposure_mean = NA,
      passed_screening = FALSE
    ))
  })
}
# Process each outcome ID
num_batches <- ceiling(length(exposure_ids) / config$batch_size) # Calculate total batches
# Pre-create all exposure_dir to avoid parallel competition
for (exposure_ID in exposure_ids) {
  exposure_dir <- file.path(config$output_base_dir, paste0("exposure_", exposure_ID))
  if (!dir.exists(exposure_dir)) dir.create(exposure_dir, recursive = TRUE)
  log_info("Pre-created exposure ID directory: %s", exposure_dir)
}
# Global summary table paths (using sapply + USE.NAMES, more robust)
global_summary_paths <- setNames(
  sapply(exposure_ids, function(exp_id) {
    file.path(config$output_base_dir, paste0("exposure_", exp_id), "analysis_summary.csv")
  }, USE.NAMES = TRUE),
  exposure_ids
)
for (outcome_ID in outcome_ids) {
  log_info("Starting to process outcome ID: %s", outcome_ID)
  # Read outcome data (only once)
  outcome_file <- file.path(config$base_dir, "GWAS", paste0(outcome_ID, ".vcf.gz"))
  if (!file.exists(outcome_file)) {
    log_error("Outcome data file does not exist: %s", outcome_file)
    next
  }
  outcome_raw <- read_and_convert_vcf(outcome_file, outcome_ID)
  log_info("Outcome %s data read, containing %d SNPs", outcome_ID, nrow(outcome_raw))
  # Process exposure IDs in batches
  outcome_results <- list() # Temporary results for each outcome
  for (batch in 1:num_batches) {
    batch_start <- (batch - 1) * config$batch_size + 1
    batch_end <- min(batch * config$batch_size, length(exposure_ids))
    batch_ids <- exposure_ids[batch_start:batch_end]
    log_info("Processing batch %d/%d, containing %d exposure IDs", batch, num_batches, length(batch_ids))
    
    # Parallel process current batch of exposure IDs
    results_exp <- foreach(exposure_ID = batch_ids,
                           .packages = c("TwoSampleMR", "data.table", "MRPRESSO", "ggplot2", "fs", "logger", "vcfR", "dplyr"),
                           .export = c("process_exposure_outcome", "read_and_convert_vcf", "config"),
                           .errorhandling = "pass") %dopar% {
                             # Set independent log file
                             exposure_dir <- file.path(config$output_base_dir, paste0("exposure_", exposure_ID))
                             exposure_log_dir <- file.path(exposure_dir, "logs")
                             if (!dir.exists(exposure_log_dir)) dir.create(exposure_log_dir, recursive = TRUE)
                             exposure_log_file <- file.path(exposure_log_dir, paste0("mr_analysis_", exposure_ID, "_", outcome_ID, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
                             log_appender(appender_file(exposure_log_file))
                             log_threshold(INFO)
                             log_formatter(formatter_sprintf) # Set sprintf formatter in parallel environment
                             
                             log_info("Starting to process exposure ID: %s with outcome ID: %s", exposure_ID, outcome_ID)
                             
                             # Call processing function
                             result <- process_exposure_outcome(exposure_ID, outcome_ID, exposure_dir, outcome_raw)
                             return(result)
                           }
    
    # Collect current batch results, filter error objects
    valid_results <- results_exp[!sapply(results_exp, inherits, "error")]
    outcome_results <- append(outcome_results, valid_results)
    log_info("Batch %d/%d completed, collected %d valid results", batch, num_batches, length(valid_results))
    
    # Save current batch to global summary table
    if (length(valid_results) > 0) {
      summary_table <- rbindlist(lapply(valid_results, as.data.table), fill = TRUE)
      
      for (exp_id in unique(summary_table$id.exposure)) {
        # Skip invalid exp_id
        if (is.null(exp_id) || !nzchar(as.character(exp_id))) {
          log_warn("Skipping invalid exp_id: %s", exp_id)
          next
        }
        
        # Force convert to character and find path
        summary_path <- global_summary_paths[[as.character(exp_id)]]
        
        if (is.null(summary_path)) {
          log_error("exp_id not found in global_summary_paths: %s", exp_id)
          next
        }
        
        exp_summary <- summary_table[id.exposure == exp_id]
        
        if (file.exists(summary_path)) {
          fwrite(exp_summary, summary_path, append = TRUE, col.names = FALSE)
        } else {
          fwrite(exp_summary, summary_path, append = FALSE, col.names = TRUE)
        }
        log_info("Analysis summary table for exposure %s appended and saved: %s", exp_id, summary_path)
      }
    } else {
      log_info("Batch %d/%d no valid results", batch, num_batches)
    }
    
    # Release memory for current batch
    rm(results_exp, valid_results, summary_table)
    gc(full = TRUE)
    memory_usage <- format(object.size(environment()), units = "Mb")
    log_info("Batch %d/%d memory released, current memory usage: %s", batch, num_batches, memory_usage)
  }
  # Release outcome_raw and outcome_results
  rm(outcome_raw, outcome_results)
  gc(full = TRUE)
  log_info("Outcome %s processing completed, memory released", outcome_ID)
}
# Close parallel backend
stopImplicitCluster()
log_info("All outcome IDs processing completed")
cat("All analyses completed.\n")

# Close parallel backend
stopImplicitCluster()
log_info("All outcome IDs processing completed")
cat("All analyses completed.\n")

# Close parallel backend
stopImplicitCluster()
log_info("All outcome IDs processing completed")
cat("All analyses completed.\n")

# ============================================================
# BH_Correction.R
# ============================================================

library(rprojroot)
library(data.table)

# ── Path configuration ────────────────────────────────────────
project_root    <- find_rstudio_root_file()
single_gene_dir <- file.path(project_root, "outputs",
                              "06_UVMR_Analysis",
                              "02_UVMR_Single_Gene")

if (!dir.exists(single_gene_dir)) {
  stop("Directory not found: ", single_gene_dir,
       "\nPlease verify the path matches your actual output directory.")
}

# ── Parameters ────────────────────────────────────────────────
TARGET_METHOD <- "Inverse variance weighted"  # Primary MR method
BH_ALPHA      <- 0.05                         # Significance threshold
F_MIN         <- 10                           # Minimum F-statistic

# ── Step 1: Scan all exposure directories ─────────────────────
exposure_dirs <- list.dirs(single_gene_dir,
                            full.names = TRUE,
                            recursive  = FALSE)

# Keep only directories named exposure_XXX
exposure_dirs <- exposure_dirs[
  grepl("^exposure_", basename(exposure_dirs))
]

cat("Exposure directories found:", length(exposure_dirs), "\n")
if (length(exposure_dirs) == 0) {
  stop("No exposure_XXX directories found. Please check the path.")
}

# ── Step 2: Read mr_results.csv for every exposure × outcome ──
all_results <- list()

for (exp_dir in exposure_dirs) {
  exp_id <- sub("^exposure_", "", basename(exp_dir))

  # Find all results_XXX subdirectories under this exposure
  result_dirs <- list.dirs(exp_dir, full.names = TRUE, recursive = FALSE)
  result_dirs <- result_dirs[
    grepl("^results_", basename(result_dirs))
  ]

  if (length(result_dirs) == 0) {
    cat("  [SKIP]", exp_id, "— no results subdirectories found\n")
    next
  }

  exp_rows <- list()

  for (res_dir in result_dirs) {

    mr_file <- file.path(res_dir, "mr_results.csv")
    if (!file.exists(mr_file)) next

    mr_dat <- tryCatch(fread(mr_file), error = function(e) NULL)
    if (is.null(mr_dat) || nrow(mr_dat) == 0) next

    # Extract the IVW row
    ivw_row <- mr_dat[method == TARGET_METHOD]
    if (nrow(ivw_row) == 0) next

    # ── Read mean F-statistic from filtered_data.csv ──────────
    mean_f  <- NA_real_
    fd_file <- file.path(res_dir, "filtered_data.csv")
    if (file.exists(fd_file)) {
      fd <- tryCatch(fread(fd_file), error = function(e) NULL)
      if (!is.null(fd) && "F.exposure" %in% colnames(fd))
        mean_f <- mean(fd$F.exposure, na.rm = TRUE)
    }

    # ── Read heterogeneity Q-test p-value ─────────────────────
    het_Q_pval <- NA_real_
    het_file   <- file.path(res_dir, "heterogeneity_results.csv")
    if (file.exists(het_file)) {
      het <- tryCatch(fread(het_file), error = function(e) NULL)
      if (!is.null(het) && "Q_pval" %in% colnames(het)) {
        ivw_het <- het[method == TARGET_METHOD]
        if (nrow(ivw_het) > 0)
          het_Q_pval <- ivw_het$Q_pval[1]
      }
    }

    # ── Read MR-Egger intercept p-value ───────────────────────
    egger_intercept_pval <- NA_real_
    egger_file           <- file.path(res_dir, "egger_intercept.csv")
    if (file.exists(egger_file)) {
      eg <- tryCatch(fread(egger_file), error = function(e) NULL)
      if (!is.null(eg) && "pval" %in% colnames(eg))
        egger_intercept_pval <- eg$pval[1]
    }

    # ── Parse SNP count and outcome ID ───────────────────────
    n_snp <- if ("nsnp" %in% colnames(ivw_row)) ivw_row$nsnp[1] else NA

    # Directory format: results_<exposure_id>_<outcome_id>
    dir_name   <- basename(res_dir)
    outcome_id <- sub(paste0("^results_", exp_id, "_"), "", dir_name)

    exp_rows[[length(exp_rows) + 1]] <- data.table(
      exposure_id          = exp_id,
      outcome_id           = outcome_id,
      method               = TARGET_METHOD,
      ivw_b                = ivw_row$b[1],
      ivw_se               = ivw_row$se[1],
      ivw_pval             = ivw_row$pval[1],
      ivw_nsnp             = n_snp,
      mean_F               = mean_f,
      het_Q_pval           = het_Q_pval,
      egger_intercept_pval = egger_intercept_pval,
      result_dir           = res_dir
    )
  }

  if (length(exp_rows) == 0) {
    cat("  [SKIP]", exp_id, "— all result files are empty or missing\n")
    next
  }

  exp_dt <- rbindlist(exp_rows, fill = TRUE)

  # ── Step 3: Per-exposure BH correction across outcomes ───────
  exp_dt[, ivw_pval_BH    := p.adjust(ivw_pval, method = "BH")]
  exp_dt[, significant_BH := !is.na(ivw_pval_BH) & ivw_pval_BH < BH_ALPHA]

  # ── QC flags ─────────────────────────────────────────────────
  exp_dt[, F_pass     := !is.na(mean_F) & mean_F > F_MIN]
  exp_dt[, het_pass   := is.na(het_Q_pval) | het_Q_pval > 0.05]
  exp_dt[, egger_pass := is.na(egger_intercept_pval) |
                           egger_intercept_pval > 0.05]
  exp_dt[, all_QC_pass := F_pass & het_pass & egger_pass]

  # ── OR conversion ─────────────────────────────────────────────
  exp_dt[, OR          := exp(ivw_b)]
  exp_dt[, OR_lower95  := exp(ivw_b - 1.96 * ivw_se)]
  exp_dt[, OR_upper95  := exp(ivw_b + 1.96 * ivw_se)]

  # ── Save per-exposure BH file ─────────────────────────────────
  out_csv <- file.path(exp_dir, paste0(exp_id, "_BH_corrected.csv"))
  fwrite(exp_dt, out_csv)

  n_sig <- sum(exp_dt$significant_BH,                    na.rm = TRUE)
  n_qc  <- sum(exp_dt$significant_BH & exp_dt$all_QC_pass, na.rm = TRUE)
  cat(sprintf(
    "  [OK] %-15s | outcomes: %2d | BH sig: %d | BH+QC: %d | saved: %s\n",
    exp_id, nrow(exp_dt), n_sig, n_qc, basename(out_csv)
  ))

  all_results[[exp_id]] <- exp_dt
}

# ── Step 4: Merge all exposures and apply global BH correction ─
if (length(all_results) == 0) {
  stop("No valid MR results were loaded. ",
       "Please check that mr_results.csv files exist.")
}

all_dt <- rbindlist(all_results, fill = TRUE)

# Global BH correction across all exposure × outcome combinations
all_dt[, ivw_pval_BH_global     := p.adjust(ivw_pval, method = "BH")]
all_dt[, significant_BH_global  := !is.na(ivw_pval_BH_global) &
                                      ivw_pval_BH_global < BH_ALPHA]

# Sort by per-exposure BH p-value ascending
setorder(all_dt, ivw_pval_BH)

# Save global summary table
global_out <- file.path(single_gene_dir, "ALL_exposures_BH_corrected.csv")
fwrite(all_dt, global_out)

# ── Step 5: Print summary ──────────────────────────────────────
cat("\n", strrep("=", 60), "\n")
cat("BH correction complete\n")
cat(strrep("=", 60), "\n")
cat(sprintf("Total analyses:                     %d\n",  nrow(all_dt)))
cat(sprintf("Per-exposure BH significant (p<0.05): %d\n",
            sum(all_dt$significant_BH, na.rm = TRUE)))
cat(sprintf("Global BH significant (p<0.05):       %d\n",
            sum(all_dt$significant_BH_global, na.rm = TRUE)))
cat(sprintf("BH significant + all QC passed:       %d\n",
            sum(all_dt$significant_BH & all_dt$all_QC_pass, na.rm = TRUE)))
cat("\nGlobal summary table saved to:\n ", global_out, "\n")

# ── Step 6: Print significant results preview ─────────────────
sig_dt <- all_dt[significant_BH == TRUE & all_QC_pass == TRUE]

if (nrow(sig_dt) > 0) {
  cat("\n── Results passing per-exposure BH correction + all QC ──\n")
  print(sig_dt[, .(
    exposure_id,
    outcome_id,
    ivw_b                = round(ivw_b, 6),
    ivw_se               = round(ivw_se, 6),
    ivw_pval             = formatC(ivw_pval,    format = "e", digits = 3),
    ivw_pval_BH          = formatC(ivw_pval_BH, format = "e", digits = 3),
    OR                   = round(OR, 4),
    OR_lower95           = round(OR_lower95, 4),
    OR_upper95           = round(OR_upper95, 4),
    mean_F               = round(mean_F, 1),
    het_Q_pval           = round(het_Q_pval, 3),
    egger_intercept_pval = round(egger_intercept_pval, 3)
  )])
} else {
  cat("\n── No results passed both BH correction and all QC filters ──\n")
  cat("Tip: results passing BH correction only (QC filters relaxed):\n")
  sig_any <- all_dt[significant_BH == TRUE]
  if (nrow(sig_any) > 0) {
    print(sig_any[, .(
      exposure_id,
      outcome_id,
      ivw_pval    = formatC(ivw_pval,    format = "e", digits = 3),
      ivw_pval_BH = formatC(ivw_pval_BH, format = "e", digits = 3),
      OR          = round(OR, 4),
      F_pass, het_pass, egger_pass
    )])
  } else {
    cat("No results passed BH correction under any QC condition.\n")
  }
}