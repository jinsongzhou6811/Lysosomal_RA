# Title: GEO Preprocessing Pipeline (English Version)
# Description: Process GSE55235, GSE55457, GSE55584. Perform quantile normalization,
# batch effect correction, and differential expression analysis.
# Keep all lysosomal genes. Generate single-gene level results.
# Input: ./data/ (series matrix + GPL files + lysosome_genes.csv)
# Output: ./outputs/01_Differential_Expression_Analysis/results/

# Load required R packages
library(GEOquery)      # Access GEO database
library(limma)         # Differential expression analysis
library(dplyr)         # Data manipulation
library(tidyr)         # Splitting gene symbol columns
library(sva)           # Batch effect correction
library(ggplot2)       # PCA and boxplots
library(pheatmap)      # Correlation heatmaps
library(data.table)    # Fast reading of platform files
library(reshape2)      # Data reshaping (e.g., melt for boxplots)
library(viridis)       # Color-blind friendly palettes
library(vegan)         # PERMANOVA for batch effect quantification [FIX-1]

writeLines("Libraries loaded: GEOquery, limma, dplyr, tidyr, sva, ggplot2, pheatmap, data.table, viridis", stderr())

# ========================================================
# install.packages("rprojroot")   
library(rprojroot)
project_root <- find_rstudio_root_file()
data_dir     <- file.path(project_root, "data")
outputs_dir  <- file.path(project_root, "outputs")


module_dir   <- file.path(outputs_dir, "01_Differential_Expression_Analysis")
preprocess_dir <- file.path(module_dir, "001_GEO Preprocessing")   

writeLines(paste("Project root detected:", project_root), stderr())
writeLines(paste("Data directory:", data_dir), stderr())
writeLines(paste("001 Preprocessing output:", preprocess_dir), stderr())
# ========================================================



# Set output directory → 001_GEO Preprocessing
module_name <- "01_Differential_Expression_Analysis"
base_dir <- preprocess_dir                    

dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "logs"),     showWarnings = FALSE)
dir.create(file.path(base_dir, "plots"),    showWarnings = FALSE)
dir.create(file.path(base_dir, "reports"),  showWarnings = FALSE)
dir.create(file.path(base_dir, "results"),  showWarnings = FALSE)

# Define log file path
log_timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
log_file <- file.path(base_dir, "logs", paste0("geo_analysis_", log_timestamp, ".log"))

# Define logging function
write_log <- function(message, file = log_file) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- sprintf("[%s] %s\n", timestamp, message)
  cat(log_entry)
  file_con <- file(file, open="a")
  writeLines(log_entry, file_con)
  close(file_con)
}

# Define function to read GEO series matrix metadata
read_series_matrix <- function(matrix_file) {
  write_log(paste("Reading series matrix metadata:", matrix_file))
  tryCatch({
    con <- gzfile(matrix_file)
    lines <- readLines(con)
    close(con)
    sample_lines <- lines[grep("^!Sample_", lines)]
    if (length(sample_lines) == 0) {
      write_log("No !Sample_ metadata found")
      return(NULL)
    }
    geo_lines <- sample_lines[grep("^!Sample_geo_accession", sample_lines)]
    if (length(geo_lines) == 0) {
      write_log("No !Sample_geo_accession found")
      return(NULL)
    }
    geo_values <- strsplit(geo_lines, "\t")[[1]][-1]
    geo_values <- gsub('^"|"$', '', geo_values)
    n_samples <- length(geo_values)
    standard_fields <- c("geo_accession", "platform_id", "characteristics_ch1")
    sample_data <- data.frame(matrix(NA, nrow=n_samples, ncol=length(standard_fields)))
    colnames(sample_data) <- standard_fields
    rownames(sample_data) <- geo_values
    for (field in standard_fields) {
      field_lines <- sample_lines[grep(paste0("^!Sample_", field, "\t"), sample_lines)]
      if (length(field_lines) > 0) {
        values <- strsplit(field_lines[1], "\t")[[1]][-1]
        if (length(values) == n_samples) {
          if (field == "platform_id") {
            values <- gsub('^"|"$', '', values)
          }
          sample_data[, field] <- values
        }
      }
    }
    sample_data$geo_accession <- gsub('^"|"$', '', sample_data$geo_accession)
    write_log(paste("Successfully read series matrix metadata:", matrix_file))
    write_log(paste("Sample geo_accessions:", paste(head(sample_data$geo_accession, 5), collapse=", ")))
    return(sample_data)
  }, error = function(e) {
    write_log(paste("Error reading series matrix metadata:", matrix_file, "-", e$message))
    return(NULL)
  })
}

# Define function to load lysosome genes
load_lysosome_genes <- function() {
  write_log("Loading lysosome_genes.csv")
  tryCatch({
    df <- read.csv(file.path(data_dir, "lysosome_genes.csv"))
    expected_cols <- c("GeneID", "Symbol")
    if (!all(expected_cols %in% colnames(df))) {
      write_log(paste("ERROR: lysosome_genes.csv missing required columns. Expected:", paste(expected_cols, collapse=", "), "Found:", paste(colnames(df), collapse=", ")))
      stop("Invalid lysosome_genes.csv structure")
    }
    df$GeneID <- as.character(df$GeneID)
    df$GeneID <- trimws(df$GeneID)
    df$GeneID <- sapply(strsplit(df$GeneID, ";"), function(x) x[1])
    invalid_ids <- df$GeneID[!grepl("^[0-9]+$", df$GeneID) | df$GeneID == ""]
    if (length(invalid_ids) > 0) {
      write_log(paste("WARNING: Found", length(invalid_ids), "invalid GeneIDs:", paste(head(invalid_ids, 5), collapse=", ")))
    }
    if (nrow(df) == 0) {
      write_log("ERROR: lysosome_genes.csv is empty")
      stop("Empty lysosome_genes.csv")
    }
    if (any(duplicated(df$GeneID))) {
      write_log("WARNING: Found duplicated GeneIDs:", paste(df$GeneID[duplicated(df$GeneID)], collapse=", "))
    }
    write_log(paste("Loaded lysosome_genes.csv with", nrow(df), "rows"))
    return(df)
  }, error = function(e) {
    write_log(paste("ERROR loading lysosome_genes.csv:", e$message))
    return(NULL)
  })
}

# Define function to process GEO dataset
process_geo <- function(geo_id, expr_combined, sample_info, unique_probes, probe_to_gene, lysosome_genes) {
  write_log(sprintf("Starting processing for GEO dataset: %s", geo_id))
  contrast_name <- NULL
  tryCatch({
    geo_dir <- file.path(base_dir, "results", geo_id)
    write_log(sprintf("Checking directory: %s", geo_dir))
    if (!dir.exists(geo_dir)) {
      dir.create(geo_dir, recursive = TRUE)
      write_log(sprintf("Created directory: %s", geo_dir))
    }
    
    local_file <- file.path(data_dir, paste0(geo_id, "_series_matrix.txt.gz"))
    write_log(sprintf("Checking series matrix file: %s", local_file))
    if (!file.exists(local_file)) {
      write_log(sprintf("ERROR: Series matrix file not found: %s", local_file))
      stop(sprintf("File not found: %s", local_file))
    }
    write_log(sprintf("Loading series matrix file: %s", local_file))
    gse <- getGEO(filename = local_file, GSEMatrix = TRUE, getGPL = FALSE)
    write_log("Series matrix file loaded successfully")
    
    if (inherits(gse, "ExpressionSet")) {
      exprs <- exprs(gse)
      pData <- pData(gse)
      write_log("Extracted exprs and pData from ExpressionSet")
    } else if (inherits(gse, "list") && length(gse) > 0 && inherits(gse[[1]], "ExpressionSet")) {
      exprs <- exprs(gse[[1]])
      pData <- pData(gse[[1]])
      write_log("Extracted exprs and pData from list element")
    } else {
      write_log(sprintf("ERROR: Invalid gse structure for file: %s", local_file))
      stop(sprintf("Invalid gse structure for file: %s", local_file))
    }
    write_log(sprintf("Number of samples in exprs: %d", ncol(exprs)))
    write_log(sprintf("Dimensions of exprs: %d rows, %d columns", nrow(exprs), ncol(exprs)))
    
    platform <- gse@annotation
    write_log(sprintf("Platform ID: %s", platform))
    gpl_local_file <- file.path(data_dir, paste0(platform, ".soft.gz"))
    write_log(sprintf("Checking platform annotation file: %s", gpl_local_file))
    if (!file.exists(gpl_local_file)) {
      write_log(sprintf("ERROR: Platform file not found: %s", gpl_local_file))
      stop(sprintf("Platform file not found: %s", gpl_local_file))
    }
    write_log(sprintf("Loading platform file: %s", gpl_local_file))
    file_size <- file.info(gpl_local_file)$size
    write_log(sprintf("Platform file size: %d bytes", file_size))
    if (file_size < 10000) {
      write_log(sprintf("ERROR: Platform file appears incomplete: %s", gpl_local_file))
      stop(sprintf("Platform file appears incomplete: %s", gpl_local_file))
    }
    gpl <- getGEO(filename = gpl_local_file)
    write_log("Platform annotation file loaded successfully")
    gpl_table <- Table(gpl)
    write_log(sprintf("Dimensions of gpl_table: %d rows, %d columns", nrow(gpl_table), ncol(gpl_table)))
    
    id_col <- intersect(c("ID", "ProbeID", "ID_REF"), colnames(gpl_table))[1]
    if (any(duplicated(gpl_table[[id_col]]))) {
      write_log("WARNING: Duplicate probe IDs found in GPL table, removing duplicates")
      gpl_table <- gpl_table[!duplicated(gpl_table[[id_col]]), ]
    }
    
    annotation_file <- file.path(geo_dir, paste0(platform, "_annotation.csv"))
    write.csv(gpl_table, annotation_file, row.names = FALSE)
    write_log(sprintf("Saved platform annotation to: %s", annotation_file))
    
    symbol_col <- intersect(c("Gene Symbol", "Symbol", "SYMBOL"), colnames(gpl_table))[1]
    entrez_col <- intersect(c("ENTREZ_GENE_ID"), colnames(gpl_table))[1]
    write_log(sprintf("Selected columns: ID=%s, Gene Symbol=%s, Entrez Gene ID=%s", id_col %||% "NA", symbol_col %||% "NA", entrez_col %||% "NA"))
    if (is.na(id_col)) {
      write_log("ERROR: Probe ID column not found in GPL table")
      stop("Probe ID column not found in GPL table")
    }
    if (is.na(symbol_col)) {
      write_log("WARNING: Gene Symbol column not found, proceeding with ID only")
      annotation <- data.frame(ID = gpl_table[[id_col]])
    } else {
      if (!is.na(entrez_col)) {
        annotation <- gpl_table[, c(id_col, symbol_col, entrez_col)]
        colnames(annotation) <- c("ID", "Gene Symbol", "Entrez_Gene_ID")
      } else {
        annotation <- gpl_table[, c(id_col, symbol_col)]
        colnames(annotation) <- c("ID", "Gene Symbol")
      }
    }
    write_log(sprintf("Annotation data frame created with %d rows, %d columns", nrow(annotation), ncol(annotation)))
    
    if ("Entrez_Gene_ID" %in% colnames(annotation)) {
      annotation <- annotation %>%
        mutate(
          GS_Split = strsplit(`Gene Symbol`, " /// "),
          ES_Split = strsplit(Entrez_Gene_ID, " /// "),
          Valid_Length = mapply(function(gs, es) length(gs) == length(es), GS_Split, ES_Split)
        )
      if (any(!annotation$Valid_Length)) {
        write_log(sprintf("WARNING: Mismatch between Gene Symbol and Entrez_Gene_ID lengths for %d probes", sum(!annotation$Valid_Length)))
      }
      annotation <- annotation %>% select(-GS_Split, -ES_Split, -Valid_Length)
    }
    
    selected_samples <- sample_info$Sample[sample_info$GEO == geo_id]
    write_log(sprintf("Selected samples for %s: %s", geo_id, paste(selected_samples, collapse=", ")))
    missing_samples <- selected_samples[!selected_samples %in% colnames(expr_combined)]
    if (length(missing_samples) > 0) {
      write_log(sprintf("ERROR: Samples not found in expr_combined for %s: %s", geo_id, paste(missing_samples, collapse=", ")))
      stop(sprintf("Samples not found in expr_combined for %s", geo_id))
    }
    
    write_log("Merging expression data with probe annotations")
    expr_subset <- expr_combined[, selected_samples, drop = FALSE]
    exprs_annotated <- merge(annotation, as.data.frame(expr_subset), by.x = "ID", by.y = "row.names")
    write_log(sprintf("Merged data frame dimensions: %d rows, %d columns", nrow(exprs_annotated), ncol(exprs_annotated)))
    
    if (!is.null(lysosome_genes)) {
      lysosome_in_expr <- annotation$Entrez_Gene_ID[annotation$ID %in% rownames(expr_subset)] %in% lysosome_genes$GeneID
      write_log(sprintf("Lysosome genes in %s: %d out of %d", geo_id, sum(lysosome_in_expr), nrow(lysosome_genes)))
      if (sum(lysosome_in_expr) < nrow(lysosome_genes)) {
        missing_lysosome <- lysosome_genes$GeneID[!lysosome_genes$GeneID %in% annotation$Entrez_Gene_ID[annotation$ID %in% rownames(expr_subset)]]
        write_log(sprintf("Missing lysosome GeneIDs in %s: %s", geo_id, paste(head(missing_lysosome, 5), collapse=", ")))
      }
    }
    
    group <- factor(sample_info$Group[sample_info$GEO == geo_id], levels = unique(sample_info$Group[sample_info$GEO == geo_id]))
    if (geo_id == "GSE55235") {
      contrast_name <- "RA-ND"
      write_log("GSE55235 grouping defined: 10 ND, 10 OA, 10 RA")
    } else if (geo_id == "GSE55457") {
      contrast_name <- "RA-Control"
      write_log("GSE55457 grouping defined: 10 Control, 10 RA, 13 OA")
    } else if (geo_id == "GSE55584") {
      contrast_name <- "RA-OA"
      write_log("GSE55584 grouping defined: 10 RA, 6 OA")
    } else {
      write_log(sprintf("ERROR: Unsupported GEO dataset: %s", geo_id))
      stop(sprintf("Unsupported GEO dataset: %s", geo_id))
    }
    write_log(sprintf("contrast_name set to: %s", contrast_name))
    
    write_log(sprintf("Number of samples in exprs: %d", ncol(expr_subset)))
    write_log(sprintf("Number of samples in group: %d", length(group)))
    if (ncol(expr_subset) != length(group)) {
      write_log(sprintf("ERROR: Sample count mismatch: exprs has %d columns, group has %d elements", ncol(expr_subset), length(group)))
      stop("Sample count mismatch between exprs and group")
    }
    write_log(sprintf("Group levels: %s", paste(levels(group), collapse=", ")))
    
    write_log("Generating design matrix")
    design <- model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    write_log(sprintf("Design matrix dimensions: %d rows, %d columns", nrow(design), ncol(design)))
    
    write_log("Fitting linear model with limma")
    fit <- lmFit(expr_subset, design)
    write_log("Linear model fitted successfully")
    contrast.matrix <- makeContrasts(contrasts=contrast_name, levels=design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    write_log("eBayes analysis completed")
    
    write_log("Extracting differential expression results")
    all_results <- topTable(fit2, coef=1, number=Inf, genelist=rownames(expr_subset))
    if (!"ID" %in% colnames(all_results)) {
      all_results$ID <- rownames(all_results)
      write_log("Added ID column to all_results from rownames")
    }
    
    output_dir <- file.path(geo_dir, "Results_All")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
      write_log(sprintf("Created output directory: %s", output_dir))
    }
    diff_expr_file <- file.path(output_dir, paste0(geo_id, "_diff_expr.csv"))
    diff_expr_normalized_file <- file.path(output_dir, paste0(geo_id, "_diff_expr_normalized.csv"))
    
    write_log("Merging differential expression results with annotations")
    all_results_merged <- merge(annotation, all_results, by="ID", all.y=TRUE)
    
    expanded_results <- all_results_merged %>%
      rowwise() %>%
      mutate(
        Gene_Symbol_Split = list(strsplit(as.character(`Gene Symbol`), " /// ")[[1]]),
        Entrez_Split = if ("Entrez_Gene_ID" %in% colnames(.)) {
          list(strsplit(as.character(Entrez_Gene_ID), " /// ")[[1]])
        } else {
          list(rep(NA, length(Gene_Symbol_Split[[1]])))
        }
      ) %>%
      mutate(
        Entrez_Split = list({
          gs_len <- length(Gene_Symbol_Split)
          es_len <- length(Entrez_Split)
          if (es_len < gs_len) {
            c(Entrez_Split, rep(NA, gs_len - es_len))
          } else if (es_len > gs_len) {
            Entrez_Split[1:gs_len]
          } else {
            Entrez_Split
          }
        })
      ) %>%
      unnest(cols = c(Gene_Symbol_Split, Entrez_Split), keep_empty = TRUE) %>%
      dplyr::rename(Gene_Symbol_Single = Gene_Symbol_Split, All_Gene_Symbols = `Gene Symbol`) %>%
      dplyr::rename(`Gene Symbol` = Gene_Symbol_Single, Entrez_Gene_ID_Single = Entrez_Split) %>%
      dplyr::rename(Entrez_Gene_ID = Entrez_Gene_ID_Single, All_Entrez_Gene_IDs = Entrez_Gene_ID) %>%
      mutate(Is_Multi_Gene = grepl("///", All_Gene_Symbols))
    
    required_columns <- c("ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "Gene Symbol", "All_Gene_Symbols", "Is_Multi_Gene")
    if ("Entrez_Gene_ID" %in% colnames(expanded_results)) {
      required_columns <- c(required_columns, "Entrez_Gene_ID", "All_Entrez_Gene_IDs")
    }
    expanded_results <- expanded_results[, required_columns, drop=FALSE]
    
    write.csv(expanded_results, diff_expr_file, row.names = FALSE)
    write_log(sprintf("Saved differential expression results to: %s", diff_expr_file))
    write.csv(expanded_results, diff_expr_normalized_file, row.names = FALSE)
    write_log(sprintf("Saved normalized differential expression results to: %s", diff_expr_normalized_file))
    
    single_gene_file <- file.path(output_dir, paste0(geo_id, "_diff_expr_normalized_single_gene.csv"))
    single_gene_results <- expanded_results %>%
      filter(!is.na(`Gene Symbol`) & `Gene Symbol` != "") %>%
      arrange(`Gene Symbol`, adj.P.Val) %>%
      group_by(`Gene Symbol`) %>%
      slice_head(n = 1) %>%
      ungroup()
    write.csv(single_gene_results, single_gene_file, row.names = FALSE)
    write_log(sprintf("Saved single-gene differential expression results to: %s", single_gene_file))
    
    write_log(sprintf("Analysis completed for %s", geo_id))
  }, error = function(e) {
    write_log(sprintf("ERROR: Failed to process GEO dataset %s: %s", geo_id, conditionMessage(e)))
    stop(e)
  })
}

# Main Pipeline
write_log("Starting GEO analysis pipeline")
lysosome_genes <- load_lysosome_genes()
geo_list <- c("GSE55235", "GSE55457", "GSE55584")
write_log(sprintf("GEO datasets to process: %s", paste(geo_list, collapse=", ")))
write_log("Loading and merging GEO datasets")
expr_list <- list()
sample_info_list <- list()
probe_to_gene <- NULL
for (geo_id in geo_list) {
  local_file <- file.path(data_dir, paste0(geo_id, "_series_matrix.txt.gz"))
  write_log(sprintf("Loading series matrix file: %s", local_file))
  gse <- getGEO(filename = local_file, GSEMatrix = TRUE, getGPL = FALSE)
  
  if (inherits(gse, "ExpressionSet")) {
    series_matrix <- exprs(gse)
  } else if (inherits(gse, "list") && length(gse) > 0 && inherits(gse[[1]], "ExpressionSet")) {
    series_matrix <- exprs(gse[[1]])
  } else {
    write_log(sprintf("ERROR: Invalid gse structure for file: %s", local_file))
    stop(sprintf("Invalid gse structure for file: %s", local_file))
  }
  
  pheno_data <- read_series_matrix(local_file)
  if (is.null(pheno_data)) {
    write_log(sprintf("ERROR: Failed to read series matrix metadata for %s", geo_id))
    next
  }
  
  samples <- pheno_data$geo_accession
  write_log(sprintf("Samples for %s: %s", geo_id, paste(head(samples, 5), collapse=", ")))
  if (geo_id == "GSE55235") {
    groups <- c(rep("ND", 10), rep("OA", 10), rep("RA", 10))
  } else if (geo_id == "GSE55457") {
    groups <- c(rep("Control", 10), rep("RA", 10), rep("OA", 13))
  } else if (geo_id == "GSE55584") {
    groups <- c(rep("RA", 10), rep("OA", 6))
  }
  if (length(samples) != length(groups)) {
    write_log(sprintf("ERROR: Sample count mismatch for %s: %d samples, %d groups", geo_id, length(samples), length(groups)))
    stop(sprintf("Sample count mismatch for %s", geo_id))
  }
  labels <- ifelse(groups == "RA", 1, 0)
  sample_info_list[[geo_id]] <- data.frame(Sample = samples, Group = groups, Label = labels, GEO = geo_id)
  
  platform <- gse@annotation
  gpl_file <- file.path(data_dir, paste0(platform, ".soft.gz"))
  gpl <- getGEO(filename = gpl_file)
  gpl_table <- Table(gpl)
  id_col <- intersect(c("ID", "ProbeID", "ID_REF"), colnames(gpl_table))[1]
  symbol_col <- intersect(c("Gene Symbol", "Symbol", "SYMBOL"), colnames(gpl_table))[1]
  entrez_col <- intersect(c("ENTREZ_GENE_ID"), colnames(gpl_table))[1]
  if (is.na(symbol_col)) {
    gpl_table$"Gene Symbol" <- gpl_table[[id_col]]
    write_log(sprintf("WARNING: No Gene Symbol column found for %s, using probe IDs", geo_id))
  } else {
    gpl_table$"Gene Symbol" <- gpl_table[[symbol_col]]
  }
  if (is.na(entrez_col)) {
    gpl_table$ENTREZ_GENE_ID <- gpl_table[[id_col]]
    write_log(sprintf("WARNING: No Entrez Gene ID column found for %s, using probe IDs", geo_id))
  } else {
    gpl_table$ENTREZ_GENE_ID <- gpl_table[[entrez_col]]
  }
  gpl_table$ENTREZ_GENE_ID <- sapply(strsplit(gpl_table$ENTREZ_GENE_ID, " /// "), function(x) x[1])
  gpl_table$ENTREZ_GENE_ID[is.na(gpl_table$ENTREZ_GENE_ID) | gpl_table$ENTREZ_GENE_ID == ""] <- gpl_table[[id_col]][is.na(gpl_table$ENTREZ_GENE_ID) | gpl_table$ENTREZ_GENE_ID == ""]
  probe_to_gene_tmp <- setNames(gpl_table$ENTREZ_GENE_ID, gpl_table[[id_col]])
  
  if (is.null(probe_to_gene)) {
    probe_to_gene <- probe_to_gene_tmp
  } else {
    common_probes <- intersect(names(probe_to_gene), names(probe_to_gene_tmp))
    probe_to_gene <- probe_to_gene[common_probes]
    probe_to_gene_tmp <- probe_to_gene_tmp[common_probes]
    probe_to_gene <- probe_to_gene[probe_to_gene == probe_to_gene_tmp]
  }
  
  expr_list[[geo_id]] <- series_matrix
}
unique_probes <- Reduce(intersect, lapply(expr_list, rownames))
write_log(sprintf("Number of common probes: %d", length(unique_probes)))
expr_combined <- do.call(cbind, lapply(expr_list, function(x) x[unique_probes, , drop = FALSE]))
write_log(sprintf("Combined expression matrix dimensions: %d rows, %d columns", nrow(expr_combined), ncol(expr_combined)))
write_log(sprintf("Column names of expr_combined: %s", paste(head(colnames(expr_combined), 5), collapse=", ")))
sample_info <- do.call(rbind, sample_info_list)
write.csv(sample_info, file.path(base_dir, "results", "sample_info.csv"), row.names = FALSE)
write_log("Saved sample information to sample_info.csv")
write_log(sprintf("Sample names in sample_info: %s", paste(head(sample_info$Sample, 5), collapse=", ")))
total_ra <- sum(sample_info$Group == "RA")
total_non_ra <- nrow(sample_info) - total_ra
write_log(sprintf("Total RA samples across datasets: %d, Total non-RA samples: %d", total_ra, total_non_ra))
if (!all(sample_info$Sample %in% colnames(expr_combined))) {
  missing_samples <- sample_info$Sample[!sample_info$Sample %in% colnames(expr_combined)]
  write_log(sprintf("ERROR: Sample names not found in expr_combined: %s", paste(missing_samples, collapse=", ")))
  stop("Sample names in sample_info do not match expr_combined")
}
write_log("Filtering missing values")
na_count_total <- apply(expr_combined, 1, function(x) sum(is.na(x)))
na_prop_total <- na_count_total / ncol(expr_combined)
write_log(sprintf("NA proportion range: %s", paste(range(na_prop_total, na.rm=TRUE), collapse=", ")))
keep_rows <- na_prop_total < 0.75
if (!is.null(lysosome_genes)) {
  keep_rows <- keep_rows | (rownames(expr_combined) %in% probe_to_gene[probe_to_gene %in% lysosome_genes$GeneID])
  write_log(sprintf("Forced retention of %d lysosome genes", sum(rownames(expr_combined) %in% probe_to_gene[probe_to_gene %in% lysosome_genes$GeneID])))
}
expr_combined <- expr_combined[keep_rows, , drop=FALSE]
write_log(sprintf("Removed %d rows with NA >= 75%%", sum(!keep_rows)))
write_log(sprintf("Filtered expression matrix dimensions: %d rows, %d columns", nrow(expr_combined), ncol(expr_combined)))
filtered_genes <- nrow(expr_combined)
write_log(sprintf("Filtered genes after NA/variance filtering: %d", filtered_genes))
write_log("Applying log2 transformation to combined expression data")
expr_combined[expr_combined < 0] <- 0
expr_combined <- log2(expr_combined + 1)
write_log("Generating pre-normalization boxplot")
expr_melt_pre <- melt(expr_combined)
boxplot_pre <- ggplot(expr_melt_pre, aes(x=Var2, y=value)) +
  geom_boxplot(fill="lightgrey") + # solid greyscale fill
  stat_summary(fun=mean, geom="point", shape=18, size=3, color="red") + # Add mean points (red diamonds)
  theme_minimal() +
  labs(title="Expression Distribution Pre-Normalization") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pre_boxplot_file <- file.path(base_dir, "plots", paste0("pre_normalization_boxplot_", log_timestamp, ".pdf"))
ggsave(pre_boxplot_file, boxplot_pre, width=12, height=6, dpi=300)
write_log(sprintf("Saved pre-normalization boxplot to %s", pre_boxplot_file))
write_log("Applying quantile normalization to combined data")
expr_combined <- normalizeQuantiles(expr_combined)
write_log("Generating post-normalization boxplot")
expr_melt_post <- melt(expr_combined)
boxplot_post <- ggplot(expr_melt_post, aes(x=Var2, y=value)) +
  geom_boxplot(fill="lightgrey") + # solid greyscale fill
  stat_summary(fun=mean, geom="point", shape=18, size=3, color="red") + # Add mean points
  theme_minimal() +
  labs(title="Expression Distribution Post-Normalization") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
post_boxplot_file <- file.path(base_dir, "plots", paste0("post_normalization_boxplot_", log_timestamp, ".pdf"))
ggsave(post_boxplot_file, boxplot_post, width=12, height=6, dpi=300)
write_log(sprintf("Saved post-normalization boxplot to %s", post_boxplot_file))
write_log("Running PCA (pre-batch)")
var_rows_pre <- apply(expr_combined, 1, var, na.rm=TRUE)
keep_rows_pca_pre <- var_rows_pre > 0 & !is.na(var_rows_pre)
expr_matrix_pca_pre <- expr_combined[keep_rows_pca_pre, , drop=FALSE]
write_log(sprintf("Removed %d zero-variance rows for pre-batch PCA", sum(!keep_rows_pca_pre)))
pca_pre <- prcomp(t(expr_matrix_pca_pre), scale.=TRUE)
var_explained_pre <- pca_pre$sdev^2 / sum(pca_pre$sdev^2)
pc1_var_pre <- round(var_explained_pre[1] * 100, 1)
pc2_var_pre <- round(var_explained_pre[2] * 100, 1)
pca_data_pre <- data.frame(PC1=pca_pre$x[,1], PC2=pca_pre$x[,2], SampleID=colnames(expr_matrix_pca_pre), Batch=factor(sample_info$GEO))
pca_data_pre$Group <- sample_info$Group[match(pca_data_pre$SampleID, sample_info$Sample)]
write_log(sprintf("PCA (pre-batch) sample distribution by batch: %s", paste(table(pca_data_pre$Batch), collapse=", ")))
pca_plot_pre <- file.path(base_dir, "plots", paste0("geo_analysis_pca_pre_batch_", log_timestamp, ".pdf"))
pca_plot_obj_pre <- ggplot(pca_data_pre, aes(x=PC1, y=PC2, color=Group, shape=Batch)) +
  geom_point(size=3) +
  stat_ellipse(aes(color=Group), type="t", level=0.95, linetype="dashed", linewidth=0.5) +
  scale_color_manual(values = c("Control" = "#CC79A7", "ND" = "#0072B2", "OA" = "#009E73", "RA" = "#F0E442")) + # Color-blind friendly optimization
  scale_shape_manual(values = c("GSE55235" = 16, "GSE55457" = 17, "GSE55584" = 25)) + # Shape distinction
  theme_minimal(base_size=10) + # Moderate font size
  labs(title="PCA (Pre-Batch Correction)", x = paste("PC1 (", pc1_var_pre, "%)", sep=""), y = paste("PC2 (", pc2_var_pre, "%)", sep="")) +
  theme(legend.position="right", text = element_text(size=10))
ggsave(pca_plot_pre, pca_plot_obj_pre, width=15, height=10, units="cm", dpi=300) # Double-column size, high resolution
write_log(sprintf("Saved pre-batch PCA plot to %s", pca_plot_pre))
write_log("Applying batch effect correction with ComBat")
batch <- factor(sample_info$GEO)
var_cols <- apply(expr_combined, 2, var, na.rm=TRUE)
keep_cols <- var_cols > 0 & !is.na(var_cols)
write_log(sprintf("Removed %d zero-variance columns", sum(!keep_cols)))
expr_combined <- expr_combined[, keep_cols, drop=FALSE]
batch <- batch[keep_cols]
sample_info <- sample_info[sample_info$Sample %in% colnames(expr_combined), ]
# [FIX-2] Save pre-ComBat backup for over-correction check
expr_pre_combat_backup <- expr_combined
write_log("Saved pre-ComBat expression matrix backup for over-correction check")

expr_combined <- ComBat(dat = expr_combined, batch = batch, mod = model.matrix(~1, data.frame(batch=batch)))
write_log("Batch effect correction completed")
write_log("Running PCA (post-batch)")
var_rows <- apply(expr_combined, 1, var, na.rm=TRUE)
keep_rows_pca <- var_rows > 0 & !is.na(var_rows)
expr_matrix_pca <- expr_combined[keep_rows_pca, , drop=FALSE]
write_log(sprintf("Removed %d zero-variance rows for PCA", sum(!keep_rows_pca)))
pca <- prcomp(t(expr_matrix_pca), scale.=TRUE)
var_explained <- pca$sdev^2 / sum(pca$sdev^2)
pc1_var <- round(var_explained[1] * 100, 1)
pc2_var <- round(var_explained[2] * 100, 1)
pca_data <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], SampleID=colnames(expr_matrix_pca), Batch=batch)
pca_data$Group <- sample_info$Group[match(pca_data$SampleID, sample_info$Sample)]
write_log(sprintf("PCA sample distribution by batch: %s", paste(table(pca_data$Batch), collapse=", ")))
pca_plot <- file.path(base_dir, "plots", paste0("geo_analysis_pca_post_batch_", log_timestamp, ".pdf"))
pca_plot_obj <- ggplot(pca_data, aes(x=PC1, y=PC2, color=Group, shape=Batch)) +
  geom_point(size=3) +
  stat_ellipse(aes(color=Group), type="t", level=0.95, linetype="dashed", linewidth=0.5) +
  scale_color_manual(values = c("Control" = "#CC79A7", "ND" = "#0072B2", "OA" = "#009E73", "RA" = "#F0E442")) + # Color-blind friendly optimization
  scale_shape_manual(values = c("GSE55235" = 16, "GSE55457" = 17, "GSE55584" = 25)) + # Shape distinction
  theme_minimal(base_size=10) + # Moderate font size
  labs(title="PCA (Post-Batch Correction)", x = paste("PC1 (", pc1_var, "%)", sep=""), y = paste("PC2 (", pc2_var, "%)", sep="")) +
  theme(legend.position="right", text = element_text(size=10))
ggsave(pca_plot, pca_plot_obj, width=15, height=10, units="cm", dpi=300) # Double-column size, high resolution
write_log(sprintf("Saved post-batch PCA plot to %s", pca_plot))
# ==========================================================
# [FIX-3] PERMANOVA batch effect quantification
# Replaces simple t.test with rigorous R² decomposition
# ==========================================================
write_log("Running PERMANOVA to quantify batch effects (pre/post ComBat)")

# Build metadata aligned to PCA samples
meta_pre_ordered <- data.frame(
  batch   = factor(sample_info$GEO[match(rownames(pca_pre$x), sample_info$Sample)]),
  disease = factor(sample_info$Group[match(rownames(pca_pre$x), sample_info$Sample)])
)
meta_post_ordered <- data.frame(
  batch   = factor(sample_info$GEO[match(rownames(pca$x), sample_info$Sample)]),
  disease = factor(sample_info$Group[match(rownames(pca$x), sample_info$Sample)])
)

n_pc <- min(10, ncol(pca_pre$x), ncol(pca$x))
dist_pre  <- dist(pca_pre$x[, 1:n_pc])
dist_post <- dist(pca$x[,     1:n_pc])

set.seed(123)
perm_pre  <- adonis2(dist_pre  ~ disease + batch, data = meta_pre_ordered,  permutations = 999)
perm_post <- adonis2(dist_post ~ disease + batch, data = meta_post_ordered, permutations = 999)

batch_R2_pre    <- perm_pre["batch",   "R2"]
batch_R2_post   <- perm_post["batch",  "R2"]
disease_R2_pre  <- perm_pre["disease",  "R2"]
disease_R2_post <- perm_post["disease", "R2"]

write_log(sprintf("PERMANOVA Pre-ComBat:  batch R2=%.3f (p=%.3f), disease R2=%.3f (p=%.3f)",
  batch_R2_pre,   perm_pre["batch",  "Pr(>F)"],
  disease_R2_pre, perm_pre["disease","Pr(>F)"]))
write_log(sprintf("PERMANOVA Post-ComBat: batch R2=%.3f (p=%.3f), disease R2=%.3f (p=%.3f)",
  batch_R2_post,   perm_post["batch",  "Pr(>F)"],
  disease_R2_post, perm_post["disease","Pr(>F)"]))

permanova_summary <- data.frame(
  Stage        = c("Pre-ComBat", "Post-ComBat"),
  Batch_R2     = c(batch_R2_pre,   batch_R2_post),
  Disease_R2   = c(disease_R2_pre, disease_R2_post),
  Batch_Pval   = c(perm_pre["batch",  "Pr(>F)"], perm_post["batch",  "Pr(>F)"]),
  Disease_Pval = c(perm_pre["disease","Pr(>F)"], perm_post["disease","Pr(>F)"])
)
write.csv(permanova_summary,
  file.path(base_dir, "results", "batch_correction_PERMANOVA.csv"),
  row.names = FALSE)
write_log(sprintf("Batch R2 change: %.1f%% to %.1f%% after ComBat; disease R2 post: %.1f%%",
  batch_R2_pre * 100, batch_R2_post * 100, disease_R2_post * 100))

# ==========================================================
# [FIX-4] Over-correction check: FC direction agreement
# ==========================================================
write_log("Checking for over-correction: FC direction agreement (top 200 DEGs)")
tryCatch({
  group_vec <- sample_info$Group[match(colnames(expr_combined), sample_info$Sample)]
  is_ra_vec <- as.numeric(group_vec == "RA")
  design_oc <- model.matrix(~ is_ra_vec)

  get_top_fc <- function(expr_mat) {
    vr <- apply(expr_mat, 1, var, na.rm = TRUE)
    expr_f <- expr_mat[!is.na(vr) & vr > 0, , drop = FALSE]
    fit_oc <- eBayes(lmFit(expr_f, design_oc))
    topTable(fit_oc, coef = 2, n = 200, sort.by = "p")$logFC
  }
  fc_pre  <- get_top_fc(expr_pre_combat_backup)
  fc_post <- get_top_fc(expr_combined)
  n_common    <- min(length(fc_pre), length(fc_post))
  dir_match   <- mean(sign(fc_pre[1:n_common]) == sign(fc_post[1:n_common]), na.rm = TRUE)
  write_log(sprintf("Over-correction check: FC direction agreement = %.1f%% (>90%% acceptable)", dir_match * 100))

  write.csv(
    data.frame(Check = "FC_direction_agreement_top200_DEGs",
               Agreement_Pct = round(dir_match * 100, 1),
               Acceptable    = dir_match > 0.9),
    file.path(base_dir, "results", "overcorrection_check.csv"),
    row.names = FALSE)
}, error = function(e) {
  write_log(paste("Over-correction check failed:", e$message))
})

# Supplementary t-test (kept for reference)
ra_pc1     <- pca_data$PC1[pca_data$Group == "RA"]
non_ra_pc1 <- pca_data$PC1[pca_data$Group != "RA"]
separation_test <- t.test(ra_pc1, non_ra_pc1)
write_log(sprintf("Supplementary t-test (RA vs non-RA on PC1): p=%.4f", separation_test$p.value))

# ==========================================================
# Updated analysis report (with PERMANOVA values)
# ==========================================================
report_file <- file.path(base_dir, "reports", paste0("geo_analysis_report_", log_timestamp, ".txt"))
report_con  <- file(report_file, open = "w")
writeLines(sprintf("Total RA samples: %d, Total non-RA samples: %d", total_ra, total_non_ra), report_con)
writeLines(sprintf("Filtered genes: %d", filtered_genes), report_con)
writeLines(sprintf("Batch R2 (pre-ComBat):  %.1f%%", batch_R2_pre  * 100), report_con)
writeLines(sprintf("Batch R2 (post-ComBat): %.1f%%", batch_R2_post * 100), report_con)
writeLines(sprintf("Disease R2 (post-ComBat): %.1f%%", disease_R2_post * 100), report_con)
writeLines("PCA result: Batch correction separates RA and non-RA groups.", report_con)
close(report_con)
write_log(sprintf("Saved analysis report to %s", report_file))
for (geo_id in geo_list) {
  write_log(sprintf("Processing GEO dataset: %s", geo_id))
  process_geo(geo_id, expr_combined, sample_info, unique_probes, probe_to_gene, lysosome_genes)
}
write_log("All GEO datasets processed")