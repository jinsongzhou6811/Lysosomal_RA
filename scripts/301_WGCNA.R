# Load necessary R packages
library(affy)
library(WGCNA)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(ggplot2)
library(pathview)
library(hgu133a.db)
library(pheatmap)
library(R.utils)
library(tools)
library(dplyr)
library(limma)
library(sva)
library(UpSetR)
library(VennDiagram)
library(gplots)
library(pROC)
library(GSEABase)
library(fgsea)
library(KEGGREST)
library(AnnotationDbi)
library(hgu133acdf)
library(msigdbr)

# ── Checkpoint  ───────
save_checkpoint <- function(obj, file_path, log_file_path, description) {
  saveRDS(obj, file_path)
  log_file(paste("✓ saved checkpoint:", description, "→", basename(file_path)), log_file_path)
}

load_checkpoint <- function(file_path, log_file_path, description) {
  obj <- readRDS(file_path)
  log_file(paste("✓ loaded checkpoint:", description, "from", basename(file_path)), log_file_path)
  return(obj)
}

# ========================================================
# Standard Path Management Template
# ========================================================
library(rprojroot)
project_root <- find_rstudio_root_file()
data_dir     <- file.path(project_root, "data")
outputs_dir  <- file.path(project_root, "outputs")
module_name  <- "03_WGCNA"
module_dir   <- file.path(outputs_dir, module_name)
dir.create(module_dir, showWarnings = FALSE, recursive = TRUE)
base_dir <- module_dir

# ── Modifiable parameters ────────────────────────────────────────────────────
min_cluster_size    <- 30
soft_power_default  <- 6
high_cor_threshold  <- 0.48
high_cor_pval       <- 0.05
core_mm_threshold   <- 0.5
core_gs_threshold   <- 0.2
tom_threshold       <- 0.1
de_logfc_threshold  <- 1
de_padj_threshold   <- 0.05
enrich_qvalue_cutoff <- 1
enrich_padj_cutoff  <- 0.05
gsea_min_size       <- 15
gsea_max_size       <- 500
gsea_pvalue_cutoff  <- 0.05
hub_top_percent     <- 0.1

# ── Log helper ───────────────────────────────────────────────────────────────
log_file <- function(message, file_path) {
  cat(paste0(Sys.time(), ": ", message, "\n"), file = file_path, append = TRUE)
  cat(paste0(Sys.time(), ": ", message, "\n"))
}

# ── [STANDARDISED] Triple-format figure export helper ────────────────────────
# Outputs: .pdf (vector), _300.png (submission), _600.tiff (production)
save_fig <- function(p, stem, w_in, h_in) {
  ggsave(paste0(stem, ".pdf"), plot = p, width = w_in, height = h_in,
         units = "in", device = "pdf", useDingbats = FALSE)
  ggsave(paste0(stem, "_300.png"), plot = p, width = w_in, height = h_in,
         units = "in", device = "png", dpi = 300)
  ggsave(paste0(stem, "_600.tiff"), plot = p, width = w_in, height = h_in,
         units = "in", device = "tiff", dpi = 600, compression = "lzw")
}

# Helper for base R plots (pheatmap, corrplot, base plot, etc.)
save_base_fig <- function(plot_expr, stem, w_in, h_in) {
  pdf(paste0(stem, ".pdf"), width = w_in, height = h_in)
  eval(plot_expr); dev.off()
  png(paste0(stem, "_300.png"), width = w_in, height = h_in, units = "in", res = 300)
  eval(plot_expr); dev.off()
  tiff(paste0(stem, "_600.tiff"), width = w_in, height = h_in, units = "in", res = 600, compression = "lzw")
  eval(plot_expr); dev.off()
}

# ── Dataset definitions ──────────────────────────────────────────────────────
datasets <- list(
  GSE55235 = list(
    expected_samples = paste0("GSM13322", sprintf("%02d", 1:30)),
    sample_counts    = c(ND = 10, RA = 10, OA = 10),
    dataset_groups   = c("ND", "RA", "OA")
  ),
  GSE55457 = list(
    expected_samples = paste0("GSM133", 7304:7336),
    sample_counts    = c(ND = 10, RA = 13, OA = 10),
    dataset_groups   = c("ND", "RA", "OA")
  ),
  GSE55584 = list(
    expected_samples = c(paste0("GSM13396", 18:27), paste0("GSM13396", 28:33)),
    sample_counts    = c(ND = 0, RA = 10, OA = 6),
    dataset_groups   = c("RA", "OA")
  )
)

seed <- 123

target_genes_file <- file.path(data_dir, "lysosome_target_genes.txt")
if (!file.exists(target_genes_file)) stop(paste("Target genes file not found:", target_genes_file))
target_genes <- readLines(target_genes_file)
target_genes <- trimws(target_genes[target_genes != ""])
if (length(target_genes) == 0) stop("Target genes file is empty")
core_lysosome_genes <- unique(target_genes)

global_log <- file.path(base_dir, "global_log.txt")
log_file(paste("Loaded", length(target_genes), "target genes from file:", target_genes_file), global_log)
log_file(paste("Core lysosome genes:", paste(core_lysosome_genes, collapse = ", ")), global_log)

# ============================================================
# Main loop over datasets
# ============================================================
for (dataset_id in names(datasets)) {
  ds               <- datasets[[dataset_id]]
  expected_samples <- ds$expected_samples
  sample_counts    <- ds$sample_counts
  dataset_groups   <- ds$dataset_groups

  results_dir <- file.path(base_dir, dataset_id, "Results_WGCNA_GO-KEGG-Reactome")
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  main_log <- file.path(results_dir, "analysis_log.txt")
  log_file(paste("Starting WGCNA analysis for", dataset_id, "with seed", seed), main_log)

  temp_dir <- file.path(results_dir, "temp")
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
  Sys.setenv(TMPDIR = temp_dir)
  log_file(paste("Set temporary directory to:", temp_dir), main_log)

  set.seed(seed)
  log_file(paste("Random seed set to", seed, "for reproducibility"), main_log)

  # ── Step 1: Read & normalise CEL files ─────────────────────────────────────
  cel_dir      <- file.path(data_dir, paste0(dataset_id, "_RAW"))
  log_file("Starting to process CEL files. This may take a few minutes...", main_log)
  cel_gz_files <- list.files(cel_dir, pattern = "\\.[cC][eE][lL]\\.gz$", full.names = TRUE)
  log_file(paste("Found", length(cel_gz_files), ".cel.gz files"), main_log)

  if (length(cel_gz_files) > 0) {
    for (gz_file in cel_gz_files) {
      cel_file <- sub("\\.gz$", "", gz_file)
      tryCatch({
        gunzip(gz_file, destname = cel_file, overwrite = TRUE, remove = FALSE)
        log_file(paste("Decompressed:", cel_file), main_log)
      }, error = function(e) {
        log_file(paste("Error decompressing", gz_file, ":", e$message), main_log)
      })
    }
  }

  cel_files <- list.files(cel_dir, pattern = "\\.[cC][eE][lL]$", full.names = TRUE)
  log_file(paste("Found", length(cel_files), ".CEL files"), main_log)
  if (length(cel_files) == 0) stop("No .cel or .CEL files found in directory: ", cel_dir)

  log_file("Reading Affymetrix data", main_log)
  raw_data <- ReadAffy(filenames = cel_files)
  log_file("Performing RMA normalization", main_log)
  eset     <- rma(raw_data)
  exprData <- t(exprs(eset))
  log_file(paste("Expression data dimensions:", nrow(exprData), "samples,", ncol(exprData), "probes"), main_log)

  # ── Step 2: Phenotype data ─────────────────────────────────────────────────
  log_file("Defining phenoData", main_log)
  pheno_status <- c()
  for (group in names(sample_counts)) {
    if (sample_counts[group] > 0) pheno_status <- c(pheno_status, rep(group, sample_counts[group]))
  }
  phenoData            <- data.frame(RA_status = pheno_status, row.names = expected_samples)
  phenoData$RA_status  <- as.factor(phenoData$RA_status)
  log_file(paste("RA_status levels:", paste(levels(phenoData$RA_status), collapse = ", ")), main_log)

  if (!is.data.frame(phenoData) || !"RA_status" %in% colnames(phenoData))
    stop("phenoData is not a data frame or RA_status column not found")

  rownames(exprData) <- sub("^(GSM\\d+)_.*", "\\1", rownames(exprData))
  log_file("Sample IDs in exprData:", main_log)
  log_file(paste(rownames(exprData), collapse = ", "), main_log)

  log_file("Checking sample ID consistency", main_log)
  missing_in_exprData  <- expected_samples[!expected_samples %in% rownames(exprData)]
  missing_in_phenoData <- rownames(exprData)[!rownames(exprData) %in% expected_samples]
  if (length(missing_in_exprData) > 0 || length(missing_in_phenoData) > 0) {
    log_file("Adjusting phenoData to match exprData samples", main_log)
    phenoData <- phenoData[rownames(phenoData) %in% rownames(exprData), , drop = FALSE]
    log_file(paste("Adjusted phenoData to", nrow(phenoData), "samples"), main_log)
  } else {
    log_file("All sample IDs match between exprData and phenoData", main_log)
  }

  log_file("Reordering exprData to match phenoData", main_log)
  match_indices <- match(rownames(phenoData), rownames(exprData))
  if (any(is.na(match_indices))) stop("Some samples could not be matched to phenoData.")
  exprData <- exprData[match_indices, ]
  log_file("exprData reordered to match phenoData", main_log)

  log_file("Checking for NA values in exprData", main_log)
  if (any(is.na(exprData))) {
    log_file("Warning: exprData contains NA values, removing affected genes", main_log)
    exprData <- exprData[, !apply(is.na(exprData), 2, any)]
    log_file(paste("exprData dimensions after NA removal:", nrow(exprData), "samples,", ncol(exprData), "probes"), main_log)
  }

  # ── Step 3: PCA ────────────────────────────────────────────────────────────
  log_file("Checking normalization with PCA", main_log)
  pca_pre <- prcomp(t(exprData), scale. = TRUE)

  # [MODIFIED] Triple-format output for PCA pre-batch
  save_base_fig(quote({
    plot(pca_pre$x[,1], pca_pre$x[,2], col = as.numeric(phenoData$RA_status),
         main = "PCA of RMA Normalized Data", xlab = "PC1", ylab = "PC2", pch = 16)
    legend("topright", legend = levels(phenoData$RA_status),
           col = 1:length(levels(phenoData$RA_status)), pch = 16)
  }), file.path(results_dir, "pca_pre_batch_correction"), w_in = 8, h_in = 6)
  log_file("PCA plot before batch correction generated (.pdf/_300.png/_600.tiff)", main_log)

  log_file("Checking for NA values in phenoData$RA_status", main_log)
  if (any(is.na(phenoData$RA_status))) stop("phenoData$RA_status contains NA values")

  log_file("No obvious batch effect observed, skipping batch correction", main_log)

  pca_post <- prcomp(t(exprData), scale. = TRUE)
  # [MODIFIED] Triple-format output for PCA post-batch
  save_base_fig(quote({
    plot(pca_post$x[,1], pca_post$x[,2], col = as.numeric(phenoData$RA_status),
         main = "PCA after No Batch Correction", xlab = "PC1", ylab = "PC2", pch = 16)
    legend("topright", legend = levels(phenoData$RA_status),
           col = 1:length(levels(phenoData$RA_status)), pch = 16)
  }), file.path(results_dir, "pca_post_no_batch_correction"), w_in = 8, h_in = 6)
  log_file("PCA plot after no batch correction generated (.pdf/_300.png/_600.tiff)", main_log)

  # ── Step 4: Probe → gene symbol mapping ────────────────────────────────────
  probe_ids    <- colnames(exprData)
  log_file("Mapping probes to gene symbols", main_log)
  gene_symbols <- mapIds(hgu133a.db, keys = probe_ids, column = "SYMBOL",
                         keytype = "PROBEID", multiVals = "first")
  probe_to_gene <- data.frame(ProbeID = probe_ids, OriginalSymbol = gene_symbols)
  if (any(is.na(probe_to_gene$OriginalSymbol))) {
    log_file(paste(sum(is.na(probe_to_gene$OriginalSymbol)),
                   "probes could not be mapped to gene symbols"), main_log)
    probe_to_gene$OriginalSymbol[is.na(probe_to_gene$OriginalSymbol)] <-
      probe_to_gene$ProbeID[is.na(probe_to_gene$OriginalSymbol)]
  }

  log_file("Aggregating probe expression by mean", main_log)
  exprData_agg         <- aggregate(t(exprData),
                                    by  = list(Group.1 = probe_to_gene$OriginalSymbol),
                                    FUN = mean)
  rownames(exprData_agg) <- exprData_agg$Group.1
  exprData_agg           <- t(exprData_agg[, -1])
  colnames(exprData_agg) <- make.unique(colnames(exprData_agg))
  exprData               <- exprData_agg
  gene_symbol_map        <- data.frame(ColumnName     = colnames(exprData),
                                       OriginalSymbol = colnames(exprData))
  log_file(paste("Expression data after probe aggregation:", ncol(exprData), "genes"), main_log)

  # ── Step 5: goodSamplesGenes ────────────────────────────────────────────────
  log_file("Performing goodSamplesGenes check", main_log)
  gsg <- goodSamplesGenes(exprData, verbose = 3)
  sink(file.path(results_dir, "preprocessing_results.txt"))
  print("Good samples and genes check:"); print(gsg)
  sink()
  file.copy(file.path(results_dir, "preprocessing_results.txt"),
            file.path(results_dir, "preprocessing_results.csv"))
  log_file(paste("Good samples:", sum(gsg$goodSamples), "Good genes:", sum(gsg$goodGenes)), main_log)
  if (!gsg$allOK) exprData <- exprData[gsg$goodSamples, gsg$goodGenes]

  # ── Step 6: Soft threshold (serial, no parallel, no checkpoint) ─────────────
  # NOTE: registerDoParallel removed — it hangs indefinitely on Windows.
  # WGCNA::allowWGCNAThreads() uses OpenMP threads instead, which are safe.
  WGCNA::allowWGCNAThreads()

  powers <- c(1:10, seq(12, 20, 2))
  log_file("Selecting soft threshold (serial)", main_log)
  sft <- pickSoftThreshold(exprData, powerVector = powers, verbose = 0)

  write.csv(sft$fitIndices,
            file.path(results_dir, "soft_threshold_results.csv"),
            row.names = FALSE)

  softPower <- sft$powerEstimate
  if (is.na(softPower)) {
    softPower <- soft_power_default
    log_file("Soft power is NA, defaulting to 6", main_log)
  }
  log_file(paste("Selected soft power:", softPower), main_log)

  # Scale-free R² summary
  tryCatch({
    fi  <- sft$fitIndices
    idx <- which(fi[, 1] == softPower)
    if (length(idx) > 0) {
      r2_col <- grep("SFT.R.sq|R.sq", colnames(fi), value = TRUE)[1]
      mk_col <- grep("mean.k",         colnames(fi), value = TRUE)[1]
      sft_summary <- data.frame(
        Dataset           = dataset_id,
        Soft_Power        = softPower,
        Scale_Free_R2     = round(abs(fi[idx, r2_col]), 3),
        Mean_Connectivity = round(fi[idx, mk_col], 1)
      )
      write.csv(sft_summary,
                file.path(results_dir, "sft_selected_power_summary.csv"),
                row.names = FALSE)
      log_file(sprintf("Dataset %s: soft_power=%d, scale-free R2=%.3f, mean.k=%.1f",
                       dataset_id, softPower,
                       sft_summary$Scale_Free_R2, sft_summary$Mean_Connectivity), main_log)
    }
  }, error = function(e) log_file(paste("SFT summary extraction failed:", e$message), main_log))

  # ── Step 7: Adjacency → TOM (no checkpoint) ─────────────────────────────────
  log_file("Validating exprData colnames before adjacency matrix construction", main_log)
  if (any(is.na(colnames(exprData))) || any(duplicated(colnames(exprData))))
    stop("Invalid colnames in exprData")

  log_file("Building adjacency matrix", main_log)
  adjacency           <- adjacency(exprData, power = softPower)
  dimnames(adjacency) <- list(colnames(exprData), colnames(exprData))

  log_file("Computing TOM similarity", main_log)
  TOM           <- TOMsimilarity(adjacency)
  dimnames(TOM) <- dimnames(adjacency)
  dissTOM       <- 1 - TOM
  log_file(paste("TOM matrix dimensions:", nrow(TOM), "x", ncol(TOM)), main_log)

  if (nrow(dissTOM) < 2 || ncol(dissTOM) < 2)
    stop("TOM distance matrix has fewer than 2 genes")

  # ── Step 8: Module detection ────────────────────────────────────────────────
  log_file("Performing hierarchical clustering", main_log)
  geneTree    <- hclust(as.dist(dissTOM), method = "average")
  log_file("Identifying dynamic modules", main_log)
  dynamicMods <- cutreeDynamic(dendro    = geneTree,
                               distM     = dissTOM,
                               minClusterSize = min_cluster_size)
  moduleColors <- labels2colors(dynamicMods)

  module_sizes    <- table(moduleColors)
  grey_proportion <- module_sizes["grey"] / sum(module_sizes)
  log_file(paste("Module sizes:", paste(names(module_sizes), module_sizes,
                                        sep = "=", collapse = ", ")), main_log)
  log_file(paste("Grey module proportion:", grey_proportion), main_log)
  if (grey_proportion > 0.5)
    log_file("Warning: Grey module > 50% — consider adjusting minClusterSize", main_log)
  log_file(paste("Identified", length(unique(moduleColors)), "modules"), main_log)

  moduleColors_df <- data.frame(Gene = colnames(exprData), ModuleColor = moduleColors)
  write.csv(moduleColors_df, file.path(results_dir, "moduleColors.txt"),  row.names = FALSE)
  write.csv(moduleColors_df, file.path(results_dir, "moduleColors.csv"),  row.names = FALSE)
  log_file("moduleColors.txt and .csv saved successfully", main_log)

  # ── Step 9: Module eigengenes ───────────────────────────────────────────────
  log_file("Calculating module eigengenes", main_log)
  valid_modules <- names(module_sizes)[module_sizes >= 2]
  if (length(valid_modules) == 0) stop("No modules with >= 2 genes found")
  MEs <- moduleEigengenes(exprData, colors = moduleColors)$eigengenes
  if (is.null(MEs)) stop("moduleEigengenes failed")
  log_file(paste("Module eigengenes calculated for", ncol(MEs), "modules"), main_log)

  module_gene_counts <- data.frame(Module    = names(module_sizes),
                                   GeneCount = as.numeric(module_sizes))
  write.csv(module_gene_counts, file.path(results_dir, "module_gene_counts.txt"), row.names = FALSE)
  write.csv(module_gene_counts, file.path(results_dir, "module_gene_counts.csv"), row.names = FALSE)
  log_file("module_gene_counts saved", main_log)

  # ── Step 10: Module-trait correlations ──────────────────────────────────────
  log_file("Calculating module-trait correlations for RA_status", main_log)
  pheno_numeric <- ifelse(phenoData$RA_status == "RA", 1, 0)
  if (length(pheno_numeric) != nrow(exprData))
    stop("Length of pheno_numeric does not match nrow(exprData)")

  ME_pheno_cor  <- cor(MEs, pheno_numeric, use = "p")
  ME_pheno_pval <- corPvalueStudent(ME_pheno_cor, nrow(exprData))

  module_trait_results <- data.frame(Module      = rownames(ME_pheno_cor),
                                     Correlation = ME_pheno_cor[, 1],
                                     P_value     = ME_pheno_pval[, 1])
  write.csv(module_trait_results,
            file.path(results_dir, "module_trait_results.csv"), row.names = FALSE)

  sink(file.path(results_dir, "module_trait_correlation.txt"))
  print("Module-trait correlations:"); print(ME_pheno_cor)
  print("Module-trait p-values:");     print(ME_pheno_pval)
  sink()
  file.copy(file.path(results_dir, "module_trait_correlation.txt"),
            file.path(results_dir, "module_trait_correlation.csv"))
  log_file("Module-trait correlations saved", main_log)

  # Heatmap
  log_file("Generating module-trait heatmap", main_log)
  module_trait_matrix <- as.matrix(ME_pheno_cor)
  pval_matrix         <- as.matrix(ME_pheno_pval)

  # [MODIFIED] Triple-format output for module-trait heatmap
  save_base_fig(quote({
    pheatmap(module_trait_matrix,
             display_numbers = TRUE, number_format = "%.2f", number_color = "black",
             fontsize_number = 8, cellwidth = 40, cellheight = 20,
             main  = "Module-Trait Correlation Heatmap",
             color = colorRampPalette(c("blue", "white", "red"))(50),
             breaks = seq(-1, 1, length.out = 51),
             cluster_rows = TRUE,
             cluster_cols = ncol(module_trait_matrix) > 1)
  }), file.path(results_dir, "module_trait_heatmap"), w_in = 8, h_in = 6)
  write.csv(data.frame(Module      = rownames(module_trait_matrix),
                       Correlation = module_trait_matrix[, 1],
                       P_value     = pval_matrix[, 1]),
            file.path(results_dir, "module_trait_heatmap_data.csv"), row.names = FALSE)
  log_file("Module-trait heatmap generated", main_log)

  # High-correlation modules
  high_cor_modules <- sub("ME", "", rownames(ME_pheno_cor)[
    abs(ME_pheno_cor[, 1]) > high_cor_threshold &
      ME_pheno_pval[, 1] < high_cor_pval])
  write.csv(data.frame(HighCorModules = high_cor_modules),
            file.path(results_dir, "high_cor_modules.csv"), row.names = FALSE)
  log_file(paste("Identified", length(high_cor_modules), "high-correlation modules"), main_log)

  # ── Step 11: Module membership & gene significance ──────────────────────────
  log_file("Calculating module membership and gene significance", main_log)
  geneModuleMembership <- tryCatch(
    cor(exprData, MEs, use = "p"),
    error = function(e) { log_file(paste("Error:", e$message), main_log); stop(e) }
  )
  log_file(paste("geneModuleMembership dimensions:",
                 nrow(geneModuleMembership), "rows,",
                 ncol(geneModuleMembership), "columns"), main_log)

  MMPvalue <- corPvalueStudent(as.matrix(geneModuleMembership), nrow(exprData))

  geneTraitSignificance <- tryCatch(
    cor(exprData, pheno_numeric, use = "p"),
    error = function(e) { log_file(paste("Error:", e$message), main_log); stop(e) }
  )
  GSPvalue <- corPvalueStudent(as.matrix(geneTraitSignificance), nrow(exprData))
  log_file("Module membership and gene significance calculated", main_log)

  # ── Step 12: GS RA vs ND ───────────────────────────────────────────────────
  keep_samples_RA_ND <- phenoData$RA_status %in% c("ND", "RA")
  if (sum(keep_samples_RA_ND) < 2) {
    log_file("Error: Fewer than 2 samples for RA vs ND analysis", main_log)
    geneTraitSignificance_RA_ND <- NULL; GSPvalue_RA_ND <- NULL
  } else {
    exprData_RA_ND  <- exprData[keep_samples_RA_ND, , drop = FALSE]
    phenoData_RA_ND <- phenoData[keep_samples_RA_ND, , drop = FALSE]
    phenoData_RA_ND$RA_status <- factor(phenoData_RA_ND$RA_status, levels = c("ND","RA"))
    datTraits_RA_ND <- as.numeric(phenoData_RA_ND$RA_status) - 1
    geneTraitSignificance_RA_ND <- tryCatch(
      cor(exprData_RA_ND, datTraits_RA_ND, use = "p"),
      error = function(e) { log_file(paste("Error:", e$message), main_log); stop(e) }
    )
    GSPvalue_RA_ND <- corPvalueStudent(as.matrix(geneTraitSignificance_RA_ND),
                                       nrow(exprData_RA_ND))
    write.csv(GSPvalue_RA_ND,
              file.path(results_dir, "wgcna_gs_pvalues_RA_vs_ND.csv"), row.names = TRUE)
    log_file("Gene significance for RA vs ND saved", main_log)
  }

  # ── Step 13: GS RA vs OA ───────────────────────────────────────────────────
  keep_samples_RA_OA <- phenoData$RA_status %in% c("OA", "RA")
  if (sum(keep_samples_RA_OA) < 2) {
    log_file("Error: Fewer than 2 samples for RA vs OA analysis", main_log)
    geneTraitSignificance_RA_OA <- NULL; GSPvalue_RA_OA <- NULL
  } else {
    exprData_RA_OA  <- exprData[keep_samples_RA_OA, , drop = FALSE]
    phenoData_RA_OA <- phenoData[keep_samples_RA_OA, , drop = FALSE]
    phenoData_RA_OA$RA_status <- factor(phenoData_RA_OA$RA_status, levels = c("OA","RA"))
    datTraits_RA_OA <- as.numeric(phenoData_RA_OA$RA_status) - 1
    geneTraitSignificance_RA_OA <- tryCatch(
      cor(exprData_RA_OA, datTraits_RA_OA, use = "p"),
      error = function(e) { log_file(paste("Error:", e$message), main_log); stop(e) }
    )
    GSPvalue_RA_OA <- corPvalueStudent(as.matrix(geneTraitSignificance_RA_OA),
                                       nrow(exprData_RA_OA))
    write.csv(GSPvalue_RA_OA,
              file.path(results_dir, "wgcna_gs_pvalues_RA_vs_OA.csv"), row.names = TRUE)
    log_file("Gene significance for RA vs OA saved", main_log)
  }

  # ── Step 14: Lysosome gene list ─────────────────────────────────────────────
  lysosome_genes_file <- file.path(data_dir, "lysosome_genes.csv")
  if (!file.exists(lysosome_genes_file)) stop("lysosome_genes.csv not found")
  lysosome_df    <- read.csv(lysosome_genes_file)
  lysosome_genes <- unique(c("GZMB","MREG","PRF1","LAMP3","NKG7",
                              "SLC39A8","TRAF3IP3","VOPP1",
                              lysosome_df$Symbol))

  # ── Step 15: Differential expression ───────────────────────────────────────
  log_file("Running differential expression analysis", main_log)
  lysosome_cols       <- gene_symbol_map$ColumnName[gene_symbol_map$OriginalSymbol %in% lysosome_genes]
  valid_lysosome_cols <- lysosome_cols[lysosome_cols %in% colnames(exprData)]
  lysosome_expr       <- exprData[, valid_lysosome_cols, drop = FALSE]
  log_file(paste("Lysosome expression data:", nrow(lysosome_expr),
                 "samples,", ncol(lysosome_expr), "genes"), main_log)
  if (any(is.na(lysosome_expr))) stop("lysosome_expr contains NA values")

  # RA vs ND
  keep_samples <- phenoData$RA_status %in% c("ND","RA")
  if (sum(keep_samples) < 2) {
    log_file("Error: <2 samples for RA vs ND DE", main_log)
    top_table_RA_ND                <- data.frame()
    significant_lysosome_genes_RA_ND <- data.frame()
  } else {
    lysosome_expr_RA_ND  <- lysosome_expr[keep_samples, , drop = FALSE]
    phenoData_binary     <- phenoData[keep_samples, , drop = FALSE]
    phenoData_binary$RA_status <- factor(phenoData_binary$RA_status, levels = c("ND","RA"))
    log_file(paste("RA vs ND sample counts:",
                   paste(table(phenoData_binary$RA_status), collapse = ", ")), main_log)
    design      <- model.matrix(~ phenoData_binary$RA_status)
    fit         <- lmFit(t(lysosome_expr_RA_ND), design)
    fit         <- eBayes(fit)
    coef_name   <- colnames(design)[grepl("RA_statusRA", colnames(design))]
    top_table_RA_ND <- topTable(fit, coef = coef_name, number = Inf, adjust.method = "BH")
    top_table_RA_ND$GeneSymbol <- gene_symbol_map$OriginalSymbol[
      match(rownames(top_table_RA_ND), gene_symbol_map$ColumnName)]
    significant_lysosome_genes_RA_ND <- top_table_RA_ND[
      abs(top_table_RA_ND$logFC) > de_logfc_threshold &
        top_table_RA_ND$adj.P.Val < de_padj_threshold, ]
    write.csv(top_table_RA_ND,
              file.path(results_dir, "lysosome_diff_expression_RA_vs_ND.csv"), row.names = FALSE)
    write.csv(significant_lysosome_genes_RA_ND,
              file.path(results_dir, "significant_lysosome_genes_RA_vs_ND.csv"), row.names = FALSE)
    log_file(paste("Found", nrow(significant_lysosome_genes_RA_ND),
                   "significant lysosome genes for RA vs ND"), main_log)
  }

  # RA vs OA
  keep_samples_RA_OA2 <- phenoData$RA_status %in% c("OA","RA")
  if (sum(keep_samples_RA_OA2) < 2) {
    log_file("Error: <2 samples for RA vs OA DE", main_log)
    top_table_RA_OA                <- data.frame()
    significant_lysosome_genes_RA_OA <- data.frame()
  } else {
    lysosome_expr_RA_OA2 <- lysosome_expr[keep_samples_RA_OA2, , drop = FALSE]
    phenoData_RA_OA2     <- phenoData[keep_samples_RA_OA2, , drop = FALSE]
    phenoData_RA_OA2$RA_status <- factor(phenoData_RA_OA2$RA_status, levels = c("OA","RA"))
    log_file(paste("RA vs OA sample counts:",
                   paste(table(phenoData_RA_OA2$RA_status), collapse = ", ")), main_log)
    design_RA_OA   <- model.matrix(~ phenoData_RA_OA2$RA_status)
    fit_RA_OA      <- lmFit(t(lysosome_expr_RA_OA2), design_RA_OA)
    fit_RA_OA      <- eBayes(fit_RA_OA)
    coef_name_RA_OA <- colnames(design_RA_OA)[grepl("RA_statusRA", colnames(design_RA_OA))]
    top_table_RA_OA <- topTable(fit_RA_OA, coef = coef_name_RA_OA,
                                number = Inf, adjust.method = "BH")
    top_table_RA_OA$GeneSymbol <- gene_symbol_map$OriginalSymbol[
      match(rownames(top_table_RA_OA), gene_symbol_map$ColumnName)]
    significant_lysosome_genes_RA_OA <- top_table_RA_OA[
      abs(top_table_RA_OA$logFC) > de_logfc_threshold &
        top_table_RA_OA$adj.P.Val < de_padj_threshold, ]
    write.csv(top_table_RA_OA,
              file.path(results_dir, "lysosome_diff_expression_RA_vs_OA.csv"), row.names = FALSE)
    write.csv(significant_lysosome_genes_RA_OA,
              file.path(results_dir, "significant_lysosome_genes_RA_vs_OA.csv"), row.names = FALSE)
    log_file(paste("Found", nrow(significant_lysosome_genes_RA_OA),
                   "significant lysosome genes for RA vs OA"), main_log)
  }

  # ── Step 16: Core lysosome gene enrichment ──────────────────────────────────
  
  entrez_core <- bitr(core_lysosome_genes, fromType = "SYMBOL",
                      toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID
  
  # ── GO enrichment for core lysosome genes ──────────────────────────────────
  go_enrich_core <- tryCatch({
    enrichGO(gene          = core_lysosome_genes,
             OrgDb         = "org.Hs.eg.db",
             keyType       = "SYMBOL",
             ont           = "ALL",
             pAdjustMethod = "BH",
             qvalueCutoff  = 1)
  }, error = function(e) {
    log_file(paste("GO core enrichment error:", e$message), main_log)
    NULL
  })
  go_enrich_core_df <- if (!is.null(go_enrich_core)) as.data.frame(go_enrich_core) else data.frame()
  write.csv(go_enrich_core_df,
            file.path(results_dir, "go_enrichment_core_lysosome.csv"), row.names = FALSE)
  log_file(paste("GO enrichment for core lysosome genes:", nrow(go_enrich_core_df), "terms"), main_log)
  
  # ── Local KEGG enrichment using the downloaded file  ─────────────────────
  kegg_file <- file.path(data_dir, "kegg_jp_link_hsa_pathway.txt")
  if (!file.exists(kegg_file)) stop("Local KEGG file not found: ", kegg_file)
  
  df <- read.table(kegg_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("Pathway", "Gene")
  df$Gene <- gsub("^hsa:", "", df$Gene)
  pathway_list <- split(df$Gene, df$Pathway)
  names(pathway_list) <- gsub("path:", "", names(pathway_list))
  
  term2gene <- data.frame(
    term = rep(names(pathway_list), lengths(pathway_list)),
    gene = unlist(pathway_list)
  )
  
  kegg_enrich_core <- tryCatch({
    enricher(gene          = entrez_core,
             TERM2GENE     = term2gene,
             pvalueCutoff  = 1,
             qvalueCutoff  = 1)
  }, error = function(e) {
    log_file(paste("Local KEGG enrichment error:", e$message), main_log)
    NULL
  })
  
  kegg_enrich_core_df <- if (!is.null(kegg_enrich_core)) as.data.frame(kegg_enrich_core) else data.frame()
  
  if (nrow(kegg_enrich_core_df) > 0) {
    kegg_categories <- lapply(kegg_enrich_core_df$ID, function(id) {
      info <- KEGGREST::keggGet(id)
      cat  <- info[[1]]$CLASS
      if (is.null(cat)) c("Unknown", "Unknown") else strsplit(cat, "; ")[[1]]
    })
    kegg_enrich_core_df$category    <- sapply(kegg_categories, `[`, 1)
    kegg_enrich_core_df$subcategory <- sapply(kegg_categories, `[`, 2)
    kegg_enrich_core_df$RichFactor  <- as.numeric(gsub("/.*", "", kegg_enrich_core_df$GeneRatio)) /
      as.numeric(gsub("/.*", "", kegg_enrich_core_df$BgRatio))
    kegg_enrich_core_df$FoldEnrichment <- kegg_enrich_core_df$RichFactor *
      (as.numeric(gsub(".*/", "", kegg_enrich_core_df$BgRatio)) /
         as.numeric(gsub(".*/", "", kegg_enrich_core_df$GeneRatio)))
    kegg_enrich_core_df$zScore <- sign(kegg_enrich_core_df$FoldEnrichment) *
      sqrt(-2 * log(kegg_enrich_core_df$pvalue))
    kegg_enrich_core_df$geneID <- sapply(kegg_enrich_core_df$geneID, function(x) {
      ids     <- strsplit(x, "/")[[1]]
      symbols <- bitr(ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")$SYMBOL
      paste(symbols, collapse = "/")
    })
    kegg_enrich_core_df <- kegg_enrich_core_df[, c("category", "subcategory", "ID", "Description",
                                                   "GeneRatio", "BgRatio", "RichFactor",
                                                   "FoldEnrichment", "zScore", "pvalue",
                                                   "p.adjust", "qvalue", "geneID", "Count")]
  } else {
    kegg_enrich_core_df <- data.frame(
      category = character(0), subcategory = character(0), ID = character(0),
      Description = character(0), GeneRatio = character(0), BgRatio = character(0),
      RichFactor = numeric(0), FoldEnrichment = numeric(0), zScore = numeric(0),
      pvalue = numeric(0), p.adjust = numeric(0), qvalue = numeric(0),
      geneID = character(0), Count = integer(0)
    )
  }
  write.csv(kegg_enrich_core_df,
            file.path(results_dir, "kegg_enrichment_core_lysosome.csv"), row.names = FALSE)
  log_file(paste("KEGG enrichment for core lysosome genes:", nrow(kegg_enrich_core_df), "pathways"), main_log)
  
  # ── Reactome enrichment for core lysosome genes ────────────────────────────
  reactome_enrich_core <- tryCatch({
    enrichPathway(gene          = entrez_core,
                  organism      = "human",
                  pAdjustMethod = "BH",
                  qvalueCutoff  = 1)
  }, error = function(e) {
    log_file(paste("Reactome core enrichment error:", e$message), main_log)
    NULL
  })
  reactome_enrich_core_df <- if (!is.null(reactome_enrich_core)) as.data.frame(reactome_enrich_core) else data.frame()
  write.csv(reactome_enrich_core_df,
            file.path(results_dir, "reactome_enrichment_core_lysosome.csv"), row.names = FALSE)
  log_file(paste("Reactome enrichment for core lysosome genes:", nrow(reactome_enrich_core_df), "pathways"), main_log)
  
  
  # ── Step 17: Core gene selection ────────────────────────────────────────────
  log_file("Selecting core genes", main_log)
  core_genes <- list()
  log_file(paste("geneModuleMembership columns:", paste(colnames(geneModuleMembership), collapse = ", ")), main_log)
  log_file(paste("MEs columns:", paste(colnames(MEs), collapse = ", ")), main_log)

  for (mod in unique(moduleColors)) {
    mod_genes <- colnames(exprData)[moduleColors == mod]
    if (length(mod_genes) >= 2) {
      mod_col <- paste0("ME", mod)
      if (!mod_col %in% colnames(geneModuleMembership)) {
        log_file(paste("Error: Column", mod_col, "not found in geneModuleMembership"), main_log)
        core_genes[[mod]] <- character(0); next
      }
      valid_mod_genes <- mod_genes[mod_genes %in% rownames(geneModuleMembership)]
      if (length(valid_mod_genes) == 0) {
        log_file(paste("Error: No valid genes for module", mod), main_log)
        core_genes[[mod]] <- character(0); next
      }
      mod_MM <- geneModuleMembership[valid_mod_genes, mod_col]
      mod_GS <- abs(geneTraitSignificance[valid_mod_genes, 1])
      core_genes[[mod]] <- valid_mod_genes[
        abs(mod_MM) > core_mm_threshold & mod_GS > core_gs_threshold &
          !is.na(mod_MM) & !is.na(mod_GS)]
      log_file(paste("Module", mod, ":", length(core_genes[[mod]]), "core genes"), main_log)
    } else {
      log_file(paste("Skipping core gene selection for module", mod, "(<2 genes)"), main_log)
      core_genes[[mod]] <- character(0)
    }
  }
  sink(file.path(results_dir, "core_genes.txt"))
  print(lapply(core_genes, function(x) if (length(x) > 0) unique(x) else "No core genes"))
  sink()
  file.copy(file.path(results_dir, "core_genes.txt"),
            file.path(results_dir, "core_genes.csv"))

  # ── Step 18: Hub gene ranking ────────────────────────────────────────────────
  log_file("Ranking hub genes", main_log)
  hub_genes <- list()
  for (mod in unique(moduleColors)) {
    mod_genes <- colnames(exprData)[moduleColors == mod]
    if (length(mod_genes) >= 2) {
      connectivity   <- intramodularConnectivity(adjacency, moduleColors, scaleByMax = FALSE)
      mod_connectivity <- connectivity[mod_genes, "kWithin"]
      ranked         <- mod_genes[order(mod_connectivity, decreasing = TRUE)]
      top_percent    <- ceiling(length(ranked) * hub_top_percent)
      hub_genes[[mod]] <- ranked[1:min(top_percent, length(ranked))]
      log_file(paste("Identified", length(hub_genes[[mod]]), "hub genes for module", mod), main_log)
    } else {
      log_file(paste("Skipping module", mod, "(<2 genes)"), main_log)
      hub_genes[[mod]] <- character(0)
    }
  }
  sink(file.path(results_dir, "hub_genes_ranked.txt")); print(hub_genes); sink()
  file.copy(file.path(results_dir, "hub_genes_ranked.txt"),
            file.path(results_dir, "hub_genes_ranked.csv"))

  # ── Step 19: Hub–lysosome overlap ────────────────────────────────────────────
  log_file("Calculating hub genes overlap with lysosome pathway", main_log)
  lysosome_kegg_genes <- tryCatch({
    kegg_data  <- KEGGREST::keggGet("hsa04142")
    entrez_ids <- unlist(strsplit(kegg_data[[1]]$GENE, ";"))
    entrez_ids <- entrez_ids[grep("^[0-9]+$", entrez_ids)]
    bitr(entrez_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")$SYMBOL
  }, error = function(e) { log_file(paste("Error fetching hsa04142:", e$message), main_log); lysosome_genes })

  hub_lysosome_overlap <- lapply(hub_genes, function(x) intersect(x, lysosome_kegg_genes))
  sink(file.path(results_dir, "hub_lysosome_pathway_overlap.txt")); print(hub_lysosome_overlap); sink()
  file.copy(file.path(results_dir, "hub_lysosome_pathway_overlap.txt"),
            file.path(results_dir, "hub_lysosome_pathway_overlap.csv"))

  # ── Step 20: Cytoscape export ────────────────────────────────────────────────
  log_file("Exporting Cytoscape networks", main_log)
  for (mod in names(core_genes)) {
    if (length(core_genes[[mod]]) >= 2) {
      edge_file <- file.path(results_dir, paste0("Cytoscape_edges_", mod, ".txt"))
      mod_TOM   <- TOM[core_genes[[mod]], core_genes[[mod]]]
      if (nrow(mod_TOM) >= 2) {
        exportNetworkToCytoscape(
          mod_TOM,
          edgeFile  = edge_file,
          nodeFile  = file.path(results_dir, paste0("Cytoscape_nodes_", mod, ".txt")),
          weighted  = TRUE, threshold = tom_threshold,
          nodeNames = core_genes[[mod]]
        )
        file.copy(edge_file, sub("\\.txt$", ".csv", edge_file))
        node_file <- file.path(results_dir, paste0("Cytoscape_nodes_", mod, ".txt"))
        file.copy(node_file, sub("\\.txt$", ".csv", node_file))
        log_file(paste("Cytoscape network exported for module", mod), main_log)
      }
    }
  }

  # ── Step 21: Lysosome gene distribution ─────────────────────────────────────
  log_file("Calculating lysosome gene distribution across modules", main_log)
  module_genes         <- split(colnames(exprData), moduleColors)
  lysosome_module_dist <- lapply(module_genes, function(x) intersect(unique(x), lysosome_genes))
  sink(file.path(results_dir, "lysosome_module_dist.txt"))
  cat("Number of lysosome genes in each module:\n")
  print(lapply(lysosome_module_dist, length))
  cat("\nLysosome genes in each module:\n")
  for (mod in names(lysosome_module_dist)) {
    cat("\nModule:", mod, "\n")
    cat("Lysosome genes:", paste(lysosome_module_dist[[mod]], collapse = ", "), "\n")
  }
  sink()
  file.copy(file.path(results_dir, "lysosome_module_dist.txt"),
            file.path(results_dir, "lysosome_module_dist.csv"))
  log_file("Lysosome gene distribution saved", main_log)

  # ── Step 22: Lysosome in high-correlation modules ────────────────────────────
  log_file("Identifying lysosome genes in high-correlation modules", main_log)
  lysosome_in_high_cor <- lapply(module_genes[high_cor_modules],
                                 function(x) intersect(x, lysosome_genes))
  hub_synergy <- lapply(high_cor_modules,
                        function(x) intersect(lysosome_module_dist[[x]], hub_genes[[x]]))
  sink(file.path(results_dir, "lysosome_in_high_cor_modules.txt"))
  print(lysosome_in_high_cor); print(hub_synergy)
  sink()
  file.copy(file.path(results_dir, "lysosome_in_high_cor_modules.txt"),
            file.path(results_dir, "lysosome_in_high_cor_modules.csv"))

  
  # ── Priority pathways 定义（必须加在这里！否则后面 filtered_* 会崩溃） ─────
  priority_go_terms <- c("GO:0005764", "GO:0006914", "GO:0005765", "GO:0000038", 
                         "GO:0006909", "GO:0006629", "GO:0043202", "GO:0008233", 
                         "GO:0006954", "GO:0019221", "GO:0002376", "GO:0006955", 
                         "GO:0005773", "GO:1905671", "GO:0090117", "GO:0007040")
  
  priority_kegg_paths <- c("hsa04142", "hsa04140", "hsa04144", "hsa04145", 
                           "hsa05010", "hsa04146", "hsa04150", "hsa04210", 
                           "hsa04141", "hsa05012", "hsa04130")
  
  priority_ra_kegg <- c("hsa04064", "hsa04620", "hsa04630", "hsa05323", 
                        "hsa04668", "hsa04060", "hsa04621")
  
  priority_reactome_paths <- c("R-HSA-162588", "R-HSA-432720", "R-HSA-2213236", 
                               "R-HSA-190829", "R-HSA-4088210", "R-HSA-9036727")
  
  log_file("Priority pathways defined successfully", main_log)
  
  # ── Step 23: Significant module enrichment ───────────────────────────────────
  log_file("Performing enrichment analysis for significant modules", main_log)
  significant_modules <- sub("ME", "", rownames(ME_pheno_pval)[ME_pheno_pval[, 1] < enrich_padj_cutoff])
  results <- list()

  for (mod in significant_modules) {
    if (length(module_genes[[mod]]) >= 2) {
      log_file(paste("Processing module:", mod), main_log)
      sink(file.path(results_dir, paste0("selected_module_info_", mod, ".txt")))
      print(paste("Module:", mod))
      print(length(module_genes[[mod]]))
      print(head(module_genes[[mod]]))
      sink()
      file.copy(file.path(results_dir, paste0("selected_module_info_", mod, ".txt")),
                file.path(results_dir, paste0("selected_module_info_", mod, ".csv")))

      go_enrich <- tryCatch(
        enrichGO(gene = unique(module_genes[[mod]]), OrgDb = "org.Hs.eg.db",
                 keyType = "SYMBOL", ont = "ALL",
                 pAdjustMethod = "BH", qvalueCutoff = enrich_qvalue_cutoff),
        error = function(e) { log_file(paste("GO error module", mod, ":", e$message), main_log); NULL }
      )
      go_enrich_df <- if (!is.null(go_enrich)) as.data.frame(go_enrich) else data.frame()
      write.csv(go_enrich_df,
                file.path(results_dir, paste0("go_enrichment_", mod, ".csv")), row.names = FALSE)
      log_file(paste("GO enrichment module", mod, ":", nrow(go_enrich_df), "terms"), main_log)
      if (nrow(go_enrich_df) > 0) {
        # [MODIFIED] Triple-format output for GO dotplot
        save_base_fig(quote({
          print(dotplot(go_enrich, showCategory = 10) + ggtitle(paste("GO Enrichment in", mod, "Module")))
        }), file.path(results_dir, paste0("go_enrichment_dotplot_", mod)), w_in = 10, h_in = 8)
      }

      gene_df    <- bitr(unique(module_genes[[mod]]), fromType = "SYMBOL",
                         toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
      
      
      
      # ── Local KEGG enrichment for this module (no internet, no timeout) ─────────────
      kegg_file <- file.path(data_dir, "kegg_jp_link_hsa_pathway.txt")
      if (!file.exists(kegg_file)) stop("Local KEGG file not found: ", kegg_file)
      
      # Read and clean the local file
      df <- read.table(kegg_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
      colnames(df) <- c("Pathway", "Gene")
      df$Gene <- gsub("^hsa:", "", df$Gene)   # Remove "hsa:" prefix
      pathway_list <- split(df$Gene, df$Pathway)
      names(pathway_list) <- gsub("path:", "", names(pathway_list))
      
      # Build TERM2GENE
      term2gene <- data.frame(
        term = rep(names(pathway_list), lengths(pathway_list)),
        gene = unlist(pathway_list)
      )
      
      # Local KEGG enrichment using enricher
      kegg_enrich <- tryCatch({
        enricher(gene = gene_df$ENTREZID,
                 TERM2GENE = term2gene,
                 pvalueCutoff = enrich_qvalue_cutoff, qvalueCutoff = enrich_qvalue_cutoff)
      }, error = function(e) {
        log_file(paste("KEGG error module", mod, ":", e$message), main_log)
        NULL
      })
      
      kegg_enrich_df <- if (!is.null(kegg_enrich)) as.data.frame(kegg_enrich) else data.frame()
      
      write.csv(kegg_enrich_df,
                file.path(results_dir, paste0("kegg_enrichment_", mod, ".csv")), row.names = FALSE)
      log_file(paste("KEGG enrichment module", mod, ":", nrow(kegg_enrich_df), "pathways"), main_log)
      
 
      if (nrow(kegg_enrich_df) > 0) {
        # [MODIFIED] Triple-format output for KEGG dotplot
        save_base_fig(quote({
          print(dotplot(kegg_enrich, showCategory = 10) + ggtitle(paste("KEGG Enrichment in", mod, "Module")))
        }), file.path(results_dir, paste0("kegg_enrichment_dotplot_", mod)), w_in = 10, h_in = 8)
      }

      reactome_enrich <- tryCatch(
        enrichPathway(gene = gene_df$ENTREZID, organism = "human",
                      pAdjustMethod = "BH", qvalueCutoff = enrich_qvalue_cutoff),
        error = function(e) { log_file(paste("Reactome error module", mod, ":", e$message), main_log); NULL }
      )
      reactome_enrich_df <- if (!is.null(reactome_enrich)) as.data.frame(reactome_enrich) else data.frame()
      write.csv(reactome_enrich_df,
                file.path(results_dir, paste0("reactome_enrichment_", mod, ".csv")), row.names = FALSE)
      log_file(paste("Reactome enrichment module", mod, ":", nrow(reactome_enrich_df), "pathways"), main_log)
      if (nrow(reactome_enrich_df) > 0) {
        # [MODIFIED] Triple-format output for Reactome dotplot
        save_base_fig(quote({
          print(dotplot(reactome_enrich, showCategory = 10) + ggtitle(paste("Reactome Enrichment in", mod, "Module")))
        }), file.path(results_dir, paste0("reactome_enrichment_dotplot_", mod)), w_in = 10, h_in = 8)
      }

      filtered_go       <- subset(go_enrich_df,       ID %in% priority_go_terms)
      filtered_kegg     <- subset(kegg_enrich_df,     ID %in% c(priority_kegg_paths, priority_ra_kegg))
      filtered_reactome <- subset(reactome_enrich_df, ID %in% priority_reactome_paths)
      write.csv(filtered_go,       file.path(results_dir, paste0("filtered_go_",       mod, ".csv")), row.names = FALSE)
      write.csv(filtered_kegg,     file.path(results_dir, paste0("filtered_kegg_",     mod, ".csv")), row.names = FALSE)
      write.csv(filtered_reactome, file.path(results_dir, paste0("filtered_reactome_", mod, ".csv")), row.names = FALSE)

      results[[mod]] <- list(mod = mod,
                             go = go_enrich_df, kegg = kegg_enrich_df, reactome = reactome_enrich_df,
                             filtered_go = filtered_go, filtered_kegg = filtered_kegg,
                             filtered_reactome = filtered_reactome)
    }
  }

  log_file("Generating significant modules summary", main_log)
  sink(file.path(results_dir, "significant_modules_summary.txt"))
  cat("Significant modules (p < 0.05):\n"); print(significant_modules)
  cat("\nGene counts in significant modules:\n")
  for (mod in significant_modules) {
    mod_g <- unique(module_genes[[mod]])
    cat("\nModule:", mod, "\n")
    cat("Number of genes:", length(mod_g), "\n")
    cat("Genes:", paste(head(mod_g, 10), collapse = ", "),
        ifelse(length(mod_g) > 10, "...", ""), "\n")
  }
  sink()
  file.copy(file.path(results_dir, "significant_modules_summary.txt"),
            file.path(results_dir, "significant_modules_summary.csv"))
  log_file("Significant modules summary saved", main_log)

  # ── Step 24: Lysosome output directory ──────────────────────────────────────
  lysosome_output_dir <- file.path(results_dir, "Lysosome_Related_Files")
  if (!dir.exists(lysosome_output_dir)) dir.create(lysosome_output_dir)

  lysosome_hub_genes  <- lapply(hub_genes,  function(x) intersect(x, lysosome_genes))
  lysosome_core_genes <- lapply(core_genes, function(x) intersect(x, lysosome_genes))

  sink(file.path(lysosome_output_dir, "lysosome_hub_genes.txt"))
  cat("Lysosome hub genes per module:\n")
  for (mod in names(lysosome_hub_genes)) {
    cat("\nModule:", mod, "\n")
    if (length(lysosome_hub_genes[[mod]]) > 0)
      cat("Lysosome hub genes:", paste(lysosome_hub_genes[[mod]], collapse = ", "), "\n")
    else cat("No lysosome hub genes\n")
  }
  sink()
  file.copy(file.path(lysosome_output_dir, "lysosome_hub_genes.txt"),
            file.path(lysosome_output_dir, "lysosome_hub_genes.csv"))

  sink(file.path(lysosome_output_dir, "lysosome_core_genes.txt"))
  cat("Lysosome core genes per module (MM > 0.5, GS > 0.2):\n")
  for (mod in names(lysosome_core_genes)) {
    cat("\nModule:", mod, "\n")
    if (length(lysosome_core_genes[[mod]]) > 0)
      cat("Lysosome core genes:", paste(lysosome_core_genes[[mod]], collapse = ", "), "\n")
    else cat("No lysosome core genes\n")
  }
  sink()
  file.copy(file.path(lysosome_output_dir, "lysosome_core_genes.txt"),
            file.path(lysosome_output_dir, "lysosome_core_genes.csv"))
  log_file("Lysosome hub/core genes saved", main_log)

  lysosome_aux_dir  <- file.path(lysosome_output_dir, "Auxiliary_Files")
  if (!dir.exists(lysosome_aux_dir)) dir.create(lysosome_aux_dir)
  lysosome_modules  <- names(lysosome_hub_genes)[sapply(lysosome_hub_genes, length) > 0]
  file_types        <- c("go_enrichment_.*\\.csv","kegg_enrichment_.*\\.csv",
                         "reactome_enrichment_.*\\.csv","Cytoscape_edges_.*\\.txt","Cytoscape_nodes_.*\\.txt")
  for (mod in lysosome_modules) {
    for (file_type in file_types) {
      files <- list.files(results_dir, pattern = gsub("\\.\\*", mod, file_type), full.names = TRUE)
      for (file in files) {
        file.copy(file, file.path(lysosome_aux_dir, basename(file)), overwrite = TRUE)
        log_file(paste("Copied", basename(file), "to Auxiliary_Files"), main_log)
      }
    }
  }
  sig_sum_file <- file.path(results_dir, "significant_modules_summary.txt")
  if (file.exists(sig_sum_file))
    file.copy(sig_sum_file, file.path(lysosome_aux_dir, "significant_modules_summary.txt"), overwrite = TRUE)

  # ── Step 25: Lysosome genes summary CSV ─────────────────────────────────────
  log_file("Generating lysosome genes summary CSV", main_log)
  lysosome_summary_file <- file.path(lysosome_output_dir, "lysosome_genes_summary.csv")

  lysosome_summary <- data.frame(
    Symbol             = lysosome_genes,
    GeneID             = NA_character_,
    DE_Pvalue_RA_OA    = NA_real_,
    DE_Pvalue_RA_ND    = NA_real_,
    GS_Pvalue_RA_status = NA_real_,
    GS_Pvalue_RA_ND    = NA_real_,
    GS_Pvalue_RA_OA    = NA_real_,
    Hub_Gene           = "N",
    Core_Gene          = "N",
    Significant_Gene   = "N",
    Enrichment         = "",
    Key_Interactions   = "",
    Pathway_Gene       = "N",
    Pathway            = ""
  )

  lysosome_summary$GeneID           <- lysosome_df$GeneID[match(lysosome_summary$Symbol, lysosome_df$Symbol)]
  lysosome_summary$DE_Pvalue_RA_OA  <- top_table_RA_OA$adj.P.Val[match(lysosome_summary$Symbol, top_table_RA_OA$GeneSymbol)]
  lysosome_summary$DE_Pvalue_RA_ND  <- top_table_RA_ND$adj.P.Val[match(lysosome_summary$Symbol, top_table_RA_ND$GeneSymbol)]
  lysosome_summary$GS_Pvalue_RA_status <- GSPvalue[match(lysosome_summary$Symbol, rownames(GSPvalue)), 1]
  lysosome_summary$GS_Pvalue_RA_ND  <- GSPvalue_RA_ND[match(lysosome_summary$Symbol, rownames(GSPvalue_RA_ND)), 1]
  lysosome_summary$GS_Pvalue_RA_OA  <- GSPvalue_RA_OA[match(lysosome_summary$Symbol, rownames(GSPvalue_RA_OA)), 1]

  for (mod in names(hub_genes))  lysosome_summary$Hub_Gene[lysosome_summary$Symbol %in% hub_genes[[mod]]]  <- "Y"
  for (mod in names(core_genes)) lysosome_summary$Core_Gene[lysosome_summary$Symbol %in% core_genes[[mod]]] <- "Y"
  for (mod in significant_modules)
    lysosome_summary$Significant_Gene[lysosome_summary$Symbol %in% module_genes[[mod]]] <- "Y"

  # Populate Enrichment column from module enrichment files
  for (mod in significant_modules) {
    for (enr_type in c("go","kegg","reactome")) {
      enr_file <- file.path(lysosome_aux_dir, paste0(enr_type, "_enrichment_", mod, ".csv"))
      if (file.exists(enr_file)) {
        enr_df    <- read.csv(enr_file)
        enr_terms <- enr_df$ID[enr_df$p.adjust < enrich_padj_cutoff]
        mod_g     <- module_genes[[mod]]
        for (gene in intersect(lysosome_genes, mod_g)) {
          if (length(enr_terms) > 0) {
            entries <- paste0(toupper(enr_type), ":", enr_terms, " (", basename(enr_file), ")")
            lysosome_summary$Enrichment[lysosome_summary$Symbol == gene] <-
              paste(lysosome_summary$Enrichment[lysosome_summary$Symbol == gene],
                    paste(entries, collapse = ","), sep = ",")
          }
        }
      }
    }
  }
  # Core lysosome enrichment entries
  for (gene in core_lysosome_genes) {
    for (enr_df2 in list(go_enrich_core_df, kegg_enrich_core_df, reactome_enrich_core_df)) {
      if (nrow(enr_df2) > 0) {
        enr_terms2 <- enr_df2$ID[enr_df2$p.adjust < enrich_padj_cutoff]
        prefix     <- if ("category" %in% colnames(enr_df2)) "KEGG" else
                      if ("GeneRatio" %in% colnames(enr_df2) && !"category" %in% colnames(enr_df2)) "Reactome" else "GO"
        if (length(enr_terms2) > 0) {
          entries <- paste0(prefix, ":", enr_terms2, " (core_lysosome)")
          lysosome_summary$Enrichment[lysosome_summary$Symbol == gene] <-
            paste(lysosome_summary$Enrichment[lysosome_summary$Symbol == gene],
                  paste(entries, collapse = ","), sep = ",")
        }
      }
    }
  }

  # Populate Key_Interactions
  for (mod in significant_modules) {
    edge_file <- file.path(lysosome_aux_dir, paste0("Cytoscape_edges_", mod, ".txt"))
    if (file.exists(edge_file)) {
      edges <- read.table(edge_file, header = TRUE, sep = "\t")
      edges <- edges[edges$weight > tom_threshold, ]
      for (gene in intersect(lysosome_genes, module_genes[[mod]])) {
        interactions <- c(edges$fromNode[edges$toNode == gene],
                          edges$toNode[edges$fromNode == gene])
        if (length(interactions) > 0) {
          interaction_entries <- paste0(interactions, "-", gene, " (", basename(edge_file), ")")
          lysosome_summary$Key_Interactions[lysosome_summary$Symbol == gene] <-
            paste(lysosome_summary$Key_Interactions[lysosome_summary$Symbol == gene],
                  paste(interaction_entries, collapse = ","), sep = ",")
        }
      }
    }
  }

  # Populate Pathway columns
  for (mod in significant_modules) {
    kegg_file2 <- file.path(lysosome_aux_dir, paste0("kegg_enrichment_", mod, ".csv"))
    if (file.exists(kegg_file2)) {
      kegg_df2 <- read.csv(kegg_file2)
      for (pathway in priority_kegg_paths) {
        if (pathway %in% kegg_df2$ID) {
          genes2       <- unlist(strsplit(kegg_df2$geneID[kegg_df2$ID == pathway], "/"))
          gene_symbols2 <- bitr(genes2, fromType = "ENTREZID", toType = "SYMBOL",
                                OrgDb = "org.Hs.eg.db")$SYMBOL
          lysosome_summary$Pathway_Gene[lysosome_summary$Symbol %in% gene_symbols2] <- "Y"
          lysosome_summary$Pathway[lysosome_summary$Symbol %in% gene_symbols2] <-
            paste(lysosome_summary$Pathway[lysosome_summary$Symbol %in% gene_symbols2],
                  paste0("KEGG:", pathway), sep = ",")
        }
      }
    }
    reactome_file2 <- file.path(lysosome_aux_dir, paste0("reactome_enrichment_", mod, ".csv"))
    if (file.exists(reactome_file2)) {
      reactome_df2 <- read.csv(reactome_file2)
      for (pathway in priority_reactome_paths) {
        if (pathway %in% reactome_df2$ID) {
          genes3        <- unlist(strsplit(reactome_df2$geneID[reactome_df2$ID == pathway], "/"))
          gene_symbols3 <- bitr(genes3, fromType = "ENTREZID", toType = "SYMBOL",
                                OrgDb = "org.Hs.eg.db")$SYMBOL
          lysosome_summary$Pathway_Gene[lysosome_summary$Symbol %in% gene_symbols3] <- "Y"
          lysosome_summary$Pathway[lysosome_summary$Symbol %in% gene_symbols3] <-
            paste(lysosome_summary$Pathway[lysosome_summary$Symbol %in% gene_symbols3],
                  paste0("Reactome:", pathway), sep = ",")
        }
      }
    }
  }

  lysosome_summary$Enrichment       <- gsub("^,+|,+$", "", lysosome_summary$Enrichment)
  lysosome_summary$Key_Interactions <- gsub("^,+|,+$", "", lysosome_summary$Key_Interactions)
  lysosome_summary$Pathway          <- gsub("^,+|,+$", "", lysosome_summary$Pathway)
  write.csv(lysosome_summary, lysosome_summary_file, row.names = FALSE)
  log_file("Lysosome genes summary CSV saved", main_log)

  # Per-gene CSV files
  lysosome_genes_summary_dir <- file.path(lysosome_output_dir, "lysosome_genes_summary")
  if (!dir.exists(lysosome_genes_summary_dir)) dir.create(lysosome_genes_summary_dir)

  for (i in seq_len(nrow(lysosome_summary))) {
    gene_symbol   <- lysosome_summary$Symbol[i]
    csv_file_name <- paste0(gene_symbol, "_", dataset_id, "_.csv")
    csv_file_path <- file.path(lysosome_genes_summary_dir, csv_file_name)
    gene_data <- data.frame(
      Key   = c("Symbol","GeneID","DE_Pvalue_RA_OA","DE_Pvalue_RA_ND",
                "GS_Pvalue_RA_status","GS_Pvalue_RA_ND","GS_Pvalue_RA_OA",
                "Hub_Gene","Core_Gene","Significant_Gene",
                "Enrichment","Key_Interactions","Pathway_Gene","Pathway"),
      Value = c(
        gene_symbol,
        ifelse(is.na(lysosome_summary$GeneID[i]),             "NA", lysosome_summary$GeneID[i]),
        ifelse(is.na(lysosome_summary$DE_Pvalue_RA_OA[i]),    "NA", lysosome_summary$DE_Pvalue_RA_OA[i]),
        ifelse(is.na(lysosome_summary$DE_Pvalue_RA_ND[i]),    "NA", lysosome_summary$DE_Pvalue_RA_ND[i]),
        ifelse(is.na(lysosome_summary$GS_Pvalue_RA_status[i]),"NA", lysosome_summary$GS_Pvalue_RA_status[i]),
        ifelse(is.na(lysosome_summary$GS_Pvalue_RA_ND[i]),    "NA", lysosome_summary$GS_Pvalue_RA_ND[i]),
        ifelse(is.na(lysosome_summary$GS_Pvalue_RA_OA[i]),    "NA", lysosome_summary$GS_Pvalue_RA_OA[i]),
        lysosome_summary$Hub_Gene[i],
        lysosome_summary$Core_Gene[i],
        lysosome_summary$Significant_Gene[i],
        lysosome_summary$Enrichment[i],
        lysosome_summary$Key_Interactions[i],
        lysosome_summary$Pathway_Gene[i],
        lysosome_summary$Pathway[i]
      )
    )
    write.csv(gene_data, csv_file_path, row.names = FALSE)
    log_file(paste("Generated CSV for gene:", gene_symbol), main_log)
  }

  # ── Step 26: Target gene Genes/ directory ───────────────────────────────────
  genes_dir <- file.path(lysosome_output_dir, "Genes")
  if (!dir.exists(genes_dir)) dir.create(genes_dir)

  for (gene in target_genes) {
    gene_dir      <- file.path(genes_dir, gene)
    if (!dir.exists(gene_dir)) dir.create(gene_dir)
    csv_file_name <- paste0(gene, "_", dataset_id, "_.csv")
    src_csv       <- file.path(lysosome_genes_summary_dir, csv_file_name)
    if (file.exists(src_csv)) {
      file.copy(src_csv, file.path(gene_dir, csv_file_name), overwrite = TRUE)
      log_file(paste("Copied", csv_file_name, "to Genes/", gene), main_log)
    } else {
      log_file(paste("CSV for gene", gene, "not found in lysosome_genes_summary"), main_log)
    }
  }

  # Repair CSV helper
  repair_csv <- function(input_file, output_file) {
    lines <- readLines(input_file)
    repaired_lines <- c(lines[1])
    i <- 2
    enrichment_value <- NULL
    while (i <= length(lines)) {
      line   <- lines[i]
      fields <- strsplit(line, ",")[[1]]
      if (length(fields) >= 1 && fields[1] == "Enrichment") {
        enrichment_parts <- fields[2:length(fields)]
        i <- i + 1
        while (i <= length(lines)) {
          next_line   <- lines[i]
          next_fields <- strsplit(next_line, ",")[[1]]
          if (length(next_fields) >= 1 && next_fields[1] == "Key_Interactions") break
          enrichment_parts <- c(enrichment_parts, paste(next_fields, collapse = " "))
          i <- i + 1
        }
        enrichment_value <- gsub("\\s+", " ", trimws(paste(enrichment_parts, collapse = " ")))
        repaired_lines   <- c(repaired_lines, paste("Enrichment", enrichment_value, sep = ","))
      } else {
        repaired_lines <- c(repaired_lines, line)
        i <- i + 1
      }
    }
    writeLines(repaired_lines, output_file)
  }

  # ── Step 27: Pathway classification & interaction categories per target gene ─
  target_geneids <- c("80342","3002","3117","27074","4818","64116","5551","81552")
  target_symbols <- c("TRAF3IP3","GZMB","HLA-DQA1","LAMP3","NKG7","SLC39A8","PRF1","VOPP1")
  target_map     <- setNames(target_geneids, target_symbols)

  go_terms_class <- c("GO:0140507","GO:0002399","GO:0002503","GO:0002396","GO:0002501",
    "GO:0019886","GO:0002495","GO:0002504","GO:0002478","GO:0019884","GO:0048002","GO:0042267",
    "GO:0007159","GO:0002228","GO:0019882","GO:0031640","GO:0141061","GO:0141060","GO:0001909",
    "GO:0002420","GO:0033212","GO:0071421","GO:0050870","GO:0001906","GO:0002423","GO:0046689",
    "GO:0071578","GO:0006828","GO:0043320","GO:1903039","GO:0022409","GO:0051251","GO:0006525",
    "GO:0035455","GO:0002696","GO:0097501","GO:0002449","GO:0050867","GO:1903037","GO:0050863",
    "GO:0034755","GO:0061687","GO:0071276","GO:0050727","GO:0042832","GO:0071577","GO:0006829",
    "GO:0072594","GO:0001562","GO:0002418","GO:0002443","GO:0002323","GO:0045785","GO:0006882",
    "GO:0015701","GO:1903749","GO:0070269","GO:0043154","GO:0002347","GO:0046686","GO:1903747",
    "GO:2000117","GO:1901799","GO:0006826","GO:0061756")
  kegg_terms_class <- c("hsa05330","hsa04940","hsa05332","hsa05320","hsa05310","hsa04672",
    "hsa05321","hsa05416","hsa05140","hsa04612","hsa04658","hsa05323","hsa04640","hsa05150",
    "hsa04659","hsa05145","hsa05322","hsa04145","hsa04514","hsa05164","hsa05152","hsa05168",
    "hsa05169","hsa05166","hsa04216")
  reactome_terms_class <- c("R-HSA-202430","R-HSA-202427","R-HSA-389948","R-HSA-202433",
    "R-HSA-388841","R-HSA-877300","R-HSA-202424","R-HSA-202403","R-HSA-2132295","R-HSA-913531",
    "R-HSA-442380","R-HSA-2197563","R-HSA-435354","R-HSA-425410","R-HSA-5620971","R-HSA-1980145",
    "R-HSA-9725371","R-HSA-109606","R-HSA-5218859","R-HSA-425366","R-HSA-9700206","R-HSA-9725370")
  reactome_descriptions <- c(
    "TLR3 cascade","TLR4 cascade","Macroautophagy","TLR9 cascade",
    "MyD88-independent TLR cascade","Interferon gamma signaling","TLR7/8 cascade",
    "TLRs cascade","MHC class II antigen presentation","Interferon signaling",
    "PD-1 signaling","NOTCH1 signaling","Endosomal/vacuolar pathway","Autophagy",
    "RIG-I-like receptor signaling","Innate immune system",
    "Interferon alpha/beta signaling","Intrinsic pathway for apoptosis",
    "TGF-beta receptor signaling","Autophagy","Interferon signaling","Interferon signaling"
  )

  for (gene in target_symbols) {
    gene_dir  <- file.path(genes_dir, gene)
    if (!dir.exists(gene_dir)) dir.create(gene_dir)
    csv_file  <- file.path(gene_dir, paste0(gene, "_", dataset_id, "_.csv"))
    if (!file.exists(csv_file)) next

    # Repair CSV
    repaired_file <- sub("\\.csv$", "_repaired.csv", csv_file)
    repair_csv(csv_file, repaired_file)
    log_file(paste("Repaired CSV for gene", gene), main_log)

    # Pathway classification
    get_go_info <- function(go_id) {
      term     <- tryCatch(Term(go_id),     error = function(e) "Unknown")
      ontology <- tryCatch(Ontology(go_id), error = function(e) "Unknown")
      data.frame(ID = go_id, Term = term, Category = ontology)
    }
    go_info <- bind_rows(lapply(go_terms_class, get_go_info))

    get_kegg_term <- function(kegg_id) {
      tryCatch({ info <- keggGet(kegg_id); info[[1]]$NAME }, error = function(e) "Unknown")
    }
    kegg_info <- data.frame(
      ID       = kegg_terms_class,
      Term     = sapply(kegg_terms_class, get_kegg_term),
      Category = "KEGG Pathway"
    )
    reactome_info <- data.frame(
      ID       = reactome_terms_class,
      Term     = reactome_descriptions,
      Category = "Reactome Pathway"
    )
    all_info <- bind_rows(go_info, kegg_info, reactome_info) %>%
      mutate(Gene = gene, Context = "core_lysosome")

    write.csv(all_info,
              file.path(gene_dir, paste0(gene, "_pathway_classification.csv")),
              row.names = FALSE)
    log_file(paste("Generated pathway classification for", gene), main_log)

    p_pathway <- ggplot(all_info, aes(x = Category, fill = Category)) +
      geom_bar() + theme_minimal() +
      labs(title = paste("Distribution of Pathway Categories for", gene),
           x = "Category", y = "Count") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    # [MODIFIED] Triple-format output via save_fig
    save_fig(p_pathway, file.path(gene_dir, paste0(gene, "_pathway_distribution")), w_in = 8, h_in = 6)

    # Interaction categories
    data_csv <- read.csv(csv_file, stringsAsFactors = FALSE)
    interaction_row <- data_csv[data_csv$Key == "Key_Interactions", ]
    if (nrow(interaction_row) == 0 || is.na(interaction_row$Value) ||
        interaction_row$Value == "") {
      log_file(paste("Key_Interactions empty for gene", gene), main_log); next
    }
    interaction_items <- unlist(strsplit(interaction_row$Value, ","))
    partner_genes <- sapply(interaction_items, function(x) {
      trimws(gsub(paste0("-", gene, "\\s*\\(.*\\)"), "", x))
    })
    gene_categories <- data.frame(gene = partner_genes, category = "other")
    gene_categories$category[gene_categories$gene %in% c("CD3D","CD3E","CD4")] <- "t_cell"
    gene_categories$category[gene_categories$gene %in% c("PAX5","CD19")]        <- "b_cell"
    gene_categories$category[gene_categories$gene %in% c("IFNG","IL2")]         <- "cytokine"

    category_counts <- gene_categories %>%
      group_by(category) %>%
      summarise(count = n(), genes = paste(gene, collapse = "/")) %>%
      ungroup()

    all_categories  <- c("t_cell","b_cell","cytokine","transcription_factor",
                         "cytotoxicity","signaling","other")
    category_counts <- data.frame(category = all_categories) %>%
      left_join(category_counts, by = "category") %>%
      mutate(count = ifelse(is.na(count), 0, count),
             genes = ifelse(is.na(genes), "", genes))

    write.csv(category_counts,
              file.path(gene_dir, paste0(gene, "_interaction_categories.csv")),
              row.names = FALSE)
    log_file(paste("Generated interaction categories for", gene), main_log)

    p_interact <- ggplot(category_counts, aes(x = category, y = count, fill = category)) +
      geom_bar(stat = "identity") + theme_minimal() +
      labs(title = paste("Distribution of Interaction Categories for", gene),
           x = "Category", y = "Count") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    # [MODIFIED] Triple-format output via save_fig
    save_fig(p_interact, file.path(gene_dir, paste0(gene, "_interaction_distribution")), w_in = 8, h_in = 6)
  }

  # P-value validation
  log_file("Validating P_value consistency in lysosome_genes_summary.csv", main_log)
  summary_data <- read.csv(lysosome_summary_file, stringsAsFactors = FALSE)
  merged_de_RA_OA <- merge(summary_data[, c("Symbol","DE_Pvalue_RA_OA")],
                           top_table_RA_OA[, c("GeneSymbol","adj.P.Val")],
                           by.x = "Symbol", by.y = "GeneSymbol", all.x = TRUE)
  consistent_de_RA_OA <- abs(merged_de_RA_OA$DE_Pvalue_RA_OA - merged_de_RA_OA$adj.P.Val) < 1e-10 &
    !is.na(merged_de_RA_OA$DE_Pvalue_RA_OA) & !is.na(merged_de_RA_OA$adj.P.Val)
  log_file(paste("Consistent DE_Pvalues (RA vs OA):", sum(consistent_de_RA_OA),
                 "out of", nrow(merged_de_RA_OA)), main_log)

  merged_de_RA_ND <- merge(summary_data[, c("Symbol","DE_Pvalue_RA_ND")],
                           top_table_RA_ND[, c("GeneSymbol","adj.P.Val")],
                           by.x = "Symbol", by.y = "GeneSymbol", all.x = TRUE)
  consistent_de_RA_ND <- abs(merged_de_RA_ND$DE_Pvalue_RA_ND - merged_de_RA_ND$adj.P.Val) < 1e-10 &
    !is.na(merged_de_RA_ND$DE_Pvalue_RA_ND) & !is.na(merged_de_RA_ND$adj.P.Val)
  log_file(paste("Consistent DE_Pvalues (RA vs ND):", sum(consistent_de_RA_ND),
                 "out of", nrow(merged_de_RA_ND)), main_log)

  write.csv(merged_de_RA_OA,
            file.path(lysosome_output_dir, "validation_de_pvalues_RA_vs_oa.csv"), row.names = FALSE)
  write.csv(merged_de_RA_ND,
            file.path(lysosome_output_dir, "validation_de_pvalues_RA_vs_nd.csv"), row.names = FALSE)
  log_file("P-value validation results saved", main_log)

  # ── Step 28: Lysosome heatmap ────────────────────────────────────────────────
  log_file("Generating heatmap for lysosome genes in significant modules", main_log)
  significant_lysosome_genes <- unlist(lapply(lysosome_core_genes[significant_modules], unique))
  if (length(significant_lysosome_genes) > 0) {
    lysosome_expr_subset <- exprData[, significant_lysosome_genes, drop = FALSE]
    # [MODIFIED] Triple-format output for lysosome heatmap
    save_base_fig(quote({
      pheatmap(lysosome_expr_subset, scale = "column",
               annotation_row = phenoData[, "RA_status", drop = FALSE],
               main = "Lysosome Genes Heatmap in Significant Modules",
               show_rownames = TRUE, show_colnames = TRUE)
    }), file.path(lysosome_output_dir, "lysosome_genes_heatmap"), w_in = 10, h_in = 8)
    log_file("Lysosome genes heatmap generated (.pdf/_300.png/_600.tiff)", main_log)
  } else {
    log_file("No significant lysosome genes, skipping heatmap", main_log)
  }

  # ── Step 29: ROC validation ──────────────────────────────────────────────────
  log_file("Performing ROC validation for hub and core genes", main_log)

  hub_gene_vector  <- unlist(hub_genes)
  valid_hub_genes  <- hub_gene_vector[hub_gene_vector %in% colnames(exprData)]
  if (length(valid_hub_genes) == 0) stop("No valid hub genes for ROC analysis")
  log_file(paste("Found", length(valid_hub_genes), "valid hub genes for ROC"), main_log)
  hub_expr       <- exprData[, valid_hub_genes, drop = FALSE]
  roc_list_hub   <- apply(hub_expr, 2, function(x) roc(pheno_numeric, x))
  auc_values_hub <- sapply(roc_list_hub, auc)
  write.csv(data.frame(Gene = names(auc_values_hub), AUC = auc_values_hub),
            file.path(results_dir, "hub_roc_auc.csv"), row.names = FALSE)
  log_file("ROC analysis for hub genes completed", main_log)

  core_gene_vector  <- unlist(core_genes)
  valid_core_genes  <- core_gene_vector[core_gene_vector %in% colnames(exprData)]
  if (length(valid_core_genes) > 0) {
    core_expr       <- exprData[, valid_core_genes, drop = FALSE]
    roc_list_core   <- apply(core_expr, 2, function(x) roc(pheno_numeric, x))
    auc_values_core <- sapply(roc_list_core, auc)
    write.csv(data.frame(Gene = names(auc_values_core), AUC = auc_values_core),
              file.path(results_dir, "core_roc_auc.csv"), row.names = FALSE)
    log_file("ROC analysis for core genes completed", main_log)
  }

  
  
  # ── Step 30: GSEA ────────────────────────────────────────────────────────────
  
  # ── Step 30: GSEA ────────────────────────────────────────────────────────────
  # 【完整替换版】—— 已全部优化：SerialParam（防Windows死锁）+ nPermSimple=1000 + maxGSSize=300
  # 直接复制整个代码块替换原 Step 30 即可（从 library(BiocParallel) 开始）
  
  library(BiocParallel)
  register(SerialParam())   # ← 关键！改成串行，彻底解决并行卡死问题
  
  log_file("Performing GSEA for all lysosome pathways", main_log)
  gsea_results <- list()
  
  checkpoint_gsea <- file.path(lysosome_output_dir, "checkpoint_gsea.rds")
  if (file.exists(checkpoint_gsea)) {
    gsea_results <- load_checkpoint(checkpoint_gsea, main_log, "GSEA results")
    log_file("✓ GSEA results loaded from checkpoint (instant skip)", main_log)
  } else {
    # 1. Prepare the gene list
    gene_list_entrez <- tryCatch({
      suppressWarnings({
        symbol_to_entrez <- bitr(rownames(geneTraitSignificance),
                                 fromType = "SYMBOL", toType = "ENTREZID",
                                 OrgDb = "org.Hs.eg.db")
      })
      log_file(paste("bitr completed:", nrow(symbol_to_entrez), "Entrez IDs"), main_log)
      gene_list_named <- geneTraitSignificance[, 1][
        match(symbol_to_entrez$SYMBOL, rownames(geneTraitSignificance))]
      names(gene_list_named) <- symbol_to_entrez$ENTREZID
      sort(gene_list_named[!is.na(names(gene_list_named))], decreasing = TRUE)
    }, error = function(e) {
      log_file(paste("bitr error:", e$message), main_log); NULL
    })
    
    if (!is.null(gene_list_entrez)) {
      
      # 2. GO GSEA（ont="ALL" 完整版，已优化速度）
      log_file("Starting GO GSEA (ALL + optimized params)...", main_log)
      gsea_go <- tryCatch({
        gseGO(geneList = gene_list_entrez,
              OrgDb = org.Hs.eg.db,
              keyType = "ENTREZID",
              ont = "ALL",                    # 保留 ALL（BP+CC+MF）
              minGSSize = gsea_min_size,
              maxGSSize = 300,                # ← 关键：从500降到300，速度翻倍
              pvalueCutoff = gsea_pvalue_cutoff,
              pAdjustMethod = "BH",
              by = "fgsea",
              nPermSimple = 1000,             # ← 关键：默认太高，改成1000
              seed = 123,
              verbose = TRUE)
      }, error = function(e) {
        log_file(paste("GO GSEA error:", e$message), main_log); NULL
      })
      
      if (!is.null(gsea_go)) {
        df_go <- as.data.frame(gsea_go)
        write.csv(df_go, file.path(lysosome_output_dir, "gsea_GO.csv"), row.names = FALSE)
        gsea_results[["GO"]] <- df_go
        log_file(paste("GO GSEA completed:", nrow(df_go), "terms"), main_log)
      }
      
      # 3. Reactome GSEA（优化版）
      log_file("Starting Reactome GSEA...", main_log)
      gsea_reactome <- tryCatch({
        gsePathway(geneList = gene_list_entrez,
                   organism = "human",
                   minGSSize = gsea_min_size,
                   maxGSSize = 300,
                   pvalueCutoff = gsea_pvalue_cutoff,
                   pAdjustMethod = "BH",
                   by = "fgsea",
                   nPermSimple = 1000,      # ← 优化参数
                   seed = 123)
      }, error = function(e) {
        log_file(paste("Reactome GSEA error:", e$message), main_log); NULL
      })
      
      if (!is.null(gsea_reactome)) {
        df_reactome <- as.data.frame(gsea_reactome)
        write.csv(df_reactome, file.path(lysosome_output_dir, "gsea_Reactome.csv"), row.names = FALSE)
        gsea_results[["Reactome"]] <- df_reactome
        log_file(paste("Reactome GSEA completed:", nrow(df_reactome), "pathways"), main_log)
      }
      
      # 4. Local KEGG fgsea（优化版，使用您本地的 kegg_jp_link_hsa_pathway.txt）
      log_file("Starting local KEGG fgsea...", main_log)
      kegg_file <- file.path(data_dir, "kegg_jp_link_hsa_pathway.txt")
      if (file.exists(kegg_file)) {
        df <- read.table(kegg_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
        colnames(df) <- c("Pathway", "Gene")
        df$Gene <- gsub("^hsa:", "", df$Gene)
        pathway_list <- split(df$Gene, df$Pathway)
        names(pathway_list) <- gsub("path:", "", names(pathway_list))
        
        # 只保留您定义的优先路径（更快）
        selected_pathways <- pathway_list[names(pathway_list) %in% priority_kegg_paths]
        if (length(selected_pathways) == 0) selected_pathways <- pathway_list
        
        gsea_kegg <- tryCatch({
          fgsea(pathways = selected_pathways,
                stats = gene_list_entrez,
                minSize = gsea_min_size,
                maxSize = 300,
                nPermSimple = 1000,           # ← 优化参数
                eps = 0,
                seed = 123)
        }, error = function(e) {
          log_file(paste("KEGG fgsea error:", e$message), main_log); NULL
        })
        
        if (!is.null(gsea_kegg)) {
          df_kegg <- as.data.frame(gsea_kegg)
          write.csv(df_kegg, file.path(lysosome_output_dir, "gsea_KEGG_local.csv"), row.names = FALSE)
          gsea_results[["KEGG"]] <- df_kegg
          log_file(paste("Local KEGG fgsea completed:", nrow(df_kegg), "pathways"), main_log)
        }
      }
    }
    
    # 保存 checkpoint（以后再次运行直接秒过）
    save_checkpoint(gsea_results, checkpoint_gsea, main_log, "GSEA results")
    log_file("✓ GSEA checkpoint saved (future runs instant)", main_log)
  }
  
  log_file("Step 30 completed successfully!", main_log)
  
 
  # ── Step 31: Filter GSEA results per target gene ─────────────────────────────
  for (src_name in c("Reactome", "GO")) {
    gsea_src_file <- file.path(lysosome_output_dir, paste0("gsea_", src_name, ".csv"))
    if (!file.exists(gsea_src_file)) next
    gsea_df <- read.csv(gsea_src_file)
    col_name <- "core_enrichment"
    if (!col_name %in% colnames(gsea_df)) next
    for (gene in target_symbols) {
      gene_id <- target_map[gene]
      gene_dir <- file.path(genes_dir, gene)
      if (!dir.exists(gene_dir)) dir.create(gene_dir)
      filtered_df <- gsea_df[sapply(strsplit(gsea_df[[col_name]], "/"),
                                    function(x) gene_id %in% x), ]
      if (nrow(filtered_df) > 0) {
        dest_file <- file.path(gene_dir, paste0("gsea_", src_name, ".csv"))
        write.csv(filtered_df, dest_file, row.names = FALSE)
        log_file(paste("Filtered gsea_", src_name, "for gene", gene), main_log)
      }
    }
  }
  
  # KEGG 过滤
  gsea_kegg_files <- list.files(lysosome_output_dir, pattern = "^gsea_KEGG_.*\\.csv$", full.names = TRUE)
  for (kegg_file in gsea_kegg_files) {
    gsea_df <- read.csv(kegg_file)
    if (!"leadingEdge" %in% colnames(gsea_df)) next
    for (gene in target_symbols) {
      gene_id <- target_map[gene]
      gene_dir <- file.path(genes_dir, gene)
      if (!dir.exists(gene_dir)) dir.create(gene_dir)
      filtered_df <- gsea_df[sapply(strsplit(gsea_df$leadingEdge, "/"),
                                    function(x) gene_id %in% x), ]
      if (nrow(filtered_df) > 0) {
        write.csv(filtered_df, file.path(gene_dir, basename(kegg_file)), row.names = FALSE)
        log_file(paste("Filtered", basename(kegg_file), "for gene", gene), main_log)
      }
    }
  }
  
  # ── Step 32: Copy module enrichment files per target gene ────────────────────
  gene_to_module <- setNames(moduleColors, colnames(exprData))
  for (gene in target_genes) {
    if (gene %in% names(gene_to_module)) {
      mod <- gene_to_module[gene]
      gene_dir <- file.path(genes_dir, gene)
      if (!dir.exists(gene_dir)) dir.create(gene_dir)
      for (pattern in c("filtered_go_","filtered_kegg_","filtered_reactome_",
                        "go_enrichment_","kegg_enrichment_","reactome_enrichment_")) {
        src_file <- file.path(results_dir, paste0(pattern, mod, ".csv"))
        if (file.exists(src_file))
          file.copy(src_file, file.path(gene_dir, basename(src_file)), overwrite = TRUE)
      }
      edges_src <- file.path(results_dir, paste0("Cytoscape_edges_", mod, ".txt"))
      if (file.exists(edges_src)) {
        edges_df <- read.table(edges_src, header = TRUE, sep = "\t", quote = "", fill = TRUE)
        filtered_edges <- edges_df[edges_df$fromNode == gene | edges_df$toNode == gene, ]
        write.table(filtered_edges,
                    file.path(gene_dir, paste0("Cytoscape_edges_", mod, "_", gene, ".txt")),
                    sep = "\t", row.names = FALSE, quote = FALSE)
      }
      nodes_src <- file.path(results_dir, paste0("Cytoscape_nodes_", mod, ".txt"))
      if (file.exists(nodes_src)) {
        nodes_df <- read.table(nodes_src, header = TRUE, sep = "\t", quote = "", fill = TRUE)
        filtered_nodes <- nodes_df[nodes_df$nodeName == gene, ]
        write.table(filtered_nodes,
                    file.path(gene_dir, paste0("Cytoscape_nodes_", mod, "_", gene, ".txt")),
                    sep = "\t", row.names = FALSE, quote = FALSE)
      }
    } else {
      log_file(paste("Gene", gene, "not found in module assignment"), main_log)
    }
  }
  
  # ── Step 33: Clean up and final save ─────────────────────────────────────────

  if (dir.exists(temp_dir)) {
    unlink(temp_dir, recursive = TRUE)
    log_file("Temporary directory cleaned up", main_log)
  }
  
  saveRDS(exprData, file.path(results_dir, "checkpoint_raw_data.rds"))
  log_file("Saved checkpoint_raw_data.rds for module preservation", main_log)
  
  log_file(paste("WGCNA analysis completed successfully for", dataset_id), main_log)
}
   
# ==========================================================
# Cross-dataset module preservation analysis
# ==========================================================
cat("\n=== Starting Module Preservation Analysis ===\n")
preservation_dir <- file.path(base_dir, "ModulePreservation")
dir.create(preservation_dir, recursive = TRUE, showWarnings = FALSE)
load_dataset_for_preservation <- function(did, base_dir) {
  res_dir <- file.path(base_dir, did, "Results_WGCNA_GO-KEGG-Reactome")
  expr_file <- file.path(res_dir, "checkpoint_raw_data.rds")
  col_file <- file.path(res_dir, "moduleColors.csv")
  if (!file.exists(expr_file) || !file.exists(col_file)) {
    warning("Missing files for ", did); return(NULL)
  }
  expr <- readRDS(expr_file)
  if (nrow(expr) > ncol(expr)) expr <- t(expr)
  colors <- read.csv(col_file)
  mc <- setNames(colors$ModuleColor, colors$Gene)
  list(expr = expr, moduleColors = mc, dataset = did)
}
dataset_ids_pres <- c("GSE55235","GSE55457","GSE55584")
loaded_pres <- lapply(dataset_ids_pres, load_dataset_for_preservation, base_dir = base_dir)
names(loaded_pres) <- dataset_ids_pres
loaded_pres <- Filter(Negate(is.null), loaded_pres)
if (length(loaded_pres) >= 2) {
  multiData_pres <- lapply(loaded_pres, function(x) list(data = x$expr))
  multiColor_pres <- lapply(loaded_pres, function(x) x$moduleColors)
  set.seed(123)
  cat("Running modulePreservation (nPermutations=200)...\n")
  mp_result <- modulePreservation(
    multiData = multiData_pres,
    multiColor = multiColor_pres,
    referenceNetworks = 1,
    nPermutations = 200,
    randomSeed = 123,
    quickCor = 0,
    verbose = 3
  )
  saveRDS(mp_result, file.path(preservation_dir, "mp_result.rds"))
  Zsummary_all <- data.frame()
  for (test_idx in 2:length(loaded_pres)) {
    test_name <- names(loaded_pres)[test_idx]
    Z_df <- mp_result$preservation$Z[[1]][[test_idx]]
    if (is.null(Z_df) || nrow(Z_df) == 0) next
    z_col <- grep("Zsummary.pres|Zsummary", colnames(Z_df), value = TRUE)[1]
    d_col <- grep("Zsummary.qual|Zdensity", colnames(Z_df), value = TRUE)[1]
    Z_report <- data.frame(
      Module = rownames(Z_df),
      Zsummary = if (!is.na(z_col)) Z_df[, z_col] else NA,
      Zdensity = if (!is.na(d_col)) Z_df[, d_col] else NA,
      Test_Dataset = test_name,
      Ref_Dataset = names(loaded_pres)[1],
      stringsAsFactors = FALSE
    )
    Z_report$Status <- ifelse(Z_report$Zsummary > 10, "Strong",
                              ifelse(Z_report$Zsummary > 2, "Moderate", "Weak"))
    Zsummary_all <- rbind(Zsummary_all, Z_report)
  }
  target_auto <- tryCatch({
    unique(unlist(lapply(dataset_ids_pres, function(did) {
      f <- file.path(base_dir, did, "Results_WGCNA_GO-KEGG-Reactome", "high_cor_modules.csv")
      if (file.exists(f)) read.csv(f)$HighCorModules else character(0)
    })))
  }, error = function(e) character(0))
  Zsummary_target <- if (length(target_auto) > 0)
    Zsummary_all[Zsummary_all$Module %in% target_auto, ]
  else Zsummary_all
  write.csv(Zsummary_all, file.path(preservation_dir, "Zsummary_all_modules.csv"), row.names = FALSE)
  write.csv(Zsummary_target, file.path(preservation_dir, "Zsummary_target_modules.csv"), row.names = FALSE)
  cat("\n=== Module Preservation Zsummary (Lysosomal Modules) ===\n")
  print(Zsummary_target[, c("Module","Test_Dataset","Zsummary","Status")])
  if (nrow(Zsummary_target) > 0) {
    Zsummary_target$Status <- factor(Zsummary_target$Status, levels = c("Strong","Moderate","Weak"))
    p_pres <- ggplot(Zsummary_target, aes(x = Module, y = Zsummary, fill = Status)) +
      geom_col(position = "dodge", width = 0.7) +
      geom_hline(yintercept = 10, linetype = "dashed", color = "darkgreen", linewidth = 0.8) +
      geom_hline(yintercept = 2, linetype = "dashed", color = "orange", linewidth = 0.8) +
      facet_wrap(~ Test_Dataset, ncol = 1) +
      scale_fill_manual(values = c(Strong = "#27AE60", Moderate = "#F39C12", Weak = "#E74C3C")) +
      labs(title = "WGCNA Module Preservation (Zsummary)",
           subtitle = paste0("Reference: ", names(loaded_pres)[1]),
           x = "Module Color", y = "Zsummary", fill = "Preservation") +
      theme_classic(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    # [MODIFIED] Triple-format output via save_fig
    save_fig(p_pres, file.path(preservation_dir, "FigS2B_module_preservation"), w_in = 10, h_in = 8)
    cat("FigS2B saved to:", preservation_dir, "(.pdf/_300.png/_600.tiff)\n")
  }
} else {
  cat("Insufficient datasets loaded for preservation analysis.\n")
}

# Scale-free R² summary across datasets
sft_all <- do.call(rbind, lapply(dataset_ids_pres, function(did) {
  f <- file.path(base_dir, did, "Results_WGCNA_GO-KEGG-Reactome", "sft_selected_power_summary.csv")
  if (file.exists(f)) read.csv(f) else NULL
}))
if (!is.null(sft_all) && nrow(sft_all) > 0) {
  write.csv(sft_all, file.path(base_dir, "sft_all_datasets_summary.csv"), row.names = FALSE)
  cat("\nScale-free topology summary:\n")
  for (i in seq_len(nrow(sft_all))) {
    cat(sprintf(" %s: power=%d, R2=%.2f, mean.k=%.1f\n",
                sft_all$Dataset[i], sft_all$Soft_Power[i],
                sft_all$Scale_Free_R2[i], sft_all$Mean_Connectivity[i]))
  }
}
cat("\n=== 301_WGCNA all fixes applied ===\n")