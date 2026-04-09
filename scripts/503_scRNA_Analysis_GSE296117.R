# Load libraries (with installation check)
# if (!require("remotes", quietly = TRUE)) install.packages("remotes")
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("GEOquery", "SingleR", "celldex", "clusterProfiler", "org.Hs.eg.db", "monocle3"), ask = FALSE)
# if (!require("hdf5r", quietly = TRUE)) install.packages("hdf5r")
# if (!require("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")
# if (!require("STRINGdb", quietly = TRUE)) install.packages("STRINGdb")

if (!require("Seurat", quietly = TRUE)) {   # Install latest Seurat (v5), remove fixed 4.4.0 to support Assay5.4.0
  install.packages("SeuratObject")
  install.packages("Seurat")
}

library(GEOquery)
library(utils)
library(Seurat) # Must use Seurat 5
library(ggplot2)
library(dplyr)
library(patchwork)
library(hdf5r)
library(biomaRt)
library(pheatmap)
library(SingleR)
library(celldex)
library(clusterProfiler)
library(STRINGdb)
library(monocle3)
library(org.Hs.eg.db)
library(igraph) # Explicitly load to ensure upgrade_graph is available
Sys.setenv(UV_OFFLINE = "1") # Avoid uv online download

# ========================================================
# Global plot parameters (modify here to affect all plots)
# ========================================================
PLOT_DPI       <- 300   # Resolution for all ggsave outputs
PLOT_FONTSIZE  <- 10    # Base font size for ggplot2 theme
AXIS_TEXT_SIZE <- 10    # Axis tick label size
AXIS_TEXT_X_ANGLE <- 45 # X-axis label angle for DotPlot / VlnPlot

# ========================================================
# [STANDARDISED] Triple-format figure export helper
# Outputs: .pdf (vector), _300.png (submission), _600.tiff (production)
# ========================================================
save_fig <- function(p, stem, w_in, h_in) {
  ggsave(paste0(stem, ".pdf"), plot = p, width = w_in, height = h_in,
         units = "in", device = "pdf", useDingbats = FALSE)
  ggsave(paste0(stem, "_300.png"), plot = p, width = w_in, height = h_in,
         units = "in", device = "png", dpi = 300)
  ggsave(paste0(stem, "_600.tiff"), plot = p, width = w_in, height = h_in,
         units = "in", device = "tiff", dpi = 600, compression = "lzw")
}

# ==========================================================
VIOLIN_X_TEXT_SIZE    <- 9
VIOLIN_Y_TEXT_SIZE    <- 10   
PSEUDOTIME_Y_TEXT_SIZE <- 10  
# ===========================================================

PHEATMAP_FONTSIZE     <- 10  # pheatmap base font size
PHEATMAP_FONTSIZE_ROW <-  9  # pheatmap row label size
PHEATMAP_FONTSIZE_COL <-  9  # pheatmap column label size

# ========================================================
# Standard Path Management Template (All Modules)
# ========================================================
# install.packages("rprojroot") # Only execute once
library(rprojroot)
project_root <- find_rstudio_root_file()
data_dir     <- file.path(project_root, "data")
outputs_dir  <- file.path(project_root, "outputs")
# ------------------- Only modify this line -------------------
module_name <- "05_scRNA_Analysis\\GSE296117" # ←←← Modify here!!
# ------------------------------------------------------
module_dir <- file.path(outputs_dir, module_name)
dir.create(module_dir, showWarnings = FALSE, recursive = TRUE)

# Define file path
rds_file <- file.path(data_dir, "GSE296117_RA_geo.rds")
if (file.exists(rds_file)) {
  cat("Local RDS file already exists:", rds_file, "\n")
} else {
  stop("GSE296117 RDS file not found, please ensure the file is located in the data directory")
}

# Configuration section: Quality control parameters
qc_params <- list(
  min_nFeature_RNA = 10,
  max_nFeature_RNA = 6000,
  max_percent_mt   = 25,
  min_nCount_RNA   = 10,
  min_cells        = 5
)
cat("Quality control parameters:\n")
cat(" min_nFeature_RNA:", qc_params$min_nFeature_RNA, "\n")
cat(" max_nFeature_RNA:", qc_params$max_nFeature_RNA, "\n")
cat(" max_percent_mt:",  qc_params$max_percent_mt,    "\n")
cat(" min_nCount_RNA:",  qc_params$min_nCount_RNA,    "\n")
cat(" min_cells:",       qc_params$min_cells,         "\n")

# Define gene file list
gene_files <- list(
  list(file = file.path(data_dir, "Target Genes GZMB.txt"),    ref = NULL),
  list(file = file.path(data_dir, "Target Genes LAMP3.txt"),   ref = NULL),
  list(file = file.path(data_dir, "Target Genes NKG7.txt"),    ref = NULL),
  list(file = file.path(data_dir, "Target Genes SLC39A8.txt"), ref = NULL)
)
all_genes <- unique(unlist(lapply(gene_files, function(gf) gsub("_", "-", readLines(gf$file)))))

# Define dataset configuration list
datasets <- list(
  list(
    name = "GSE296117",
    base_output_dir = module_dir,
    load_func = function() {
      data <- readRDS(rds_file)
      if (inherits(data, "Seurat")) {
        obj <- data
      } else if (inherits(data, "SingleCellExperiment")) {
        obj <- as.Seurat(data)
      } else if (inherits(data, "dgCMatrix") || inherits(data, "matrix")) {
        rownames(data) <- make.unique(gsub("_", "-", rownames(data)))
        colnames(data) <- gsub("_", "-", colnames(data))
        obj <- CreateSeuratObject(counts = data, project = "GSE296117")
      } else {
        stop("Unknown object type, cannot convert to Seurat object")
      }
      tryCatch({
        rownames(obj[["RNA"]]) <- make.unique(gsub("_", "-", rownames(obj[["RNA"]])))
        colnames(obj) <- gsub("_", "-", colnames(obj))
      }, error = function(e) cat("Failed to set rownames/colnames:", e$message, "\n"))
      obj@meta.data$sample <- "GSE296117"
      if (inherits(obj[["RNA"]], "Assay")) {
        obj[["RNA"]] <- as(object = obj[["RNA"]], Class = "Assay5")
        cat("Converted Assay to Assay5 format\n")
      }
      if (inherits(obj[["RNA"]], "Assay5")) {
        if (!"data" %in% Layers(obj[["RNA"]])) {
          obj <- SetAssayData(obj, assay = "RNA", layer = "data",
                              new.data = GetAssayData(obj, assay = "RNA", layer = "counts"))
        }
      }
      return(obj)
    }
  )
)

# ========================================================
# Preprocessing function: execute common steps only once
# ========================================================
preprocess_dataset <- function(obj, dataset_name, base_output_dir, qc_params) {
  preprocessed_file <- file.path(base_output_dir, "preprocessed_seurat_obj.rds")
  if (file.exists(preprocessed_file)) {
    cat("Preprocessed Seurat object already exists, loading...\n")
    obj <- readRDS(preprocessed_file)
    return(obj)
  }
  
  cat("Starting preprocessing dataset:", dataset_name, "\n")
  DefaultAssay(obj) <- "RNA"
  cat("Current default assay:", DefaultAssay(obj), "\n")
  
  counts_data <- GetAssayData(obj, assay = "RNA", layer = "counts")
  data_data   <- GetAssayData(obj, assay = "RNA", layer = "data")
  cat("Initial RNA assay counts dimension:", dim(counts_data), "\n")
  cat("Initial RNA assay data dimension:",   dim(data_data),   "\n")
  cat("nFeature_RNA distribution:", summary(obj$nFeature_RNA), "\n")
  cat("nCount_RNA distribution:",   summary(obj$nCount_RNA),   "\n")
  
  # Ensembl ID → gene symbol conversion
  original_rownames <- rownames(counts_data)
  if (any(grepl("^ENSG", original_rownames))) {
    cat("Detected Ensembl IDs, converting to gene symbols...\n")
    tryCatch({
      ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      gene_map <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                        filters = "ensembl_gene_id",
                        values = original_rownames,
                        mart = ensembl)
      matching     <- gene_map$external_gene_name[match(original_rownames, gene_map$ensembl_gene_id)]
      new_rownames <- ifelse(is.na(matching) | matching == "", original_rownames, matching)
      new_rownames <- make.unique(new_rownames)
      rownames(obj[["RNA"]]) <- new_rownames
      cat("Gene conversion completed, number of genes:", length(new_rownames), "\n")
      cat("Example rownames after conversion:", head(new_rownames, 5), "\n")
    }, error = function(e) {
      warning("biomaRt query failed, using original gene names: ", e$message)
      cat("biomaRt query failed:", e$message, "\n")
    })
  }
  
  # Quality control
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-|^mt-")
  if (all(is.na(obj[["percent.mt"]]) | obj[["percent.mt"]] == 0)) {
    warning("No mitochondrial genes detected, percent.mt set to 0")
    obj[["percent.mt"]] <- 0
    cat("No mitochondrial genes detected, percent.mt set to 0\n")
  }
  cat("percent.mt distribution:", summary(obj$percent.mt), "\n")
  
  obj <- subset(obj,
                subset = nFeature_RNA >= qc_params$min_nFeature_RNA &
                  nFeature_RNA <= qc_params$max_nFeature_RNA &
                  percent.mt   <  qc_params$max_percent_mt   &
                  nCount_RNA   >= qc_params$min_nCount_RNA)
  
  cat("Number of cells after quality control:", ncol(obj), "\n")
  if (ncol(obj) < qc_params$min_cells) {
    warning("Dataset ", dataset_name, " has insufficient cells after QC (", ncol(obj), "), skipping")
    return(NULL)
  }
  
  # Normalization
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE)
  data_data <- GetAssayData(obj, assay = "RNA", layer = "data")
  if (all(dim(data_data) == 0)) {
    warning("data layer empty after normalization, using counts layer")
    obj <- SetAssayData(obj, assay = "RNA", layer = "data",
                        new.data = GetAssayData(obj, assay = "RNA", layer = "counts"))
  }
  cat("Normalization completed, data layer dimension:", dim(data_data), "\n")
  cat("data layer data distribution:", summary(as.vector(data_data)), "\n")
  
  # UMAP
  if (!"umap" %in% names(obj@reductions)) {
    cat("UMAP reduction not found, running UMAP...\n")
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, npcs = 20, verbose = FALSE)
    obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE, k.param = min(20, ncol(obj) - 1))
    obj <- FindClusters(obj, verbose = FALSE)
    obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
    cat("UMAP reduction completed\n")
  }
  
  # SingleR annotation
  if (!"cell_type" %in% colnames(obj@meta.data)) {
    cat("Using SingleR for automatic cell type annotation (HumanPrimaryCellAtlasData)...\n")
    tryCatch({
      hpca_ref      <- HumanPrimaryCellAtlasData()
      expr_matrix     <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data"))
      singler_results <- SingleR(test = expr_matrix, ref = hpca_ref, labels = hpca_ref$label.main)
      obj@meta.data$cell_type <- singler_results$labels
      cat("Cell type distribution:\n")
      print(table(obj@meta.data$cell_type))
    }, error = function(e) {
      warning("SingleR annotation failed: ", e$message)
      cat("SingleR annotation failed:", e$message, "\n")
    })
  } else {
    cat("cell_type annotation already exists, skipping SingleR annotation\n")
  }
  
  # ====================== Module Scores + Wilcoxon ======================
  cat("\n=== Calculating Module Scores + Wilcoxon tests ===\n")
  
  module_definitions <- list(
    "GZMB"    = gsub("_", "-", readLines(file.path(data_dir, "Target Genes GZMB.txt"))),
    "LAMP3"   = gsub("_", "-", readLines(file.path(data_dir, "Target Genes LAMP3.txt"))),
    "NKG7"    = gsub("_", "-", readLines(file.path(data_dir, "Target Genes NKG7.txt"))),
    "SLC39A8" = gsub("_", "-", readLines(file.path(data_dir, "Target Genes SLC39A8.txt")))
  )
  
  for (mod_name in names(module_definitions)) {
    genes <- module_definitions[[mod_name]][module_definitions[[mod_name]] %in% rownames(obj[["RNA"]])]
    if (length(genes) > 0) {
      obj <- AddModuleScore(obj, features = list(genes), name = mod_name, assay = "RNA")
      default_col <- paste0(mod_name, "1")
      if (default_col %in% colnames(obj@meta.data)) {
        colnames(obj@meta.data)[colnames(obj@meta.data) == default_col] <- mod_name
      }
      cat("✓ Module score added:", mod_name, "(", length(genes), "genes)\n")
    }
  }
  
  Idents(obj) <- obj$cell_type
  module_scores  <- c("GZMB", "LAMP3", "NKG7", "SLC39A8")
  wilcox_results <- data.frame()
  cat("Performing Wilcoxon rank-sum tests...\n")
  
  for (mod in module_scores) {
    if (!mod %in% colnames(obj@meta.data)) { cat("Warning:", mod, "not found, skipping.\n"); next }
    for (cell_type in unique(obj$cell_type)) {
      cells_in_type <- WhichCells(obj, idents = cell_type)
      cells_other   <- setdiff(WhichCells(obj), cells_in_type)
      score_in      <- obj@meta.data[cells_in_type, mod]
      score_other   <- obj@meta.data[cells_other,   mod]
      test_result   <- wilcox.test(score_in, score_other, alternative = "two.sided")
      median_in     <- median(score_in, na.rm = TRUE)
      p_val         <- test_result$p.value
      wilcox_results <- rbind(wilcox_results, data.frame(
        Module       = mod,
        CellType     = cell_type,
        Median_Score = round(median_in, 3),
        P_Value      = p_val,
        log10_P      = ifelse(p_val > 0, round(-log10(p_val), 4), 300)
      ))
    }
  }
  wilcox_results$adj_P <- p.adjust(wilcox_results$P_Value, method = "bonferroni")
  
  results_dir <- file.path(base_output_dir, "wilcoxon_results")
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(wilcox_results[order(wilcox_results$P_Value), ],
            file.path(results_dir, "wilcoxon_module_celltype_enrichment.csv"),
            row.names = FALSE)
  cat("Wilcoxon tests completed!\n")
  cat("Results saved to:", file.path(results_dir, "wilcoxon_module_celltype_enrichment.csv"), "\n")
  print(head(wilcox_results[order(wilcox_results$P_Value), ], 10))
  
  saveRDS(obj, preprocessed_file)
  cat("Preprocessed Seurat object saved to:", preprocessed_file, "\n")
  return(obj)
}

# ========================================================
# Gene-specific analysis function
# ========================================================
analyze_genes <- function(obj, dataset_name, output_dir, target_genes, reference_gene = NULL) {
  log_file <- file.path(output_dir, "analysis_log.txt")
  cat("Analyzing gene group:", ifelse(is.null(reference_gene), "Target_Genes", reference_gene), "\n")
  
  tryCatch({
    cat("Attempting to create and write log file:", log_file, "\n")
    cat(paste("Analysis start time:", Sys.time(), "\n"),
        paste("Dataset:", dataset_name, "\n"),
        paste("Quality control parameters:\n"),
        paste(" min_nFeature_RNA:", qc_params$min_nFeature_RNA, "\n"),
        paste(" max_nFeature_RNA:", qc_params$max_nFeature_RNA, "\n"),
        paste(" max_percent_mt:",  qc_params$max_percent_mt,    "\n"),
        paste(" min_nCount_RNA:",  qc_params$min_nCount_RNA,    "\n"),
        paste(" min_cells:",       qc_params$min_cells,         "\n"),
        paste("Total target genes:", length(target_genes), "\n"),
        paste("First 5 target genes:", paste(head(target_genes, 5), collapse = ", "), "\n"),
        if (!is.null(reference_gene)) paste("Reference gene:", reference_gene, "\n") else "",
        file = log_file, append = FALSE)
  }, error = function(e) {
    cat("Unable to write log file:", e$message, "\n")
    log_file <<- NULL
  })
  
  tryCatch({
    # ---- Available genes ----
    available_genes <- target_genes[target_genes %in% rownames(GetAssayData(obj, assay = "RNA", layer = "data"))]
    if (!is.null(reference_gene) && reference_gene %in% rownames(GetAssayData(obj, assay = "RNA", layer = "data"))) {
      available_genes <- unique(c(reference_gene, available_genes))
    }
    cat("Number of target genes present in the dataset:", length(available_genes), "\n")
    cat("Target genes present in the dataset:", paste(available_genes, collapse = ", "), "\n")
    missing_genes <- setdiff(target_genes, available_genes)
    if (length(missing_genes) > 0) {
      cat("Missing target genes:", paste(missing_genes, collapse = ", "), "\n")
      if (!is.null(log_file)) cat(paste("Missing target genes:", paste(missing_genes, collapse = ", "), "\n"), file = log_file, append = TRUE)
    }
    if (!is.null(log_file)) {
      cat(paste("Number of target genes present in the dataset:", length(available_genes), "\n"),
          paste("Target genes present in the dataset:", paste(available_genes, collapse = ", "), "\n"),
          file = log_file, append = TRUE)
    }
    
    available_genes_file <- file.path(output_dir, "available_genes.txt")
    writeLines(available_genes, available_genes_file)
    cat("Available target gene list saved to:", available_genes_file, "\n")
    
    if (length(available_genes) == 0) {
      warning("Dataset ", dataset_name, " does not contain any target genes")
      cat(paste("Dataset", dataset_name, "does not contain any target genes, skipping analysis\n"))
      return(NULL)
    }
    
    # ---- Expression matrix ----
    expr_matrix <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[available_genes, , drop = FALSE])
    write.csv(expr_matrix, file.path(output_dir, "target_genes_expression.csv"), row.names = TRUE)
    cat("Target gene expression matrix saved to:", file.path(output_dir, "target_genes_expression.csv"), "\n")
    
    # ---- Correlation heatmap ----
    if (length(available_genes) > 1) {
      tryCatch({
        expr_matrix_clean <- expr_matrix
        expr_matrix_clean[is.na(expr_matrix_clean) | is.nan(expr_matrix_clean) | is.infinite(expr_matrix_clean)] <- 0
        gene_vars   <- apply(expr_matrix_clean, 1, var, na.rm = TRUE)
        valid_genes <- rownames(expr_matrix_clean)[gene_vars > 0 & !is.na(gene_vars)]
        if (length(valid_genes) >= 2) {
          expr_matrix_clean <- expr_matrix_clean[valid_genes, , drop = FALSE]
          cor_matrix <- cor(t(expr_matrix_clean), method = "spearman", use = "pairwise.complete.obs")
          write.csv(cor_matrix, file.path(output_dir, "gene_correlation_matrix.csv"), row.names = TRUE)
          cat("Inter-gene correlation matrix (Spearman) saved to:", file.path(output_dir, "gene_correlation_matrix.csv"), "\n")
          suppressWarnings(pheatmap(
            cor_matrix,
            filename        = file.path(output_dir, "correlation_heatmap.png"),
            width           = 10, height = 8,
            main            = "Spearman Correlation of Target Genes",
            color           = colorRampPalette(c("blue", "white", "red"))(100),
            breaks          = seq(-1, 1, length.out = 101),
            fontsize        = PHEATMAP_FONTSIZE,
            fontsize_row    = PHEATMAP_FONTSIZE_ROW,
            fontsize_col    = PHEATMAP_FONTSIZE_COL,
            silent          = TRUE
          ))
          cat("Correlation heatmap saved to:", file.path(output_dir, "correlation_heatmap.png"), "\n")
        }
        cat("⚡ Skipped full-cell expression heatmap to avoid memory error.\n")
      }, error = function(e) {
        cat("Failed to compute inter-gene correlation or heatmap:", e$message, "\n")
      })
    }
    
    # ---- DotPlot (available_genes) ----
    if ("cell_type" %in% colnames(obj@meta.data)) {
      tryCatch({
        dot_plot <- DotPlot(obj, features = available_genes, group.by = "cell_type") +
          theme(
            axis.text.x  = element_text(size = AXIS_TEXT_SIZE, angle = AXIS_TEXT_X_ANGLE, hjust = 1),
            axis.text.y  = element_text(size = AXIS_TEXT_SIZE),
            text         = element_text(size = PLOT_FONTSIZE)
          ) +
          labs(title = "Expression of Target Genes by Cell Type")
        # [MODIFIED] Triple-format output via save_fig
        save_fig(dot_plot, file.path(output_dir, "dotplot_target_genes"), w_in = 12, h_in = 8)
        cat("Target gene dot plot saved (.pdf/_300.png/_600.tiff)\n")
      }, error = function(e) { cat("Failed to generate dot plot:", e$message, "\n") })
    }
    
    # ---- Expression summary ----
    tryCatch({
      expr_summary <- data.frame(
        Gene              = available_genes,
        Mean_Expression   = apply(expr_matrix, 1, mean,   na.rm = TRUE),
        Median_Expression = apply(expr_matrix, 1, median, na.rm = TRUE),
        NonZero_Prop      = apply(expr_matrix, 1, function(x) mean(x > 0, na.rm = TRUE)),
        Variance          = apply(expr_matrix, 1, var,    na.rm = TRUE)
      )
      write.table(expr_summary, file.path(output_dir, "target_genes_summary.txt"),
                  row.names = FALSE, sep = "\t", quote = FALSE)
      cat("Target gene expression statistics summary saved to:", file.path(output_dir, "target_genes_summary.txt"), "\n")
      print(head(expr_summary))
    }, error = function(e) { cat("Failed to generate expression summary:", e$message, "\n") })
    
    # ---- Sparsity + Correlation (reference gene) ----
    if (!is.null(reference_gene) && reference_gene %in% rownames(GetAssayData(obj, assay = "RNA", layer = "data"))) {
      cat(paste("Calculating sparsity of target genes and", reference_gene, "...\n"))
      valid_genes <- c(reference_gene, target_genes)[c(reference_gene, target_genes) %in%
                                                       rownames(GetAssayData(obj, assay = "RNA", layer = "data"))]
      if (length(valid_genes) > 0) {
        expr_ref <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[valid_genes, , drop = FALSE])
        sparsity_stats <- data.frame(
          Gene         = valid_genes,
          NonZero_Prop = apply(expr_ref, 1, function(x) mean(x > 0, na.rm = TRUE)),
          Variance     = apply(expr_ref, 1, var, na.rm = TRUE)
        )
        write.csv(sparsity_stats,
                  file.path(output_dir, paste0(tolower(reference_gene), "_sparsity_stats.csv")),
                  row.names = FALSE)
        cat("Sparsity statistics saved to:", file.path(output_dir, paste0(tolower(reference_gene), "_sparsity_stats.csv")), "\n")
        
        # Sparsity by cell type
        if ("cell_type" %in% colnames(obj@meta.data)) {
          cell_types <- unique(obj@meta.data$cell_type)
          sparsity_by_celltype <- lapply(cell_types, function(ct) {
            cells <- colnames(obj)[obj@meta.data$cell_type == ct]
            if (length(cells) > 10) {
              expr_ct <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[valid_genes, cells, drop = FALSE])
              data.frame(Gene = valid_genes, Cell_Type = ct,
                         NonZero_Prop = apply(expr_ct, 1, function(x) mean(x > 0, na.rm = TRUE)),
                         Variance     = apply(expr_ct, 1, var, na.rm = TRUE))
            } else NULL
          })
          sparsity_by_celltype <- do.call(rbind, sparsity_by_celltype[!sapply(sparsity_by_celltype, is.null)])
          write.csv(sparsity_by_celltype,
                    file.path(output_dir, paste0(tolower(reference_gene), "_sparsity_by_celltype.csv")),
                    row.names = FALSE)
          cat("Sparsity by cell type saved to:", file.path(output_dir, paste0(tolower(reference_gene), "_sparsity_by_celltype.csv")), "\n")
        }
        
        # Correlation analysis
        cat(paste("Starting correlation analysis between target genes and", reference_gene, "...\n"))
        expr_ref <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[valid_genes, , drop = FALSE])
        expr_ref[is.na(expr_ref) | is.nan(expr_ref) | is.infinite(expr_ref)] <- 0
        expr_ref <- expr_ref[apply(expr_ref, 1, function(x) var(x, na.rm = TRUE) > 0), , drop = FALSE]
        
        if ("cell_type" %in% colnames(obj@meta.data)) {
          for (ct in unique(obj@meta.data$cell_type)) {
            cells <- colnames(obj)[obj@meta.data$cell_type == ct]
            if (length(cells) > 10) {
              expr_ct <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[valid_genes, cells, drop = FALSE])
              expr_ct <- expr_ct[apply(expr_ct, 1, function(x) var(x, na.rm = TRUE) > 0), , drop = FALSE]
              if (nrow(expr_ct) >= 2 && reference_gene %in% rownames(expr_ct)) {
                cor_ct <- data.frame(
                  Gene        = valid_genes[valid_genes != reference_gene],
                  Correlation = sapply(valid_genes[valid_genes != reference_gene], function(g) {
                    if (g %in% rownames(expr_ct))
                      suppressWarnings(cor.test(expr_ct[reference_gene, ], expr_ct[g, ], method = "spearman")$estimate)
                    else NA
                  }),
                  P_Value     = sapply(valid_genes[valid_genes != reference_gene], function(g) {
                    if (g %in% rownames(expr_ct))
                      suppressWarnings(cor.test(expr_ct[reference_gene, ], expr_ct[g, ], method = "spearman")$p.value)
                    else NA
                  })
                )
                cor_ct_file <- file.path(output_dir, paste0(tolower(reference_gene), "_correlation_", gsub(" ", "_", ct), ".csv"))
                write.csv(cor_ct, cor_ct_file, row.names = FALSE)
                cat("Correlation (", ct, ") saved to:", cor_ct_file, "\n")
              }
            }
          }
          
          # Sparse gene correlation
          sparsity_threshold  <- 0.5
          sparsity_data       <- read.csv(file.path(output_dir, paste0(tolower(reference_gene), "_sparsity_by_celltype.csv")))
          high_ref_cell_types <- sparsity_data$Cell_Type[sparsity_data$Gene == reference_gene & sparsity_data$NonZero_Prop > 0.5]
          for (ct in high_ref_cell_types) {
            cells <- colnames(obj)[obj@meta.data$cell_type == ct]
            if (length(cells) > 10) {
              expr_ct <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[valid_genes, cells, drop = FALSE])
              expr_ct <- expr_ct[apply(expr_ct, 1, function(x) var(x, na.rm = TRUE) > 0), , drop = FALSE]
              if (nrow(expr_ct) >= 2 && reference_gene %in% rownames(expr_ct)) {
                ct_sparsity  <- sparsity_data[sparsity_data$Cell_Type == ct & sparsity_data$Gene != reference_gene, ]
                sparse_genes <- ct_sparsity$Gene[ct_sparsity$NonZero_Prop < sparsity_threshold & ct_sparsity$Gene %in% rownames(expr_ct)]
                if (length(sparse_genes) > 0) {
                  expr_sparse_ct <- expr_ct[c(reference_gene, sparse_genes), , drop = FALSE]
                  cor_ct <- data.frame(
                    Gene        = sparse_genes,
                    Correlation = sapply(sparse_genes, function(g)
                      suppressWarnings(cor.test(expr_sparse_ct[reference_gene, ], expr_sparse_ct[g, ], method = "spearman")$estimate)),
                    P_Value     = sapply(sparse_genes, function(g)
                      suppressWarnings(cor.test(expr_sparse_ct[reference_gene, ], expr_sparse_ct[g, ], method = "spearman")$p.value))
                  )
                  cor_ct_file <- file.path(output_dir, paste0(tolower(reference_gene), "_correlation_", gsub(" ", "_", ct), "_sparse.csv"))
                  write.csv(cor_ct, cor_ct_file, row.names = FALSE)
                  cat("Sparse correlation (", ct, ") saved to:", cor_ct_file, "\n")
                }
              }
            }
          }
        }
        
        # Correlation in reference-positive cells
        ref_positive_cells <- colnames(obj)[expr_ref[reference_gene, ] > 0]
        if (length(ref_positive_cells) > 10) {
          expr_ref_pos <- expr_ref[, ref_positive_cells, drop = FALSE]
          expr_ref_pos <- expr_ref_pos[apply(expr_ref_pos, 1, function(x) var(x, na.rm = TRUE) > 0), , drop = FALSE]
          if (nrow(expr_ref_pos) >= 2 && reference_gene %in% rownames(expr_ref_pos)) {
            cor_pos <- data.frame(
              Gene            = valid_genes[valid_genes != reference_gene],
              Correlation_Pos = sapply(valid_genes[valid_genes != reference_gene], function(g) {
                if (g %in% rownames(expr_ref_pos))
                  suppressWarnings(cor.test(expr_ref_pos[reference_gene, ], expr_ref_pos[g, ], method = "spearman")$estimate)
                else NA
              }),
              P_Value_Pos     = sapply(valid_genes[valid_genes != reference_gene], function(g) {
                if (g %in% rownames(expr_ref_pos))
                  suppressWarnings(cor.test(expr_ref_pos[reference_gene, ], expr_ref_pos[g, ], method = "spearman")$p.value)
                else NA
              })
            )
            write.csv(cor_pos,
                      file.path(output_dir, paste0(tolower(reference_gene), "_correlation_positive_cells.csv")),
                      row.names = FALSE)
            cat("Positive-cell correlation saved to:", file.path(output_dir, paste0(tolower(reference_gene), "_correlation_positive_cells.csv")), "\n")
          }
        }
        
        # Reference-gene DotPlot
        if ("cell_type" %in% colnames(obj@meta.data)) {
          tryCatch({
            dot_plot <- DotPlot(obj, features = valid_genes, group.by = "cell_type") +
              theme(
                axis.text.x  = element_text(size = AXIS_TEXT_SIZE, angle = AXIS_TEXT_X_ANGLE, hjust = 1),
                axis.text.y  = element_text(size = AXIS_TEXT_SIZE),
                text         = element_text(size = PLOT_FONTSIZE)
              ) +
              labs(title = paste("Expression of Genes with", reference_gene, "by Cell Type"))
            # [MODIFIED] Triple-format output via save_fig
            save_fig(dot_plot, file.path(output_dir, paste0(tolower(reference_gene), "_dotplot")), w_in = 12, h_in = 8)
            cat("Dot plot saved (.pdf/_300.png/_600.tiff)\n")
          }, error = function(e) { cat("Failed to generate dot plot:", e$message, "\n") })
        }
        
        # Reference-gene correlation heatmap
        if (nrow(expr_ref) > 1) {
          cor_ref <- cor(t(expr_ref), method = "spearman", use = "pairwise.complete.obs")
          tryCatch({
            pheatmap(
              cor_ref,
              filename     = file.path(output_dir, paste0(tolower(reference_gene), "_correlation_heatmap.png")),
              width        = 8, height = 6,
              main         = paste("Spearman Correlation: Genes vs", reference_gene),
              color        = colorRampPalette(c("blue", "white", "red"))(100),
              breaks       = seq(-1, 1, length.out = 101),
              fontsize     = PHEATMAP_FONTSIZE,
              fontsize_row = PHEATMAP_FONTSIZE_ROW,
              fontsize_col = PHEATMAP_FONTSIZE_COL,
              silent       = TRUE
            )
            cat("Correlation heatmap saved to:", file.path(output_dir, paste0(tolower(reference_gene), "_correlation_heatmap.png")), "\n")
          }, error = function(e) { cat("Failed to generate correlation heatmap:", e$message, "\n") })
        }
      }
    }
    
    # ---- Violin plots (batched) ----
    if (length(available_genes) > 0) {
      tryCatch({
        num_genes  <- length(available_genes)
        batch_size <- 9
        for (b in seq(1, num_genes, batch_size)) {
          batch_genes <- available_genes[b:min(b + batch_size - 1, num_genes)]
          vln_data <- data.frame()
          for (gene in batch_genes) {
            expr <- GetAssayData(obj, assay = "RNA", layer = "data")[gene, ]
            vln_data <- rbind(vln_data, data.frame(
              Cell       = names(expr),
              Gene       = gene,
              Expression = as.numeric(expr),
              Cell_Type  = obj@meta.data$cell_type[match(names(expr), rownames(obj@meta.data))]
            ))
          }
          write.csv(vln_data,
                    file.path(output_dir, paste0("violin_data_batch_", ceiling(b/batch_size), ".csv")),
                    row.names = FALSE)
          cat("Violin plot data saved to:", file.path(output_dir, paste0("violin_data_batch_", ceiling(b/batch_size), ".csv")), "\n")
          
          violin_plots <- VlnPlot(obj, features = batch_genes, group.by = "cell_type", pt.size = 0) &
            theme(
              axis.text.x = element_text(size = VIOLIN_X_TEXT_SIZE, angle = AXIS_TEXT_X_ANGLE, hjust = 1, vjust = 1),
              axis.text.y = element_text(size = VIOLIN_Y_TEXT_SIZE),
              axis.title.x = element_blank(),
              text         = element_text(size = PLOT_FONTSIZE)
            )
          # [MODIFIED] Triple-format output via save_fig
          save_fig(violin_plots, file.path(output_dir, paste0("violin_plots_target_genes_batch_", ceiling(b/batch_size))), w_in = 14, h_in = 10)
        }
        cat("Violin plots saved (in batches)\n")
      }, error = function(e) { cat("Failed to generate violin plots:", e$message, "\n") })
      
      # ---- UMAP feature plots (batched) ----
      tryCatch({
        num_genes  <- length(available_genes)
        batch_size <- 9
        for (b in seq(1, num_genes, batch_size)) {
          batch_genes <- available_genes[b:min(b + batch_size - 1, num_genes)]
          umap_coords <- as.data.frame(obj@reductions$umap@cell.embeddings)
          umap_coords$Cell <- rownames(umap_coords)
          expr_batch  <- as.data.frame(t(as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[batch_genes, , drop = FALSE])))
          expr_batch$Cell <- rownames(expr_batch)
          umap_data   <- merge(umap_coords, expr_batch, by = "Cell")
          write.csv(umap_data,
                    file.path(output_dir, paste0("umap_data_batch_", ceiling(b/batch_size), ".csv")),
                    row.names = FALSE)
          cat("UMAP feature plot data saved to:", file.path(output_dir, paste0("umap_data_batch_", ceiling(b/batch_size), ".csv")), "\n")
          
          feature_plots <- FeaturePlot(obj, features = batch_genes, combine = TRUE) &
            theme(
              axis.text  = element_text(size = AXIS_TEXT_SIZE),
              text       = element_text(size = PLOT_FONTSIZE)
            )
          # [MODIFIED] Triple-format output via save_fig
          save_fig(feature_plots, file.path(output_dir, paste0("UMAP_feature_plots_batch_", ceiling(b/batch_size))), w_in = 12, h_in = 8)
        }
        cat("UMAP feature plots saved (in batches)\n")
      }, error = function(e) { cat("Failed to generate UMAP feature plots:", e$message, "\n") })
    }
    
    # ---- GO/KEGG enrichment ----
    tryCatch({
      gene_entrez <- bitr(available_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
      if (length(gene_entrez) == 0) {
        cat("No genes could be mapped to Entrez IDs, skipping enrichment.\n")
      } else {
        ego <- enrichGO(gene_entrez, OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 0.05)
        
        kegg_link_file <- file.path(data_dir, "kegg_jp_link_hsa_pathway.txt")
        if (file.exists(kegg_link_file)) {
          cat("✓ Using LOCAL KEGG file:", kegg_link_file, "\n")
          kegg_raw         <- read.delim(kegg_link_file, header = FALSE, sep = "\t",
                                         stringsAsFactors = FALSE, col.names = c("Pathway", "Gene"))
          kegg_raw$Pathway <- sub("^path:", "", kegg_raw$Pathway)
          kegg_raw$Gene    <- sub("^hsa:", "", kegg_raw$Gene)
          TERM2GENE        <- kegg_raw[, c("Pathway", "Gene")]
          colnames(TERM2GENE) <- c("term", "gene")
          ekg <- enricher(gene = gene_entrez, TERM2GENE = TERM2GENE,
                          pvalueCutoff = 0.05, pAdjustMethod = "BH",
                          minGSSize = 10, maxGSSize = 500)
          cat("✓ Local KEGG enrichment completed (", if (!is.null(ekg)) nrow(ekg@result) else 0, "terms found)\n")
        } else {
          warning("Local KEGG file not found, falling back to online")
          ekg <- enrichKEGG(gene_entrez, organism = "hsa", pvalueCutoff = 0.05)
        }
        
        if (!is.null(ego) && nrow(ego) > 0) {
          write.csv(ego@result, file.path(output_dir, "go_enrichment_data.csv"), row.names = FALSE)
          cat("GO enrichment data saved to:", file.path(output_dir, "go_enrichment_data.csv"), "\n")
          dotplot(ego, showCategory = 10) +
            theme(
              axis.text.x = element_text(size = AXIS_TEXT_SIZE),
              axis.text.y = element_text(size = AXIS_TEXT_SIZE),
              text        = element_text(size = PLOT_FONTSIZE)
            )
          # [MODIFIED] Triple-format output via save_fig
          save_fig(last_plot(), file.path(output_dir, "go_bubble"), w_in = 8, h_in = 6)
          cat("GO enrichment completed and saved.\n")
        }
        if (!is.null(ekg) && nrow(ekg@result) > 0) {
          write.csv(ekg@result, file.path(output_dir, "kegg_enrichment_data.csv"), row.names = FALSE)
          cat("KEGG enrichment results saved to:", file.path(output_dir, "kegg_enrichment_data.csv"), "\n")
          dotplot(ekg, showCategory = 10) +
            theme(
              axis.text.x = element_text(size = AXIS_TEXT_SIZE),
              axis.text.y = element_text(size = AXIS_TEXT_SIZE),
              text        = element_text(size = PLOT_FONTSIZE)
            )
          # [MODIFIED] Triple-format output via save_fig
          save_fig(last_plot(), file.path(output_dir, "kegg_bubble"), w_in = 8, h_in = 6)
        }
      }
    }, error = function(e) { cat("GO/KEGG enrichment analysis failed:", e$message, "\n") })
    
    # ---- PPI ----
    tryCatch({
      string_db <- STRINGdb$new(species = 9606, score_threshold = 500, version = "12.0", input_directory = data_dir)
      if (exists("graph", envir = string_db) && is.igraph(string_db$graph)) {
        string_db$graph <- upgrade_graph(string_db$graph)
        cat("STRINGdb internal igraph object upgraded\n")
      }
      mapped <- string_db$map(data.frame(gene = available_genes), "gene", removeUnmappedRows = TRUE)
      if (nrow(mapped) < 2) {
        cat("Insufficient mapped genes for PPI analysis.\n")
      } else {
        interactions <- string_db$get_interactions(mapped$STRING_id)
        if (nrow(interactions) == 0) {
          cat("No interactions found.\n")
        } else {
          from_genes <- mapped$gene[match(interactions$from, mapped$STRING_id)]
          to_genes   <- mapped$gene[match(interactions$to,   mapped$STRING_id)]
          joint_expr <- mean(expr_matrix[from_genes, ] > 0 & expr_matrix[to_genes, ] > 0, na.rm = TRUE)
          cat("PPI joint expression rate:", joint_expr, "\n")
          write.csv(interactions, file.path(output_dir, "ppi_interactions.csv"), row.names = FALSE)
          cat("PPI interaction data saved to:", file.path(output_dir, "ppi_interactions.csv"), "\n")
        }
      }
    }, error = function(e) { cat("PPI analysis failed:", e$message, "\n") })
    
    # ================================================================
    # 5. Monocle3 pseudotime trajectory analysis
    # ================================================================
    tryCatch({
      gene_annotation <- data.frame(
        gene_short_name = rownames(obj[["RNA"]]),
        row.names       = rownames(obj[["RNA"]])
      )
      cds <- new_cell_data_set(
        expression_data = GetAssayData(obj, assay = "RNA", layer = "counts"),
        cell_metadata   = obj@meta.data,
        gene_metadata   = gene_annotation
      )
      cds <- preprocess_cds(cds, method = "PCA")
      cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")
      cds <- cluster_cells(cds, reduction_method = "UMAP")
      cat("Cluster partitions completed for Monocle3\n")
      if (is.null(partitions(cds))) stop("No partitions found after cluster_cells")
      cds <- learn_graph(cds)
      
      uvmr_genes <- intersect(available_genes, rownames(cds))
      cat("uvmr_genes (filtered):", paste(head(uvmr_genes, 5), collapse = ", "), "\n")
      
      cat("=== Pseudotime debug ===\n")
      cat("uvmr_genes count:", length(uvmr_genes), "\n")
      cat("uvmr_genes:", paste(uvmr_genes, collapse = ", "), "\n")
      cat("cell_type in colData:", "cell_type" %in% colnames(colData(cds)), "\n")
      cat("partitions null:", is.null(partitions(cds)), "\n")
      
      if (length(uvmr_genes) == 0) {
        cat("uvmr_genes is empty, skipping pseudotime analysis\n")
      } else {

        # ================================================================
        # Interactive root selection: user clicks root node on UMAP plot
        # ================================================================
        cds <- order_cells(cds)
        cat("✓ Pseudotime ordered (interactive root selection)\n")
        
        pt_check <- pseudotime(cds)
        cat("Pseudotime Inf count:",   sum(is.infinite(pt_check)), "\n")
        cat("Pseudotime valid count:", sum(!is.infinite(pt_check) & !is.na(pt_check)), "\n")
        
        # ---- Bootstrap CI ----
        tryCatch({
          pt_raw   <- pseudotime(cds)
          valid_pt <- !is.infinite(pt_raw) & !is.na(pt_raw)
          pt_norm  <- (pt_raw[valid_pt] - min(pt_raw[valid_pt])) / diff(range(pt_raw[valid_pt]))
          calc_pct <- function(ev, pv, e = c(0, 0.3), l = c(0.7, 1.0)) {
            me <- median(ev[pv >= e[1] & pv <= e[2]], na.rm = TRUE)
            ml <- median(ev[pv >= l[1] & pv <= l[2]], na.rm = TRUE)
            if (is.na(me) || abs(me) < 1e-10) return(NA)
            (ml - me) / abs(me) * 100
          }
          pt_ci_stats <- do.call(rbind, lapply(uvmr_genes, function(g) {
            if (!g %in% rownames(cds)) return(NULL)
            ev  <- as.numeric(exprs(cds[g, valid_pt]))
            obs <- calc_pct(ev, pt_norm)
            bts <- replicate(200, {
              idx <- sample(length(ev), replace = TRUE)
              calc_pct(ev[idx], pt_norm[idx])
            })
            data.frame(
              Gene       = g,
              Dataset    = "GSE296117",
              Pct_Change = round(obs, 1),
              CI_lower   = round(quantile(bts, 0.025, na.rm = TRUE), 1),
              CI_upper   = round(quantile(bts, 0.975, na.rm = TRUE), 1),
              CV_pct     = round(sd(bts, na.rm = TRUE) / abs(mean(bts, na.rm = TRUE)) * 100, 1)
            )
          }))
          write.csv(pt_ci_stats, file.path(output_dir, "pseudotime_change_CI.csv"), row.names = FALSE)
          cat("Pseudotime bootstrap CI saved to:", file.path(output_dir, "pseudotime_change_CI.csv"), "\n")
        }, error = function(e) cat("Bootstrap CI failed:", e$message, "\n"))
        
        # ---- Save pseudotime data ----
        pseudotime_vals        <- pseudotime(cds)
        expr_pseudo            <- as.data.frame(t(as.matrix(exprs(cds[uvmr_genes, ]))))
        expr_pseudo$Cell       <- rownames(expr_pseudo)
        expr_pseudo$Pseudotime <- pseudotime_vals[match(expr_pseudo$Cell, colnames(cds))]
        write.csv(expr_pseudo, file.path(output_dir, "pseudotime_data.csv"), row.names = FALSE)
        cat("Pseudotime data saved to:", file.path(output_dir, "pseudotime_data.csv"), "\n")
        if (!is.null(log_file)) cat(paste("Pseudotime data saved to:", file.path(output_dir, "pseudotime_data.csv"), "\n"), file = log_file, append = TRUE)
        
        # ---- Pseudotime plot ----
        p <- tryCatch(
          plot_genes_in_pseudotime(cds[uvmr_genes, ]) +
            theme(
              axis.text.x = element_text(size = AXIS_TEXT_SIZE),
              axis.text.y = element_text(size = PSEUDOTIME_Y_TEXT_SIZE),
              text        = element_text(size = PLOT_FONTSIZE)
            ),
          error = function(e) { cat("plot_genes_in_pseudotime failed:", e$message, "\n"); NULL }
        )
        if (!is.null(p)) {
          # [MODIFIED] Triple-format output via save_fig
          save_fig(p, file.path(output_dir, "pseudotime_plot"), w_in = 8, h_in = 6)
          cat("Pseudotime trajectory plot saved (.pdf/_300.png/_600.tiff)\n")
          if (!is.null(log_file)) cat(paste("Pseudotime trajectory plot saved to:", file.path(output_dir, "pseudotime_plot.png"), "\n"), file = log_file, append = TRUE)
        } else {
          cat("Pseudotime plot skipped due to error in plot_genes_in_pseudotime\n")
        }
      }
    }, error = function(e) {
      cat("Pseudotime analysis failed:", e$message, "\n")
      if (!is.null(log_file)) cat(paste("Pseudotime analysis failed:", e$message, "\n"), file = log_file, append = TRUE)
    })
    
    # ---- Save Seurat object ----
    saveRDS(obj, file.path(output_dir, "processed_seurat_obj.rds"))
    cat("Seurat object saved to:", file.path(output_dir, "processed_seurat_obj.rds"), "\n")
    
    cat(paste("Analysis completion time:", Sys.time(), "\n"),
        paste("Number of cells after quality control:", ncol(obj), "\n"),
        paste("Number of target genes:", length(available_genes), "\n"))
    if (!is.null(log_file)) {
      cat(paste("Analysis completion time:", Sys.time(), "\n"),
          paste("Number of cells:", ncol(obj), "\n"),
          paste("Number of target genes:", length(available_genes), "\n"),
          file = log_file, append = TRUE)
    }
    return(obj)
    
  }, error = function(e) {
    cat("Gene group", ifelse(is.null(reference_gene), "Target_Genes", reference_gene), "analysis failed:", e$message, "\n")
    if (!is.null(log_file)) {
      cat(paste("Analysis failure time:", Sys.time(), "\n"),
          paste("Error message:", e$message, "\n"),
          "Skipping this gene group\n",
          file = log_file, append = TRUE)
    }
    return(NULL)
  })
}

# ========================================================
# Main loop
# ========================================================
for (ds in datasets) {
  base_output_dir <- ds$base_output_dir
  cat("Checking base output directory:", base_output_dir, "\n")
  if (!dir.exists(base_output_dir)) {
    dir.create(base_output_dir, showWarnings = FALSE, recursive = TRUE)
    cat("Base output directory created:", base_output_dir, "\n")
  }
  cat("Reading dataset:", ds$name, "\n")
  tryCatch({
    obj_raw <- ds$load_func()
    cat("Successfully loaded", ds$name, "data, number of cells:", ncol(obj_raw), "\n")
    
    obj <- preprocess_dataset(obj = obj_raw, dataset_name = ds$name,
                              base_output_dir = base_output_dir, qc_params = qc_params)
    if (is.null(obj)) { cat("Preprocessing failed, skipping dataset\n"); next }
    
    for (gf in gene_files) {
      gene_list <- gsub("_", "-", readLines(gf$file))
      cat(paste("Processing file:", gf$file, ", total genes:", length(gene_list), "\n"))
      
      reference_gene <- if (gf$file != "Target Genes.txt" && length(gene_list) > 0) gene_list[1] else NULL
      target_genes   <- if (gf$file != "Target Genes.txt" && length(gene_list) > 1) gene_list[-1] else gene_list
      
      output_dir <- if (is.null(reference_gene)) {
        file.path(base_output_dir, "Target_Genes")
      } else {
        file.path(base_output_dir, reference_gene)
      }
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      cat("Output directory created:", output_dir, "\n")
      
      result <- analyze_genes(
        obj            = obj,
        dataset_name   = ds$name,
        output_dir     = output_dir,
        target_genes   = target_genes,
        reference_gene = reference_gene
      )
      
      if (is.null(result)) {
        cat(paste("Dataset", ds$name, "(", ifelse(is.null(reference_gene), "Target_Genes", reference_gene),
                  ") analysis failed, check log:", file.path(output_dir, "analysis_log.txt"), "\n"))
      } else {
        cat(paste("Dataset", ds$name, "(", ifelse(is.null(reference_gene), "Target_Genes", reference_gene),
                  ") analysis completed, output saved to:", output_dir, "\n"))
      }
    }
  }, error = function(e) {
    cat("Failed to read or preprocess", ds$name, "data:", e$message, "\n")
    global_log_file <- file.path(base_output_dir, "global_analysis_log.txt")
    tryCatch({
      cat(paste("Failed time:", Sys.time(), "\n"),
          paste("Error:", e$message, "\n"),
          "Skipping this dataset\n",
          file = global_log_file, append = TRUE)
    }, error = function(e2) cat("Unable to write global log:", e2$message, "\n"))
  })
}

cat("All dataset analysis scripts execution completed, time:", Sys.time(), "\n")
warnings()