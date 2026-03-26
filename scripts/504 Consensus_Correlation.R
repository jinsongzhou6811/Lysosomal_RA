# ========================================================
# Standard Path Management Template (All Modules)
# ========================================================
# install.packages("rprojroot") # Only execute once
library(rprojroot)
# Automatically detect project root directory (works on any computer, any drive)
project_root <- find_rstudio_root_file()
# Fixed paths (shared by all modules)
outputs_dir <- file.path(project_root, "outputs")
# ------------------- Input directory (key to modify) -------------------
# Location of previous analysis results: outputs/05_scRNA_Analysis/dataset/gene/
input_dir <- file.path(outputs_dir, "05_scRNA_Analysis")
# ------------------- Output directory (optional to modify) -------------------
# If you want to put the final results under outputs/05_scRNA_Analysis/Consensus_Correlation/
module_name <- "05_scRNA_Analysis\\Consensus_Correlation" # ← Keep this
# ------------------------------------------------------
# Actual module output directory (automatically created)
module_dir <- if (module_name != "") {
  file.path(outputs_dir, module_name)
} else {
  outputs_dir
}
dir.create(module_dir, showWarnings = FALSE, recursive = TRUE)
# ========================================================
# Cross-dataset consensus correlation heatmap merging script (rprojroot + new input path version)
# ========================================================
library(readr)
library(dplyr)
library(tibble)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
# Genes and datasets to process
genes <- c("GZMB", "LAMP3", "NKG7", "TRAF3IP3")
datasets <- c("GSE159117", "GSE202375", "GSE296117")
# ==================== Core function ====================
compute_merged_matrix <- function(gene, input_dir) {
  files <- file.path(input_dir, datasets, gene, "gene_correlation_matrix.csv")
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
    warning("Missing files (skip this module): ", paste(basename(missing), collapse = ", "))
    return(NULL)
  }
  mats <- lapply(files, function(f) {
    read_csv(f, col_types = cols()) %>%
      column_to_rownames(var = "...1") %>%
      as.matrix()
  })
  common_genes <- Reduce(intersect, lapply(mats, rownames))
  cat("Module", gene, "number of common genes:", length(common_genes), "\n")
  sub_mats <- lapply(mats, function(m) m[common_genes, common_genes])
  merged <- Reduce("+", sub_mats) / length(sub_mats)
  merged <- (merged + t(merged)) / 2
  diag(merged) <- 1
  return(list(matrix = merged, sub_mats = sub_mats, genes = common_genes))
}
# ==================== Main loop ====================
for (gene in genes) {
  result <- compute_merged_matrix(gene, input_dir) # ← Use the new input_dir
  if (is.null(result)) next
  mat <- result$matrix
  sub_mats <- result$sub_mats
  # Output directory: outputs/05_scRNA_Analysis/Consensus_Correlation/gene/
  out_dir <- if (module_name != "") {
    file.path(module_dir, gene)
  } else {
    file.path(outputs_dir, gene)
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  # Save consensus matrix
  write.csv(mat, file.path(out_dir, paste0(gene, "_consensus_correlation_matrix.csv")), row.names = TRUE)
  # Draw heatmap
  p <- pheatmap(mat,
                color = colorRampPalette(rev(brewer.pal(7, "RdBu")))(100),
                breaks = seq(-1, 1, length.out = 101),
                clustering_distance_rows = "correlation",
                clustering_distance_cols = "correlation",
                main = paste(gene, "Module - Consensus Correlation"),
                fontsize_row = 10, fontsize_col = 10)
  pdf(file.path(out_dir, paste0(gene, "_merged_heatmaps.pdf")), width = 7, height = 6)
  print(p); dev.off()
  png(file.path(out_dir, paste0(gene, "_merged_heatmaps.png")), width = 7, height = 6, units = "in", res = 300)
  print(p); dev.off()
  # Bar plot
  if (gene %in% rownames(mat)) {
    mean_rs <- sapply(sub_mats, function(m) mean(m[gene, names(m[gene, ]) != gene], na.rm = TRUE))
    overall_mean <- mean(mean_rs)
    sd_r <- sd(mean_rs)
    
    summary_data <- data.frame(Dataset = datasets, Mean_r = mean_rs, Module = gene)
    write.csv(summary_data, file.path(out_dir, paste0(gene, "_summary_barplot_data.csv")), row.names = FALSE)
    
    p_bar <- ggplot(summary_data, aes(x = Dataset, y = Mean_r, fill = Mean_r > 0)) +
      geom_bar(stat = "identity", width = 0.6) +
      geom_errorbar(aes(ymin = Mean_r - sd_r/3, ymax = Mean_r + sd_r/3), width = 0.2) +
      scale_fill_manual(values = c("FALSE" = "red", "TRUE" = "blue"), guide = "none") +
      labs(title = paste("Consensus Correlation of", gene, "with Module Genes"),
           x = "Dataset", y = "Average r") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(out_dir, paste0(gene, "_summary_barplot.png")), p_bar, width = 6, height = 4, dpi = 300)
    
    cat("Module", gene, "completed! Average r =", round(overall_mean, 4), "±", round(sd_r, 4), "\n")
  }
}
cat("All processing completed!\nInput from: outputs/05_scRNA_Analysis/dataset/gene/\nOutput saved to: outputs/05_scRNA_Analysis/Consensus_Correlation/gene/\n")