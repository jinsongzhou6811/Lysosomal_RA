# ========================================================
# DEG Analysis Module (002)
# ========================================================
# install.packages("rprojroot")
library(rprojroot)
library(GEOquery)
library(limma)
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(dplyr)
library(tidyr)
library(grid)
library(viridis)
library(gridExtra)
project_root <- find_rstudio_root_file()
outputs_dir <- file.path(project_root, "outputs")
module_dir <- file.path(outputs_dir, "01_Differential_Expression_Analysis")
data_source_dir <- file.path(module_dir, "001_GEO Preprocessing", "results")
analysis_output_dir <- file.path(module_dir, "002_DEG_Analysis")
# Ensure output directory exists
if (!dir.exists(analysis_output_dir)) dir.create(analysis_output_dir, recursive = TRUE)
# Log file is placed in 01_Differential_Expression_Analysis
log_file_path <- file.path(module_dir, "geo_deg_analysis.log")
log_message <- function(type, message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- sprintf("[%s] %s: %s\n", timestamp, type, message)
  cat(log_entry)
  con <- file(log_file_path, open = "a")
  cat(log_entry, file = con)
  close(con)
}

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

# ------------------------------------------------------------------
# Parameters and Initialization
# ------------------------------------------------------------------
FDR_THRESHOLD <- 0.05
LOGFC_THRESHOLD <- 1
log_message("INFO", "DEG Analysis Module (002) started")
log_message("INFO", sprintf("Thresholds: FDR < %s, |logFC| >= %s", FDR_THRESHOLD, LOGFC_THRESHOLD))
datasets <- c("GSE55235", "GSE55457", "GSE55584")
deg_summary <- data.frame(
  Dataset = character(), Up = integer(), Down = integer(),
  Total_Probe_Level = integer(), Total_Gene_Level = integer(),
  Total_Genes = integer(), DEG_Proportion = numeric(),
  stringsAsFactors = FALSE
)
all_degs <- list()
gene_sets <- list()
# ------------------------------------------------------------------
# Loop through each dataset
# ------------------------------------------------------------------
for (dataset in datasets) {
  log_message("INFO", sprintf("Processing: %s", dataset))
  
  diff_expr_file <- file.path(data_source_dir, dataset, "Results_All",
                              paste0(dataset, "_diff_expr_normalized_single_gene.csv"))
  
  if (!file.exists(diff_expr_file)) {
    log_message("WARNING", sprintf("Input file not found: %s", diff_expr_file))
    next
  }
  
  diff_expr <- read.csv(diff_expr_file)
  
  if ("Gene Symbol" %in% colnames(diff_expr)) {
    diff_expr <- diff_expr %>% rename(Gene.Symbol = `Gene Symbol`)
  }
  
  diff_expr <- diff_expr[!is.na(diff_expr$logFC) & !is.na(diff_expr$adj.P.Val) &
                           !is.na(diff_expr$Gene.Symbol) & diff_expr$Gene.Symbol != "", ]
  diff_expr$adj.P.Val <- pmax(diff_expr$adj.P.Val, 1e-300)
  
  total_genes <- nrow(diff_expr)
  
  degs <- diff_expr[diff_expr$adj.P.Val < FDR_THRESHOLD & abs(diff_expr$logFC) >= LOGFC_THRESHOLD, ]
  
  # Output directory: outputs/01_Differential_Expression_Analysis/002_DEG_Analysis/{dataset}/DEG_Results
  deg_dir <- file.path(analysis_output_dir, dataset, "DEG_Results")
  if (!dir.exists(deg_dir)) dir.create(deg_dir, recursive = TRUE)
  
  write.csv(degs %>% arrange(adj.P.Val) %>% select(Gene.Symbol, logFC, adj.P.Val),
            file.path(deg_dir, paste0(dataset, "_DEGs_genes.csv")), row.names = FALSE)
  
  # Volcano plot (keep the original style from Module 001)
  volcano_data <- diff_expr %>%
    mutate(color = case_when(
      adj.P.Val < FDR_THRESHOLD & logFC >= LOGFC_THRESHOLD ~ "Up",
      adj.P.Val < FDR_THRESHOLD & logFC <= -LOGFC_THRESHOLD ~ "Down",
      TRUE ~ "Not Significant"
    ))
  
  volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = color)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_viridis_d(option = "viridis", direction = -1) +
    theme_minimal() +
    labs(title = paste("Volcano Plot -", dataset)) +
    geom_vline(xintercept = c(-LOGFC_THRESHOLD, LOGFC_THRESHOLD), linetype = "dashed", color = "gray") +
    geom_hline(yintercept = -log10(FDR_THRESHOLD), linetype = "dashed", color = "gray") +
    geom_text_repel(data = subset(volcano_data, adj.P.Val < 1e-10 & abs(logFC) > 2),
                    aes(label = Gene.Symbol), max.overlaps = 20, size = 3)
  
  # [MODIFIED] Triple-format output via save_fig
  save_fig(volcano_plot, file.path(deg_dir, paste0(dataset, "_volcano")), w_in = 8, h_in = 6)
  
  # Summary table
  deg_summary <- bind_rows(deg_summary, data.frame(
    Dataset = dataset, Up = sum(degs$logFC > 0), Down = sum(degs$logFC < 0),
    Total_Probe_Level = nrow(degs),
    Total_Gene_Level = length(unique(degs$Gene.Symbol)),
    Total_Genes = total_genes,
    DEG_Proportion = round(nrow(degs) / total_genes, 4)
  ))
  
  all_degs[[dataset]] <- degs
  gene_sets[[dataset]] <- unique(degs$Gene.Symbol)
}
# ------------------------------------------------------------------
# Final summary (placed in 002_DEG_Analysis/Summary_Results)
# ------------------------------------------------------------------
summary_dir <- file.path(analysis_output_dir, "Summary_Results")
if (!dir.exists(summary_dir)) dir.create(summary_dir, recursive = TRUE)
if (nrow(deg_summary) > 0) {
  write.csv(deg_summary, file.path(summary_dir, "deg_summary.csv"), row.names = FALSE)
  
  # Venn diagram + detailed common genes analysis
  if (length(gene_sets) >= 2) {
    
    venn_list <- venn.diagram(x = gene_sets, category.names = names(gene_sets),
                              filename = NULL, fill = viridis(length(gene_sets)), alpha = 0.6)
    venn_grob <- gTree(children = venn_list)
    
    common_genes <- Reduce(intersect, gene_sets)
    if (length(common_genes) > 0) {
      common_df <- data.frame(Gene.Symbol = common_genes)
      for (ds in names(all_degs)) {
        tmp <- all_degs[[ds]] %>% select(Gene.Symbol, logFC, adj.P.Val) %>%
          rename(!!paste0("logFC_", ds) := logFC, !!paste0("FDR_", ds) := adj.P.Val)
        common_df <- inner_join(common_df, tmp, by = "Gene.Symbol")
      }
      common_df$avg_logFC <- rowMeans(common_df[, grep("logFC_", names(common_df))], na.rm = TRUE)
      common_df$avg_FDR <- exp(rowMeans(log(common_df[, grep("FDR_", names(common_df))]), na.rm = TRUE))
      
      write.csv(common_df, file.path(summary_dir, "Common_DEGs_Full_Details.csv"), row.names = FALSE)
      
      top10 <- common_df %>% arrange(avg_FDR) %>% head(10) %>% select(Gene.Symbol)
      t_grob <- tableGrob(top10, rows = NULL, theme = ttheme_minimal(base_size = 8))
      
      # [MODIFIED] Triple-format output for Venn diagram (base R graphics)
      venn_stem <- file.path(summary_dir, "DEGs_Venn_Diagram")
      
      pdf(paste0(venn_stem, ".pdf"), width = 8, height = 10)
      grid.arrange(venn_grob, t_grob, ncol = 1, heights = c(3, 1))
      dev.off()
      
      png(paste0(venn_stem, "_300.png"), width = 8, height = 10, units = "in", res = 300)
      grid.arrange(venn_grob, t_grob, ncol = 1, heights = c(3, 1))
      dev.off()
      
      tiff(paste0(venn_stem, "_600.tiff"), width = 8, height = 10, units = "in", res = 600, compression = "lzw")
      grid.arrange(venn_grob, t_grob, ncol = 1, heights = c(3, 1))
      dev.off()
    }
  }
}
# Merge all DEGs
bind_rows(all_degs, .id = "Dataset") %>%
  write.csv(file.path(summary_dir, "All_Datasets_DEGs_Combined.csv"), row.names = FALSE)
log_message("INFO", "===== DEG Analysis Module (002) Completed Successfully =====")
