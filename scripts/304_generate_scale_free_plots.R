# =============================================
# WGCNA Scale-free Topology Plot 
# Final optimized version - compatible with soft_threshold_results.csv format
# =============================================

library(ggplot2)
library(gridExtra)
library(grid)

# ========================================================
# Standard Path Management Template (All Modules)
# ========================================================
# install.packages("rprojroot") # Only execute once for the first time
library(rprojroot)
# Automatically detect project root directory (works on any computer, any drive)
project_root <- find_rstudio_root_file()
# Fixed paths (shared by all modules)
data_dir <- file.path(project_root, "data")
outputs_dir <- file.path(project_root, "outputs")
# ------------------- Only modify this line -------------------
module_name <- "03_WGCNA" # ←←← Modify here!!
# ------------------------------------------------------
# Module actual directory (automatically created)
module_dir <- file.path(outputs_dir, module_name)
dir.create(module_dir, showWarnings = FALSE, recursive = TRUE)

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

datasets <- c("GSE55235", "GSE55457", "GSE55584")

for (gse in datasets) {
  
  csv_path <- file.path(module_dir, gse, "Results_WGCNA_GO-KEGG-Reactome", "soft_threshold_results.csv")
  
  if (!file.exists(csv_path)) {
    cat("❌ File not found:", csv_path, "\n")
    next
  }
  
  # Key fix: skip text header
  lines <- readLines(csv_path, warn = FALSE)
  start_idx <- grep("^\\s*Power", lines)[1]
  
  if (is.na(start_idx)) {
    cat("⚠️ Could not find Power header in", gse, "\n")
    next
  }
  
  df <- read.table(text = lines[start_idx:length(lines)], 
                   header = TRUE, 
                   sep = "", 
                   fill = TRUE, 
                   stringsAsFactors = FALSE,
                   check.names = FALSE)
  
  cat("✓ Successfully read", gse, "- ", nrow(df), "rows of data\n")
  
  # Column name standardization
  colnames(df) <- trimws(colnames(df))
  if ("mean.k." %in% colnames(df)) colnames(df)[colnames(df) == "mean.k."] <- "mean.k"
  
  # Plot 1: Scale-free Topology Model Fit
  p1 <- ggplot(df, aes(x = Power, y = SFT.R.sq)) +
    geom_point(color = "#D73027", size = 3.8) +
    geom_line(color = "#D73027", linewidth = 1.2) +
    geom_hline(yintercept = 0.85, linetype = "dashed", color = "blue", linewidth = 1) +
    annotate("text", x = max(df$Power)*0.65, y = 0.92, label = "R² ≥ 0.85", 
             color = "darkblue", size = 4.8, fontface = "bold") +
    labs(title = paste("Scale-free Topology Model Fit -", gse),
         x = "Soft Threshold Power", 
         y = expression(SFT.R^2)) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Plot 2: Mean Connectivity
  p2 <- ggplot(df, aes(x = Power, y = mean.k)) +
    geom_point(color = "#4575B4", size = 3.8) +
    geom_line(color = "#4575B4", linewidth = 1.2) +
    labs(title = paste("Mean Connectivity -", gse),
         x = "Soft Threshold Power", 
         y = "Mean Connectivity") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Combine panels
  combined <- grid.arrange(p1, p2, ncol = 2, 
                           top = textGrob(paste("WGCNA Soft-threshold Power Selection:", gse),
                                          gp = gpar(fontsize = 16, fontface = "bold")))
  
  # [MODIFIED] Triple-format output via save_fig
  fig_stem <- file.path(module_dir, gse, "Results_WGCNA_GO-KEGG-Reactome",
                        paste0(gse, "_scale_free_topology_S2"))
  save_fig(combined, fig_stem, w_in = 14, h_in = 6)
  
  cat("✅ Successfully generated:", basename(fig_stem), "(.pdf/_300.png/_600.tiff)\n")
}

cat("\n=== All scale-free topology figures have been successfully generated! ===\n")
cat("Files are saved in each GSE's Results_WGCNA_GO-KEGG-Reactome folder.\n")
