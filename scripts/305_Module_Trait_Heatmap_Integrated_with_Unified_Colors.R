# ====================== 005_Module_trait_heatmap_integrated.R ======================
# Standard Path Management Template (All Modules)
library(rprojroot)
project_root <- find_rstudio_root_file()
data_dir     <- file.path(project_root, "data")
outputs_dir  <- file.path(project_root, "outputs")

module_name <- "03_WGCNA/Module_trait_heatmap"
module_dir  <- file.path(outputs_dir, module_name)
dir.create(module_dir, showWarnings = FALSE, recursive = TRUE)

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

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

# -------------------------- 1. Define file paths --------------------------
wgcna_base_path <- file.path(outputs_dir, "03_WGCNA")
deg_intersection_dir <- file.path(outputs_dir, "02_Lysosomal_DEG_enrichment_intersection")

datasets <- c("GSE55235", "GSE55457", "GSE55584")

file_paths <- list(
  module_trait = setNames(
    file.path(wgcna_base_path, datasets, "Results_WGCNA_GO-KEGG-Reactome", "module_trait_heatmap_data.csv"),
    datasets
  ),
  lysosome_dist = setNames(
    file.path(wgcna_base_path, datasets, "Results_WGCNA_GO-KEGG-Reactome", "lysosome_module_dist.csv"),
    datasets
  ),
  module_genes = setNames(
    file.path(wgcna_base_path, datasets, "Results_WGCNA_GO-KEGG-Reactome", "moduleColors.csv"),
    datasets
  )
)

# -------------------------- 2. Read lysosomal DEGs --------------------------
deg_file <- file.path(deg_intersection_dir, "Lysosome_DEGs_Intersection_C33.csv")
deg_data <- read.csv(deg_file, stringsAsFactors = FALSE)
lys_deg_genes <- unique(deg_data$Gene.Symbol)

# -------------------------- 3. Get lysosomal modules --------------------------
get_lysosomal_modules <- function(file_path, min_genes = 10) {
  lines <- readLines(file_path)
  module_lines <- lines[grep("^\\$", lines)]
  modules <- sub("^\\$(.*)", "\\1", module_lines)
  gene_counts <- as.numeric(gsub(".*\\[1\\] (\\d+).*", "\\1", lines[grep("^\\[1\\]", lines)]))
  lysosome_data <- data.frame(Module = modules, LysosomeGeneCount = gene_counts)
  lysosomal_modules <- lysosome_data %>%
    filter(LysosomeGeneCount >= min_genes) %>%
    pull(Module)
  return(paste0("ME", lysosomal_modules))
}

lysosomal_modules_list <- lapply(names(file_paths$lysosome_dist), function(dataset) {
  modules <- get_lysosomal_modules(file_paths$lysosome_dist[[dataset]])
  data.frame(Dataset = dataset, Module = modules)
}) %>% bind_rows()

# -------------------------- 4. Read module-trait data --------------------------
combined_data <- lapply(names(file_paths$module_trait), function(dataset) {
  df <- read.csv(file_paths$module_trait[[dataset]], stringsAsFactors = FALSE)
  df$Dataset <- dataset
  return(df)
}) %>% bind_rows()

# Filter significant lysosomal modules
integrated_data <- combined_data %>%
  inner_join(lysosomal_modules_list, by = c("Dataset", "Module")) %>%
  filter(P_value < 0.05)

if (nrow(integrated_data) == 0) {
  stop("No significant lysosomal modules found.")
}

# -------------------------- 5. Module gene lists & DEG sets --------------------------
module_genes_list <- lapply(names(file_paths$module_genes), function(dataset) {
  df <- read.csv(file_paths$module_genes[[dataset]], stringsAsFactors = FALSE)
  df$Dataset <- dataset
  df$Module <- paste0("ME", df$ModuleColor)
  df <- df %>% filter(Gene %in% lys_deg_genes)
  return(df)
}) %>% bind_rows()

module_deg_genes <- module_genes_list %>%
  group_by(Dataset, Module) %>%
  summarise(DEG_Genes = list(unique(Gene)), .groups = "drop")

# -------------------------- 6. Unified module mapping (no threshold) --------------------------
reference_dataset <- "GSE55235"
ref_modules <- module_deg_genes %>% filter(Dataset == reference_dataset)

jaccard_similarity <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  if (union == 0) return(0)
  return(intersection / union)
}

unified_module_map <- module_deg_genes %>%
  filter(Dataset != reference_dataset) %>%
  rowwise() %>%
  mutate(
    Best_Match = ref_modules$Module[which.max(sapply(ref_modules$DEG_Genes, function(ref_genes) jaccard_similarity(DEG_Genes, ref_genes)))]
  ) %>%
  ungroup() %>%
  select(Dataset, Original_Module = Module, Unified_Module = Best_Match)

ref_map <- ref_modules %>%
  select(Dataset, Original_Module = Module, Unified_Module = Module)

all_module_map <- bind_rows(ref_map, unified_module_map)

integrated_data_unified <- integrated_data %>%
  left_join(all_module_map, by = c("Dataset", "Module" = "Original_Module")) %>%
  mutate(Module = coalesce(Unified_Module, Module)) %>%
  select(-Unified_Module)

# -------------------------- 6.1 Summarize to handle duplicates --------------------------
integrated_data_unified <- integrated_data_unified %>%
  group_by(Dataset, Module) %>%
  summarise(Correlation = mean(Correlation), P_value = min(P_value), .groups = "drop")

# -------------------------- 7. Save CSV data --------------------------
output_dir <- module_dir
csv_path <- file.path(output_dir, "module_trait_heatmap_integrated_unified_data.csv")
write.csv(integrated_data_unified, csv_path, row.names = FALSE)
cat("CSV saved to:", csv_path, "\n")

# -------------------------- 8. Prepare heatmap data --------------------------
heatmap_data <- integrated_data_unified %>%
  mutate(
    Dataset_Label = paste(Dataset, Module),
    Fill_Color = if_else(Correlation > 0, "Positive", "Negative"),
    P_value_Label = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01 ~ "**",
      P_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  arrange(desc(Correlation))

# -------------------------- 9. Plot --------------------------
p <- ggplot(heatmap_data, aes(x = "RA Trait", y = Dataset_Label)) +
  geom_tile(aes(fill = Fill_Color), color = "white", size = 1) +
  geom_text(aes(label = P_value_Label), color = "white", size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Positive" = "#DC143C", "Negative" = "#4169E1")) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(color = "black", face = "plain"),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold", lineheight = 1.1)
  ) +
  labs(
    title = "Lysosomal Module - RA Correlation\n(Cross-Dataset Unified Modules)",
    x = "RA Trait", y = "Dataset",
    fill = "Correlation Sign\nSignificance:\n*** P < 0.001\n** P < 0.01\n* P < 0.05"
  )

# -------------------------- 10. Save figure --------------------------
# [MODIFIED] Triple-format output via save_fig + keep SVG
fig_stem <- file.path(output_dir, "module_trait_heatmap_final")
save_fig(p, fig_stem, w_in = 9, h_in = 8)
svg_path <- file.path(output_dir, "module_trait_heatmap_final.svg")
ggsave(svg_path, p, width = 9, height = 8)
cat("Heatmap saved to:", fig_stem, "(.pdf/_300.png/_600.tiff)\n")
cat("SVG saved to:", svg_path, "\n")
# =====================================================================
