# Lysosome_DEG_enrichment_intersection: Lysosomal protein gene screening, hypergeometric test and intersection analysis
# Function: Screen lysosomal-related DEGs, perform hypergeometric test and generate bar plots;
#           then calculate the intersection of three datasets and generate intersection bar plot

# Load required libraries
library(stats)
library(dplyr)
library(ggplot2)
library(readr)
library(purrr)
library(ggrepel)   
library(rprojroot) # Used to automatically locate the project root

# ========================================================
# Standard Path Management Template
# ========================================================
project_root <- find_rstudio_root_file()
data_dir <- file.path(project_root, "data")
outputs_dir <- file.path(project_root, "outputs")
module_name <- "02_Lysosomal_DEG_enrichment_intersection"
module_dir <- file.path(outputs_dir, module_name)
output_base_dir <- module_dir
dir.create(output_base_dir, showWarnings = FALSE, recursive = TRUE)

# Log file (keeps the same directory tree structure)
log_file <- file.path(output_base_dir, "Lysosome_DEG_enrichment_intersection_log.txt")
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

# Log function
log_message <- function(message, log_file, append = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- paste("[", timestamp, "] ", message, "\n", sep = "")
  cat(log_entry, file = log_file, append = append)
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

# Record module start
log_message("Merged Module 3 started", log_file, append = FALSE)
log_message(paste("Project root:", project_root), log_file)
log_message(paste("Data dir:", data_dir), log_file)
log_message(paste("Output dir:", output_base_dir), log_file)

# Define input/output paths
input_dir <- file.path(outputs_dir, "01_Differential_Expression_Analysis", "002_DEG_Analysis")
deg_input <- file.path(input_dir, "Summary_Results", "Common_DEGs_Full_Details.csv")
lysosome_genes_file <- file.path(data_dir, "lysosome_genes_metadata.csv")
summary_stats_file <- file.path(input_dir, "Summary_Results", "deg_summary.csv")
datasets <- c("GSE55235", "GSE55457", "GSE55584")
dataset_deg_files <- list(
  GSE55235 = file.path(input_dir, "GSE55235", "DEG_Results", "GSE55235_DEGs_genes.csv"),
  GSE55457 = file.path(input_dir, "GSE55457", "DEG_Results", "GSE55457_DEGs_genes.csv"),
  GSE55584 = file.path(input_dir, "GSE55584", "DEG_Results", "GSE55584_DEGs_genes.csv")
)

# Check input files
input_files <- c(deg_input, lysosome_genes_file, summary_stats_file, unlist(dataset_deg_files))
for (file in input_files) {
  if (!file.exists(file)) stop(sprintf("Input file not found: %s", file))
}
log_message("Input file existence check passed", log_file)

# Create output subfolders for each dataset
for (dataset in datasets) {
  dir.create(file.path(output_base_dir, dataset), recursive = TRUE, showWarnings = FALSE)
}
log_message("Output directories checked or created", log_file)

# Read input files
deg_list <- read.csv(deg_input, stringsAsFactors = FALSE, check.names = FALSE, quote = "\"", fill = TRUE, blank.lines.skip = TRUE)
if ("avg_FDR" %in% colnames(deg_list)) {
  colnames(deg_list)[colnames(deg_list) == "avg_FDR"] <- "adj.P.Val"
  log_message("Renamed 'avg_FDR' to 'adj.P.Val'", log_file)
}

tryCatch({
  lysosome_genes <- read.csv(lysosome_genes_file, stringsAsFactors = FALSE, check.names = FALSE, quote = "\"", fill = TRUE, blank.lines.skip = TRUE)
  log_message("Successfully read lysosome_genes_metadata.csv", log_file)
}, error = function(e) {
  log_message(sprintf("Failed to read lysosome_genes_metadata.csv: %s", conditionMessage(e)), log_file)
  tryCatch({
    lysosome_genes <- read.table(lysosome_genes_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = "\"", fill = TRUE, blank.lines.skip = TRUE)
    log_message("Successfully read lysosome_genes_metadata.csv with tab delimiter", log_file)
  }, error = function(e2) {
    log_message(sprintf("Tab-delimited read failed: %s", conditionMessage(e2)), log_file)
    stop(sprintf("Failed to read lysosome_genes_metadata.csv: %s", conditionMessage(e)))
  })
})

summary_stats <- read.csv(summary_stats_file, stringsAsFactors = FALSE, check.names = FALSE, quote = "\"", fill = TRUE, blank.lines.skip = TRUE)
log_message("Input files read completed", log_file)

# Debug information
log_message(sprintf("lysosome_genes column names: %s", paste(colnames(lysosome_genes), collapse=", ")), log_file)
log_message(sprintf("lysosome_genes first 10 rows:\n%s", paste(capture.output(print(head(lysosome_genes, 10))), collapse="\n")), log_file)

# Handle Gene Names column name variants
if ("Gene_Names" %in% colnames(lysosome_genes)) {
  colnames(lysosome_genes)[colnames(lysosome_genes) == "Gene_Names"] <- "Gene Names"
  log_message("Renamed 'Gene_Names' to 'Gene Names'", log_file)
} else if ("GeneNames" %in% colnames(lysosome_genes)) {
  colnames(lysosome_genes)[colnames(lysosome_genes) == "GeneNames"] <- "Gene Names"
  log_message("Renamed 'GeneNames' to 'Gene Names'", log_file)
}

# Validate required columns
required_cols_deg <- c("Gene.Symbol", "adj.P.Val")
if (!all(required_cols_deg %in% colnames(deg_list))) {
  stop(sprintf("deg_list is missing required columns: %s", paste(setdiff(required_cols_deg, colnames(deg_list)), collapse=", ")))
}
if (!all(c("Symbol", "Gene Names") %in% colnames(lysosome_genes))) {
  stop(sprintf("lysosome_genes is missing required columns: %s", paste(setdiff(c("Symbol", "Gene Names"), colnames(lysosome_genes)), collapse=", ")))
}
if (!"Total_Genes" %in% colnames(summary_stats)) {
  stop("summary_stats is missing the 'Total_Genes' column")
}
log_message("Data column validation passed", log_file)

deg_list$Gene.Symbol <- as.character(deg_list$Gene.Symbol)
log_message("Converted Gene.Symbol to character vector", log_file)
log_message(sprintf("deg_list structure: %d rows, %d columns", nrow(deg_list), ncol(deg_list)), log_file)
log_message(sprintf("Gene.Symbol type: %s", class(deg_list$Gene.Symbol)), log_file)

# Data preprocessing
deg_list$Gene.Symbol <- sapply(strsplit(deg_list$Gene.Symbol, " /// "), function(x) x[1])
log_message("Cleaned multiple gene symbols in Gene.Symbol", log_file)
deg_list <- deg_list[!is.na(deg_list$Gene.Symbol) & deg_list$Gene.Symbol != "Unknown", ]
log_message(sprintf("After filtering NA/Unknown, deg_list has %d rows remaining", nrow(deg_list)), log_file)
if (nrow(deg_list) == 0) stop("No valid genes left after filtering NA/Unknown")

# Preprocess lysosome_genes
lysosome_gene_names <- unique(c(
  lysosome_genes$Symbol,
  unlist(strsplit(lysosome_genes$`Gene Names`, " "))
))
lysosome_gene_names <- lysosome_gene_names[!is.na(lysosome_gene_names) & lysosome_gene_names != ""]
log_message(sprintf("Total lysosomal gene names (Symbol + Gene Names): %d", length(lysosome_gene_names)), log_file)

# Process lysosomal-related DEGs per dataset
overlap_summary_all <- data.frame(
  Dataset = character(),
  Overlap_Count = integer(),
  P_Value = numeric(),
  Expected_Overlap = numeric(),
  Enrichment_Fold = numeric(),
  stringsAsFactors = FALSE
)

total_genes <- mean(summary_stats$Total_Genes, na.rm = TRUE)
if (is.na(total_genes) || total_genes <= 0) stop(sprintf("Invalid total_genes: %s", total_genes))
lysosome_gene_count <- nrow(lysosome_genes)
log_message(sprintf("Background gene total: %.0f, lysosomal gene total: %d", total_genes, lysosome_gene_count), log_file)

for (dataset in datasets) {
  dataset_deg_file <- dataset_deg_files[[dataset]]
  dataset_deg <- read.csv(dataset_deg_file, stringsAsFactors = FALSE, check.names = FALSE, quote = "\"", fill = TRUE, blank.lines.skip = TRUE)
  
  required_cols_dataset <- c("Gene.Symbol", "adj.P.Val", "logFC")
  if (!all(required_cols_dataset %in% colnames(dataset_deg))) {
    stop(sprintf("%s is missing required columns: %s", dataset, paste(setdiff(required_cols_dataset, colnames(dataset_deg)), collapse=", ")))
  }
  
  dataset_deg$Dataset <- dataset
  dataset_deg$Direction <- ifelse(dataset_deg$logFC > 0, "Up", ifelse(dataset_deg$logFC < 0, "Down", "Neutral"))
  log_message(sprintf("%s: Added Dataset and Direction columns", dataset), log_file)
  
  dataset_deg$Gene.Symbol <- as.character(dataset_deg$Gene.Symbol)
  dataset_deg$Dataset <- as.character(dataset_deg$Dataset)
  dataset_deg$Gene.Symbol <- sapply(strsplit(dataset_deg$Gene.Symbol, " /// "), function(x) x[1])
  dataset_deg <- dataset_deg[!is.na(dataset_deg$Gene.Symbol) & dataset_deg$Gene.Symbol != "Unknown", ]
  
  if (nrow(dataset_deg) == 0) {
    log_message(sprintf("%s: No valid genes left after filtering NA/Unknown", dataset), log_file)
    next
  }
  
  lysosome_set <- dataset_deg[dataset_deg$Gene.Symbol %in% lysosome_gene_names, ]
  log_message(sprintf("%s: Found %d lysosomal-related DEGs", dataset, nrow(lysosome_set)), log_file)
  
  output_lysosome <- file.path(output_base_dir, dataset, sprintf("%s_lysosome_related_DEGs.rds", dataset))
  output_lysosome_csv <- file.path(output_base_dir, dataset, sprintf("%s_lysosome_related_DEGs.csv", dataset))
  saveRDS(lysosome_set, output_lysosome)
  write.csv(lysosome_set, output_lysosome_csv, row.names = FALSE)
  
  if (nrow(lysosome_set) == 0) {
    log_message(sprintf("%s: No lysosomal-related DEGs found, skipping hypergeometric test", dataset), log_file)
    next
  }
  
  lysosome_by_dataset <- lysosome_set %>%
    dplyr::group_by(Dataset) %>%
    dplyr::summarise(
      Lysosome_DEG_Count = n(),
      Up_Lysosome_DEG = sum(Direction == "Up"),
      Down_Lysosome_DEG = sum(Direction == "Down")
    ) %>%
    dplyr::ungroup()
  
  deg_count <- length(unique(dataset_deg$Gene.Symbol))
  overlap_count <- length(unique(dataset_deg$Gene.Symbol[dataset_deg$Gene.Symbol %in% lysosome_gene_names]))
  
  if (deg_count == 0) {
    p_value <- NA
    expected_overlap <- NA
    enrichment_fold <- NA
  } else {
    p_value <- phyper(
      q = overlap_count - 1,
      m = lysosome_gene_count,
      n = total_genes - lysosome_gene_count,
      k = deg_count,
      lower.tail = FALSE
    )
    expected_overlap <- (deg_count * lysosome_gene_count) / total_genes
    enrichment_fold <- ifelse(expected_overlap > 0, overlap_count / expected_overlap, NA)
  }
  
  overlap_summary_all <- rbind(overlap_summary_all, data.frame(
    Dataset = dataset,
    Overlap_Count = overlap_count,
    P_Value = p_value,
    Expected_Overlap = expected_overlap,
    Enrichment_Fold = enrichment_fold
  ))
  log_message(sprintf("%s: Hypergeometric test: overlap=%d, p=%.2e", dataset, overlap_count, p_value), log_file)
}

# Save comprehensive hypergeometric test results
output_summary <- file.path(output_base_dir, "overlap_summary.rds")
output_summary_csv <- file.path(output_base_dir, "overlap_summary.csv")
saveRDS(overlap_summary_all, output_summary)
write.csv(overlap_summary_all, output_summary_csv, row.names = FALSE)
log_message("Comprehensive hypergeometric test results saved", log_file)

# Visualization: Lysosome DEG overlap bar plot
if (nrow(overlap_summary_all) == 0) {
  stop("No overlap data available for visualization")
}
p <- ggplot(overlap_summary_all, aes(x = Dataset, y = Overlap_Count)) +
  geom_bar(stat = "identity", fill = "#4CAF50") +
  geom_text(aes(label = sprintf("p=%.2e", P_Value)), vjust = -0.5) +
  theme_minimal() +
  labs(title = "Lysosome-Related DEG Overlap by Dataset", x = "Dataset", y = "Overlap Count")
# [MODIFIED] Triple-format output via save_fig
save_fig(p, file.path(output_base_dir, "lysosome_deg_overlap"), w_in = 8, h_in = 6)
log_message("Comprehensive lysosomal DEG overlap bar chart saved (.pdf/_300.png/_600.tiff)", log_file)

# Generate per-dataset bar plots
for (dataset in datasets) {
  overlap_data <- overlap_summary_all %>% filter(Dataset == dataset)
  if (nrow(overlap_data) == 0) next
  p <- ggplot(overlap_data, aes(x = Dataset, y = Overlap_Count)) +
    geom_bar(stat = "identity", fill = "#4CAF50") +
    geom_text(aes(label = sprintf("p=%.2e", P_Value)), vjust = -0.5) +
    theme_minimal() +
    labs(title = sprintf("Lysosome-Related DEG Overlap for %s", dataset), x = "Dataset", y = "Overlap Count")
  # [MODIFIED] Triple-format output via save_fig
  save_fig(p, file.path(output_base_dir, dataset, "lysosome_deg_overlap"), w_in = 6, h_in = 4)
  log_message(sprintf("%s: Dataset-specific bar chart saved (.pdf/_300.png/_600.tiff)", dataset), log_file)
}

# --- Transition to intersection analysis ---
log_message("Starting intersection analysis section", log_file)

check_files <- function(ds) {
  file_path <- file.path(output_base_dir, ds, paste0(ds, "_lysosome_related_DEGs.csv"))
  if (!file.exists(file_path)) {
    stop(paste("File does not exist:", file_path))
  }
  return(file_path)
}

lysosome_degs <- tryCatch({
  map(datasets, function(ds) {
    file_path <- check_files(ds)
    log_message(sprintf("Reading: %s", file_path), log_file)
    read_csv(file_path, show_col_types = FALSE) %>%
      select(Gene.Symbol, logFC, adj.P.Val, Direction) %>%
      filter(!is.na(Gene.Symbol))
  }) %>% set_names(datasets)
}, error = function(e) {
  log_message(sprintf("Error: Unable to read files - %s", e$message), log_file)
  NULL
})

if (is.null(lysosome_degs)) {
  stop("File reading failed. Please check whether the generated lysosomal DEG files exist.")
}

combinations <- combn(datasets, 3, simplify = FALSE)
intersection_results <- list()
for (combo in combinations) {
  ds1 <- lysosome_degs[[combo[1]]]
  ds2 <- lysosome_degs[[combo[2]]]
  ds3 <- lysosome_degs[[combo[3]]]
  
  common_genes <- Reduce(intersect, list(ds1$Gene.Symbol, ds2$Gene.Symbol, ds3$Gene.Symbol))
  
  if (length(common_genes) > 0) {
    result <- bind_rows(
      ds1 %>% filter(Gene.Symbol %in% common_genes) %>% mutate(Dataset = combo[1]),
      ds2 %>% filter(Gene.Symbol %in% common_genes) %>% mutate(Dataset = combo[2]),
      ds3 %>% filter(Gene.Symbol %in% common_genes) %>% mutate(Dataset = combo[3])
    ) %>%
      mutate(Combination = paste(combo, collapse = "_vs_")) %>%
      select(Combination, Dataset, Gene.Symbol, logFC, adj.P.Val, Direction)
    
    intersection_results[[paste(combo, collapse = "_vs_")]] <- result
  }
}

final_results <- bind_rows(intersection_results)
if (nrow(final_results) > 0) {
  output_file <- file.path(output_base_dir, "Lysosome_DEGs_Intersection_C33.csv")
  write_csv(final_results, output_file)
  log_message(sprintf("Intersection analysis results saved to: %s", output_file), log_file)
} else {
  log_message("No lysosomal-related genes intersecting across the three datasets were found.", log_file)
}

if (nrow(final_results) > 0) {
  log_message("Intersection gene count statistics:", log_file)
  intersection_summary <- final_results %>%
    group_by(Combination, Gene.Symbol) %>%
    summarise(Gene_Count = n_distinct(Gene.Symbol), .groups = "drop") %>%
    group_by(Combination) %>%
    summarise(Total_Genes = n(), .groups = "drop")
  print(intersection_summary)
  
  p <- ggplot(intersection_summary, aes(x = Combination, y = Total_Genes)) +
    geom_bar(stat = "identity", fill = "#4CAF50") +
    geom_text(aes(label = Total_Genes), vjust = -0.5) +
    theme_minimal() +
    labs(title = "Lysosome DEGs Intersection Count", x = "Combination", y = "Total Intersection Genes")
  # [MODIFIED] Triple-format output via save_fig
  save_fig(p, file.path(output_base_dir, "lysosome_degs_intersection"), w_in = 8, h_in = 6)
  log_message("Intersection gene count bar chart saved (.pdf/_300.png/_600.tiff)", log_file)
}

log_message("Merged Module 3 completed", log_file)

# ==================================================================
# 4. Volcano Plot Generation for Lysosomal DEGs 
# ==================================================================
log_message("Starting Volcano Plot Generation for Lysosomal DEGs", log_file)

# Create a dedicated output directory for volcano maps
volcano_dir <- file.path(output_base_dir, "volcano_plots")
dir.create(volcano_dir, showWarnings = FALSE, recursive = TRUE)
log_message(paste("Volcano plots output directory:", volcano_dir), log_file)

# Define comparison name (improved title)
contrast_names <- c(
  "GSE55235" = "RA vs ND",
  "GSE55457" = "RA vs Control",
  "GSE55584" = "RA vs OA"
)

# Read the genes common to all three datasets (Common_DEGs) for highlighting）
intersection_file <- file.path(output_base_dir, "Lysosome_DEGs_Intersection_C33.csv")
common_genes <- character(0)
if (file.exists(intersection_file)) {
  intersection_data <- read.csv(intersection_file, stringsAsFactors = FALSE)
  common_genes <- unique(toupper(trimws(intersection_data$Gene.Symbol)))
  common_genes <- common_genes[!is.na(common_genes) & common_genes != ""]
  log_message(paste("Loaded common intersection genes:", length(common_genes)), log_file)
} else {
  log_message("WARNING: Intersection file not found. Common genes will not be highlighted.", log_file)
}

# Use the existing lysosome gene list (lysosome_gene_names has been prepared above)
lysosome_genes_upper <- toupper(lysosome_gene_names)

# Generate a volcano plot for each dataset
for (geo_id in datasets) {
  log_message(paste("Generating volcano plot for dataset:", geo_id), log_file)
  
  deg_file <- file.path(output_base_dir, geo_id, paste0(geo_id, "_lysosome_related_DEGs.csv"))
  
  if (!file.exists(deg_file)) {
    log_message(paste("WARNING: DEG file missing for", geo_id), log_file)
    next
  }
  
  deg_data <- read.csv(deg_file, stringsAsFactors = FALSE)
  
  # Data cleaning and standardization
  deg_data <- deg_data %>%
    filter(!is.na(Gene.Symbol), Gene.Symbol != "", !is.na(logFC), !is.na(adj.P.Val)) %>%
    mutate(
      Gene.Symbol = toupper(trimws(sapply(strsplit(as.character(Gene.Symbol), " /// "), function(x) x[1]))),
      Direction = ifelse(logFC > 0, "Up", "Down")
    ) %>%
    filter(Gene.Symbol != "")
  
  if (nrow(deg_data) == 0) {
    log_message(paste("WARNING: No valid data for", geo_id), log_file)
    next
  }
  
  # Marked groups (Common_DEG highlighted first)
  deg_data$group <- ifelse(deg_data$Gene.Symbol %in% common_genes, "Common_DEG",
                           ifelse(deg_data$Gene.Symbol %in% lysosome_genes_upper, "Lysosome", "Other"))
  
  # Keep only lysosomal-related genes
  volcano_data <- deg_data %>%
    filter(group %in% c("Lysosome", "Common_DEG")) %>%
    mutate(
      neg_log10_fdr = -log10(pmax(adj.P.Val, 1e-300)),
      Significance = case_when(
        adj.P.Val < 0.05 & abs(logFC) >= 1 & Direction == "Up"   ~ "Lysosome_Up",
        adj.P.Val < 0.05 & abs(logFC) >= 1 & Direction == "Down" ~ "Lysosome_Down",
        TRUE ~ "Not_Significant"
      ),
      label_text = ifelse(group == "Common_DEG",
                          paste0(Gene.Symbol, "\nlog2FC=", round(logFC, 2)),
                          NA_character_)
    )
  
  if (nrow(volcano_data) == 0) {
    log_message(paste("WARNING: No lysosomal genes available for plotting in", geo_id), log_file)
    next
  }
  
  # === Draw a volcano diagram ===
  p_volcano <- ggplot(volcano_data, aes(x = logFC, y = neg_log10_fdr, color = Significance)) +
    geom_point(size = 2.8, alpha = 0.85) +
    scale_color_manual(
      values = c("Lysosome_Up" = "#E31A1C", "Lysosome_Down" = "#1F78B4", "Not_Significant" = "grey65"),
      labels = c("Lysosome_Up" = "Lysosome Up", "Lysosome_Down" = "Lysosome Down", "Not_Significant" = "Not Significant")
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.6) +
    geom_text_repel(
      data = filter(volcano_data, group == "Common_DEG" & !is.na(label_text)),
      aes(label = label_text),
      size = 4.2, color = "black", fontface = "bold",
      box.padding = 1.3, point.padding = 0.9,
      max.overlaps = Inf, segment.color = "grey60", segment.size = 0.4
    ) +
    labs(
      title = paste("Lysosomal DEGs Volcano Plot -", geo_id),
      subtitle = contrast_names[geo_id],
      x = "Log2 Fold Change",
      y = "-Log10 (FDR)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold")
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  # [MODIFIED] Triple-format output via save_fig (replaces single PDF ggsave)
  tryCatch({
    save_fig(p_volcano, file.path(volcano_dir, paste0("volcano_lysosome_", geo_id)), w_in = 9.5, h_in = 7.5)
    log_message(paste("Volcano plot saved (.pdf/_300.png/_600.tiff):", geo_id), log_file)
  }, error = function(e) {
    log_message(paste("ERROR: Failed to save volcano for", geo_id, "-", e$message), log_file)
  })
}

log_message(paste("All lysosomal DEG volcano plots saved to:", volcano_dir), log_file)
log_message("Volcano plot generation completed successfully", log_file)
