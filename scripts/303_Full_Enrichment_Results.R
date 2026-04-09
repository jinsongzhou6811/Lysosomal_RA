# =====================================================================
# Create_CoreGene_Go_Kegg_Reactome_CSV.R
#
# Purpose: Generate data/Kegg_Reactome.csv from the 9 CSV files located in
# the data/raw_enrichment/ directory
#
# Filtering rules (strictly following project requirements):
# 1. Keep only pathways related to the 8 core lysosomal genes
# 2. For GO files, retain only BP and CC terms (completely remove MF)
# 3. When multiple rows share the same Description, keep only the one
# with the highest FoldEnrichment
# 4. Append type suffix to Description: (KEGG), (Reactome), (BP), (CC)
#
# How to use:
# 1. Ensure that 002_Copy_and_Rename_Raw_Enrichment.R has already been run
# (this generates the 9 files)
# 2. Simply run this script
# =====================================================================
library(rprojroot)
library(dplyr)
library(stringr)
select <- dplyr::select
# ====================== 1. Project Path Settings ======================
project_root <- find_rstudio_root_file()
raw_enrichment_dir <- file.path(project_root,"data", "Raw_enrichment")
data_dir <- file.path(project_root, "data")
output_csv <- file.path(data_dir, "Full_Enrichment_Results.csv")
# Create data directory (if it does not exist)
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
cat("=== Starting to generate Kegg_Reactome.csv ===\n")
cat("Input directory :", raw_enrichment_dir, "\n")
cat("Output file :", output_csv, "\n\n")
# ====================== 2. Core Lysosomal Gene List ======================
core_genes <- c("GZMB", "LAMP3", "MREG", "NKG7", "SLC2A6",
                "SLC39A8", "TRAF3IP3", "VOPP1")
# ====================== 3. Read and Process All 9 Files ======================
all_data <- list()
# --- 3.1 Process 3 KEGG files ---
kegg_files <- list.files(raw_enrichment_dir,
                         pattern = "_kegg_enrichment_core_lysosome\\.csv$",
                         full.names = TRUE)
for (f in kegg_files) {
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  df$Type <- "KEGG"
  df$Description <- paste0(trimws(df$Description), "(KEGG)")
  all_data[[length(all_data)+1]] <- df
  cat("✓ Processing KEGG :", basename(f), " →", nrow(df), " rows\n")
}
# --- 3.2 Process 3 Reactome files ---
reactome_files <- list.files(raw_enrichment_dir,
                             pattern = "_reactome_enrichment_core_lysosome\\.csv$",
                             full.names = TRUE)
for (f in reactome_files) {
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  df$Type <- "Reactome"
  df$Description <- paste0(trimws(df$Description), "(Reactome)")
  all_data[[length(all_data)+1]] <- df
  cat("✓ Processing Reactome :", basename(f), " →", nrow(df), " rows\n")
}
# --- 3.3 Process 3 GO files (keep only BP and CC) ---
go_files <- list.files(raw_enrichment_dir,
                       pattern = "_go_enrichment_core_lysosome\\.csv$",
                       full.names = TRUE)
for (f in go_files) {
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Keep only BP and CC, remove MF
  df <- df %>% filter(ONTOLOGY %in% c("BP", "CC"))
  
  df$Type <- df$ONTOLOGY # Type = BP or CC
  df$Description <- paste0(trimws(df$Description), "(", df$ONTOLOGY, ")")
  
  all_data[[length(all_data)+1]] <- df
  cat("✓ Processing GO :", basename(f), " → retained BP+CC:", nrow(df), " rows\n")
}
# ====================== 4. Combine All Data ======================
combined <- bind_rows(all_data)
cat("\nTotal records after merging (before filtering):", nrow(combined), "\n")
# ====================== 5. Core Gene Filtering ======================
# Keep rows where geneID contains any of the core genes
combined <- combined %>%
  filter(str_detect(geneID, paste(core_genes, collapse = "|")))
cat("Remaining records after core gene filtering:", nrow(combined), "\n")
# ====================== 6. Deduplication (keep the most significant one for the same pathway) ======================
final_df <- combined %>%
  group_by(Description) %>%
  slice_max(order_by = FoldEnrichment, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Type, desc(FoldEnrichment))
cat("Final record count after deduplication:", nrow(final_df), "\n")
# ====================== 7. Reorder Columns and Save ======================
final_columns <- c("Type", "Description", "GeneRatio", "BgRatio",
                   "RichFactor", "FoldEnrichment", "zScore",
                   "pvalue", "p.adjust", "qvalue", "geneID", "Count")
final_df <- final_df %>% select(all_of(intersect(final_columns, colnames(final_df))))
write.csv(final_df, output_csv, row.names = FALSE, quote = TRUE, na = "")
# ====================== 8. Result Summary ======================
cat("\n✓ Processing complete!\n")
cat("File saved to: ", output_csv, "\n")
cat("Final record count: ", nrow(final_df), " rows\n")
cat("Distribution by type:\n")
print(table(final_df$Type))
cat("\nPreview of top 10 rows:\n")
print(head(final_df[, c("Type", "Description", "FoldEnrichment", "p.adjust", "geneID")], 10))
cat("\nScript execution finished, ready for use in 004 Sankey visualization script!\n")