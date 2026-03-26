# =====================================================================
# 002_Copy_and_Rename_Raw_Enrichment.R
#
# Purpose: Automatically copy the 9 enrichment files from the three 
#          datasets under outputs/03_WGCNA, add dataset prefixes to 
#          avoid filename conflicts, and save them to the raw_enrichment/ 
#          folder at the project root.
#
# Usage:
#   1. Run 001_WGCNA.R first
#   2. Run this script directly
# =====================================================================

library(rprojroot)

# Automatically detect project root (works on any computer / drive)
project_root <- find_rstudio_root_file()

# Fixed paths (used across the whole project)
outputs_dir        <- file.path(project_root, "outputs")
raw_enrichment_dir <- file.path(project_root, "data", "Raw_enrichment")
dir.create(raw_enrichment_dir, showWarnings = FALSE, recursive = TRUE)

# ====================== Copy and Rename ======================
base_dir <- file.path(outputs_dir, "03_WGCNA")
datasets <- c("GSE55235", "GSE55457", "GSE55584")

file_list <- c(
  "kegg_enrichment_core_lysosome.csv",
  "reactome_enrichment_core_lysosome.csv",
  "go_enrichment_core_lysosome.csv"
)

cat("Starting to copy and rename raw enrichment files...\n\n")

for (ds in datasets) {
  src_path <- file.path(base_dir, ds, "Results_WGCNA_GO-KEGG-Reactome")
  
  for (f in file_list) {
    old_file <- file.path(src_path, f)
    if (file.exists(old_file)) {
      new_name <- paste0(ds, "_", f)
      new_file <- file.path(raw_enrichment_dir, new_name)
      file.copy(old_file, new_file, overwrite = TRUE)
      cat("  Copied:", basename(old_file), "→", new_name, "\n")
    }
  }
}

cat("\n✓ Done! 9 raw files have been copied to:\n")
cat("   ", raw_enrichment_dir, "\n")