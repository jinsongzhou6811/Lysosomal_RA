# data/check_data_ready.R
# Final data readiness check script (updated 2026-03-05)
# Now includes the 14 new GWAS .vcf.gz files

library(rprojroot)
project_root <- find_rstudio_root_file()
data_dir <- file.path(project_root, "data")
gwas_dir <- file.path(data_dir, "GWAS")

cat("=== Lysosomal_RA Data Readiness Check (All Modules) ===\n\n")

# 1. Required small/STRINGdb files (in data/)
stringdb_files <- c(
  "9606.protein.aliases.v11.5.txt.gz", "9606.protein.aliases.v12.0.txt.gz",
  "9606.protein.info.v11.5.txt.gz",    "9606.protein.info.v12.0.txt.gz",
  "9606.protein.links.v11.5.txt.gz",   "9606.protein.links.v12.0.txt.gz"
)

# 2. Other large files (in data/)
other_large_files <- c(
  "GSE296117_RA_geo.rds",
  "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
  "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
  "hg38_screen_v10_clust.regions_vs_motifs.rankings.feather",
  "hg38_screen_v10_clust.regions_vs_motifs.scores.feather"
)

# 3. 14 GWAS VCF files (in data/GWAS/)
gwas_files <- c(
  "bbj-a-151.vcf.gz", "bbj-a-72.vcf.gz", "bbj-a-73.vcf.gz", "bbj-a-74.vcf.gz",
  "ieu-a-831.vcf.gz", "ieu-a-832.vcf.gz", "ieu-a-833.vcf.gz", "ieu-a-834.vcf.gz",
  "ukb-a-105.vcf.gz", "ukb-b-11874.vcf.gz", "ukb-b-9125.vcf.gz",
  "ukb-d-M06.vcf.gz", "ukb-d-M13_RHEUMA.vcf.gz", "ukb-d-RHEUMA_NOS.vcf.gz"
)

# 4. Required folders (after manual extraction of .tar files)
required_folders <- c("GSE159117_RAW", "GSE55235_RAW", "GSE55457_RAW", "GSE55584_RAW")

# Check everything
missing_stringdb <- stringdb_files[!file.exists(file.path(data_dir, stringdb_files))]
missing_other    <- other_large_files[!file.exists(file.path(data_dir, other_large_files))]
missing_gwas     <- gwas_files[!file.exists(file.path(gwas_dir, gwas_files))]
missing_folders  <- required_folders[!dir.exists(file.path(data_dir, required_folders))]

# Create GWAS folder if missing
if (!dir.exists(gwas_dir)) {
  dir.create(gwas_dir, recursive = TRUE)
  cat("✅ Created missing folder: data/GWAS/\n")
}

# Final report
if (length(missing_stringdb) == 0 && 
    length(missing_other) == 0 && 
    length(missing_gwas) == 0 && 
    length(missing_folders) == 0) {
  
  cat("🎉 ALL DATA IS READY!\n")
  cat("✅ 6 STRINGdb files\n")
  cat("✅ 5 other large files (feather + rds)\n")
  cat("✅ 14 GWAS VCF.gz files (in data/GWAS/)\n")
  cat("✅ 4 extracted GSE_RAW folders\n\n")
  cat("You can now run the pipeline scripts in scripts/ (start from 101_GEO Preprocessing.R)\n")
  
} else {
  cat("⚠️  Some files/folders are still missing!\n\n")
  
  if (length(missing_stringdb) > 0) {
    cat("Missing STRINGdb files (6):\n")
    cat(paste("   -", missing_stringdb, collapse = "\n"), "\n\n")
  }
  if (length(missing_other) > 0) {
    cat("Missing other large files (5):\n")
    cat(paste("   -", missing_other, collapse = "\n"), "\n\n")
  }
  if (length(missing_gwas) > 0) {
    cat("Missing GWAS files (14) in data/GWAS/:\n")
    cat(paste("   -", missing_gwas, collapse = "\n"), "\n\n")
  }
  if (length(missing_folders) > 0) {
    cat("Missing extracted folders (4):\n")
    cat(paste("   -", missing_folders, collapse = "\n"), "\n\n")
  }
  
  cat("Please open data/DOWNLOAD_INSTRUCTIONS.md for complete download and extraction steps.\n")
  cat("After fixing, run this script again.\n")
}

cat("\nData readiness check completed at", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")