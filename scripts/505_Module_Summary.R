# ========================================================
# 005_Module_Summary_Table3.R
# Input: outputs/05_scRNA_Analysis/{dataset}/{gene}/target_genes_summary.txt
# Output: outputs/05_scRNA_Analysis/module_summary.csv
# ========================================================
library(rprojroot)
library(dplyr)
project_root <- find_rstudio_root_file()
outputs_dir <- file.path(project_root, "outputs")
modules <- c("GZMB", "LAMP3", "NKG7", "TRAF3IP3")
datasets <- c("GSE159117", "GSE202375", "GSE296117")
# Output path (fixed under 05_scRNA_Analysis)
output_path <- file.path(outputs_dir, "05_scRNA_Analysis", "module_summary.csv")
dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
# ==================== Function ====================
parse_summary_txt <- function(file_path) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  df <- tryCatch({
    read.table(file_path, header = TRUE, sep = "", stringsAsFactors = FALSE,
               fill = TRUE, strip.white = TRUE, quote = "")
  }, error = function(e) {
    warning(paste("Error reading", file_path, ":", e$message))
    return(NULL)
  })
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  colnames(df) <- trimws(tolower(colnames(df)))
  col_map <- list(
    gene = c("gene"),
    mean = c("mean_expression", "meanexpression", "mean"),
    nonzero_prop = c("nonzero_prop", "nonzeroprop", "nonzero"),
    variance = c("variance", "var")
  )
  for (target in names(col_map)) {
    for (variant in col_map[[target]]) {
      if (variant %in% colnames(df)) {
        colnames(df)[colnames(df) == variant] <- target
        break
      }
    }
  }
  required_cols <- c("gene", "mean", "nonzero_prop", "variance")
  df <- df[, intersect(required_cols, colnames(df)), drop = FALSE]
  for (col in setdiff(required_cols, colnames(df))) df[[col]] <- NA
  return(df)
}
# ==================== Main loop ====================
table3 <- data.frame(Module = character(), Gene = character(),
                     Avg_Mean = numeric(), Avg_Nonzero_Prop = numeric(),
                     Avg_Variance = numeric(), stringsAsFactors = FALSE)
for (mod in modules) {
  stats_list <- list()
  
  for (ds in datasets) {
    txt_path <- file.path(outputs_dir, "05_scRNA_Analysis", ds, mod, "target_genes_summary.txt")
    stats_df <- parse_summary_txt(txt_path)
    if (!is.null(stats_df)) stats_list[[ds]] <- stats_df
  }
  
  if (length(stats_list) == 0) next
  
  all_genes <- unique(unlist(lapply(stats_list, `[[`, "gene")))
  avg_stats <- data.frame(Gene = all_genes, stringsAsFactors = FALSE)
  
  avg_stats$Avg_Mean <- sapply(all_genes, function(g) {
    means <- sapply(stats_list, function(df) if (g %in% df$gene) df$mean[df$gene == g] else NA)
    mean(means, na.rm = TRUE)
  })
  avg_stats$Avg_Nonzero_Prop <- sapply(all_genes, function(g) {
    props <- sapply(stats_list, function(df) if (g %in% df$gene) df$nonzero_prop[df$gene == g] else NA)
    mean(props, na.rm = TRUE)
  })
  avg_stats$Avg_Variance <- sapply(all_genes, function(g) {
    vars <- sapply(stats_list, function(df) if (g %in% df$gene) df$variance[df$gene == g] else NA)
    mean(vars, na.rm = TRUE)
  })
  
  avg_stats$Module <- mod
  table3 <- rbind(table3, avg_stats)
}
# ==================== Output ====================
write.csv(table3, output_path, row.names = FALSE)
cat("Table 3 has been generated!\n")
cat("Output path: ", normalizePath(output_path), "\n")
cat("Input path: outputs/05_scRNA_Analysis/{dataset}/{gene}/target_genes_summary.txt\n")