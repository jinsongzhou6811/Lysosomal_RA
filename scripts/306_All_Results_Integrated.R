# Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(igraph)
library(VennDiagram)
library(metap)
# ========================================================
# Standard Path Management Template (Newly Added, Portable Path)
# ========================================================
library(rprojroot) # Automatically find the root directory where .Rproj is located
project_root <- find_rstudio_root_file()
outputs_dir <- file.path(project_root, "outputs")
wgcna_dir <- file.path(outputs_dir, "03_WGCNA")
# This module output directory (Integrated)
base_dir <- file.path(wgcna_dir, "Integrated")
dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)

# ========================================================
# [STANDARDISED] Triple-format figure export helpers
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

save_base_fig <- function(plot_expr, stem, w_in, h_in) {
  pdf(paste0(stem, ".pdf"), width = w_in, height = h_in)
  eval(plot_expr); dev.off()
  png(paste0(stem, "_300.png"), width = w_in, height = h_in, units = "in", res = 300)
  eval(plot_expr); dev.off()
  tiff(paste0(stem, "_600.tiff"), width = w_in, height = h_in, units = "in", res = 600, compression = "lzw")
  eval(plot_expr); dev.off()
}

# Original result paths for the three GSEs (automatically adapt to any computer)
gse_paths <- list(
  GSE55235 = file.path(wgcna_dir, "GSE55235/Results_WGCNA_GO-KEGG-Reactome/Lysosome_Related_Files/Genes"),
  GSE55457 = file.path(wgcna_dir, "GSE55457/Results_WGCNA_GO-KEGG-Reactome/Lysosome_Related_Files/Genes"),
  GSE55584 = file.path(wgcna_dir, "GSE55584/Results_WGCNA_GO-KEGG-Reactome/Lysosome_Related_Files/Genes")
)
gse_ids <- c("GSE55235", "GSE55457", "GSE55584")
# Get gene list, excluding HLA-DQA1
genes <- c("GZMB", "LAMP3", "MREG", "NKG7", "PRF1", "SLC39A8", "TRAF3IP3", "VOPP1")
cat("Gene list for analysis:\n")
print(genes)
# ====================== Get module names (Code 1 filename rules) ======================
get_edge_modules <- function(gene, gse_id, gse_path) {
  input_dir <- file.path(gse_path, gene)
  edge_files <- list.files(input_dir, pattern = paste0("Cytoscape_edges_.*_", gene, ".txt"))
  if (length(edge_files) == 0) return(NULL)
  sub(paste0("Cytoscape_edges_(.*)_", gene, ".txt"), "\\1", edge_files)
}
# Main loop: process each gene (all code below is exactly the same as original, no changes made)
for (gene in genes) {
  cat(sprintf("\n--- Processing gene: %s ---\n", gene))
 
  # Set gene directory and output directory
  gene_dir <- file.path(base_dir, gene)
  output_dir <- file.path(gene_dir, paste0(gene, "_Integrated_Analysis"))
 
  if (!dir.exists(gene_dir)) {
    dir.create(gene_dir, recursive = TRUE)
    cat(sprintf("Created gene directory: %s\n", gene_dir))
  }
 
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
 
  # Initialize log file
  log_file <- file.path(output_dir, paste0(gene, "_analysis_log.txt"))
  sink(log_file, append = TRUE)
  cat(sprintf("Starting analysis for gene: %s Time: %s\n", gene, Sys.time()))
  sink()
 
  setwd(gene_dir)
 
  # 1. Load gene statistics data (all below remains unchanged)
  df_gene_stats <- list()
  for (gse_id in gse_ids) { # Note: gse_ids variable is defined below, here use original code order
    file_name <- file.path(gse_paths[[gse_id]], gene, paste0(gene, "_", gse_id, "__repaired.csv"))
    if (file.exists(file_name)) {
      df <- read.csv(file_name, stringsAsFactors = FALSE)
      if (!all(c("Key", "Value") %in% colnames(df))) {
        warning_msg <- sprintf("Warning: File %s does not contain Key and Value columns, skipping\n", file_name)
        cat(warning_msg)
        sink(log_file, append = TRUE)
        cat(warning_msg)
        sink()
        next
      }
      df_wide <- df %>%
        pivot_wider(names_from = Key, values_from = Value, values_fn = list(Value = first)) %>%
        mutate(数据集 = gse_id)
      df_gene_stats[[gse_id]] <- df_wide
    } else {
      warning_msg <- sprintf("Warning: File does not exist: %s\n", file_name)
      cat(warning_msg)
      sink(log_file, append = TRUE)
      cat(warning_msg)
      sink()
    }
  }
 
  # Merge gene statistics data (all below is identical to Code 2)
  if (length(df_gene_stats) > 0) {
    df_gene <- bind_rows(df_gene_stats)
    cat(sprintf("Merged %s gene statistics data:\n", gene))
    print(df_gene)
   
    write.csv(df_gene, file.path(output_dir, paste0(gene, "_Integrated_Gene_Stats.csv")), row.names = FALSE)
    cat(sprintf("%s merged gene statistics data has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_Integrated_Gene_Stats.csv"))))
   
    pvalue_cols <- c("DE_Pvalue_RA_OA", "DE_Pvalue_RA_ND", "GS_Pvalue_RA_status", "GS_Pvalue_RA_ND", "GS_Pvalue_RA_OA")
    existing_pvalue_cols <- intersect(pvalue_cols, colnames(df_gene))
    if (length(existing_pvalue_cols) > 0) {
      for (col in existing_pvalue_cols) {
        cat(sprintf("%s column exists, data type:\n", col))
        print(str(df_gene[[col]]))
        df_gene[[col]] <- as.numeric(df_gene[[col]])
        if (any(is.na(df_gene[[col]]))) {
          warning_msg <- sprintf("Warning: Column %s contains NA values, possibly due to 'NA' or non-numeric content in original data\n", col)
          cat(warning_msg)
          sink(log_file, append = TRUE)
          cat(warning_msg)
          sink()
        }
      }
     
      df_gene_long <- df_gene %>%
        select(数据集, all_of(existing_pvalue_cols)) %>%
        pivot_longer(cols = all_of(existing_pvalue_cols), names_to = "Pvalue_Type", values_to = "Pvalue") %>%
        filter(!is.na(Pvalue))
     
      if (nrow(df_gene_long) > 0) {
        p_pvalue_box <- ggplot(df_gene_long, aes(x = 数据集, y = Pvalue, fill = Pvalue_Type)) +
          geom_boxplot() +
          theme_minimal() +
          labs(title = paste(gene, "P-value Distribution"), x = "Dataset", y = "P-value") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        # [MODIFIED] Triple-format output via save_fig
        save_fig(p_pvalue_box, file.path(output_dir, paste0(gene, "_Pvalue_Boxplot")), w_in = 8, h_in = 6)
        cat(sprintf("%s P-value distribution boxplot has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_Pvalue_Boxplot.png"))))
      } else {
        warning_msg <- sprintf("Warning: %s has no valid p-value data, cannot generate boxplot\n", gene)
        cat(warning_msg)
        sink(log_file, append = TRUE)
        cat(warning_msg)
        sink()
      }
    } else {
      warning_msg <- sprintf("Warning: %s gene statistics data lacks all p-value columns, skipping p-value analysis\n", gene)
      cat(warning_msg)
      sink(log_file, append = TRUE)
      cat(warning_msg)
      sink()
    }
   
    status_cols <- c("Hub_Gene", "Core_Gene", "Significant_Gene", "Pathway_Gene")
    existing_status_cols <- intersect(status_cols, colnames(df_gene))
    if (length(existing_status_cols) > 0) {
      df_gene_status <- df_gene %>%
        select(数据集, all_of(existing_status_cols)) %>%
        pivot_longer(cols = all_of(existing_status_cols), names_to = "Gene_Type", values_to = "Status")
     
      p_status_bar <- ggplot(df_gene_status, aes(x = 数据集, fill = Status)) +
        geom_bar(position = "dodge") +
        facet_wrap(~Gene_Type) +
        theme_minimal() +
        labs(title = paste0(gene, " Gene Status Distribution"), x = "Dataset", y = "Count") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      # [MODIFIED] Triple-format output via save_fig
      save_fig(p_status_bar, file.path(output_dir, paste0(gene, "_Gene_Status_Bar")), w_in = 10, h_in = 6)
      cat(sprintf("%s gene status bar chart has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_Gene_Status_Bar.png"))))
    }
   
    if ("Enrichment" %in% colnames(df_gene)) {
      enrichment_data <- df_gene %>%
        select(数据集, Enrichment) %>%
        filter(!is.na(Enrichment) & Enrichment != "")
     
      if (nrow(enrichment_data) > 0) {
        enrichment_split <- enrichment_data %>%
          separate_rows(Enrichment, sep = ",") %>%
          mutate(Enrichment = trimws(Enrichment))
       
        enrichment_count <- enrichment_split %>%
          group_by(Enrichment) %>%
          summarise(Count = n(), .groups = 'drop')
       
        write.csv(enrichment_count, file.path(output_dir, paste0(gene, "_Enrichment_P.csv")), row.names = FALSE)
        cat(sprintf("%s Enrichment P-value statistics table has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_Enrichment_P.csv"))))
       
        p_enrichment <- ggplot(enrichment_count, aes(x = reorder(Enrichment, -Count), y = Count)) +
          geom_bar(stat = "identity", fill = "purple") +
          theme_minimal() +
          labs(title = paste0(gene, " Enrichment Pathway Distribution"), x = "Pathway", y = "Count") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        # [MODIFIED] Triple-format output via save_fig
        save_fig(p_enrichment, file.path(output_dir, paste0(gene, "_Enrichment_Bar")), w_in = 8, h_in = 6)
        cat(sprintf("%s Enrichment pathway bar chart has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_Enrichment_Bar.png"))))
      }
    }
  }
 
  # 2. Load GSEA GO data (using Code 1 naming rules)
  df_go_list <- list()
  for (gse_id in gse_ids) {
    file_name <- file.path(gse_paths[[gse_id]], gene, "gsea_GO.csv")
    if (file.exists(file_name)) {
      df <- read.csv(file_name, stringsAsFactors = FALSE)
      if (nrow(df) > 0) {
        df$数据集 <- gse_id
        df_go_list[[gse_id]] <- df
      } else {
        warning_msg <- sprintf("Warning: GO file %s is empty, skipping\n", file_name)
        cat(warning_msg)
        sink(log_file, append = TRUE)
        cat(warning_msg)
        sink()
      }
    } else {
      warning_msg <- sprintf("Warning: GO file does not exist: %s\n", file_name)
      cat(warning_msg)
      sink(log_file, append = TRUE)
      cat(warning_msg)
      sink()
    }
  }
  # All GO parts below are identical to Code 2 (common_go, NES pivot, combined_p, heatmap, empty files etc.)
  if (length(df_go_list) > 0) {
    df_go <- bind_rows(df_go_list)
    cat(sprintf("Merged %s GO data:\n", gene))
    print(head(df_go))
   
    common_go <- df_go %>%
      group_by(ID) %>%
      summarise(Count = n()) %>%
      filter(Count >= 2)
    cat(sprintf("%s Common GO pathways table:\n", gene))
    print(common_go %>% left_join(df_go %>% select(ID, Description) %>% unique(), by = "ID") %>% select(ID, Description, Count))
   
    df_common_go <- df_go %>%
      filter(ID %in% common_go$ID)
    cat(sprintf("%s Common GO pathway data example:\n", gene))
    print(head(df_common_go))
   
    if (nrow(df_common_go) == 0) {
      warning_msg <- sprintf("Warning: %s has no common GO data, generating empty NES and P-value files\n", gene)
      cat(warning_msg)
      sink(log_file, append = TRUE)
      cat(warning_msg)
      sink()
      write.csv(data.frame(), file.path(output_dir, paste0(gene, "_GO_NES_Pivot.csv")), row.names = FALSE)
      cat(sprintf("%s empty GO NES table has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_GO_NES_Pivot.csv"))))
      write.csv(data.frame(), file.path(output_dir, paste0(gene, "_GO_Combined_P.csv")), row.names = FALSE)
      cat(sprintf("%s empty GO combined P-value table has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_GO_Combined_P.csv"))))
    } else if (!"pvalue" %in% colnames(df_common_go)) {
      warning_msg <- sprintf("Warning: %s df_common_go is missing 'pvalue' column\n", gene)
      cat(warning_msg)
      sink(log_file, append = TRUE)
      cat(warning_msg)
      sink()
    } else {
      cat(sprintf("%s pvalue column exists, data type:\n", gene))
      print(str(df_common_go$pvalue))
      if (!is.numeric(df_common_go$pvalue) || any(is.na(df_common_go$pvalue))) {
        warning_msg <- sprintf("Warning: %s pvalue column contains non-numeric or NA values\n", gene)
        cat(warning_msg)
        sink(log_file, append = TRUE)
        cat(warning_msg)
        sink()
      } else {
        nes_go_pivot <- df_common_go %>%
          select(ID, NES, 数据集) %>%
          pivot_wider(names_from = 数据集, values_from = NES) %>%
          mutate(mean_NES = rowMeans(select(., all_of(intersect(colnames(.), gse_ids))), na.rm = TRUE)) %>%
          arrange(desc(mean_NES))
        cat(sprintf("%s GO NES Pivot table:\n", gene))
        print(nes_go_pivot)
       
        write.csv(nes_go_pivot, file.path(output_dir, paste0(gene, "_GO_NES_Pivot.csv")), row.names = FALSE)
        cat(sprintf("%s GO NES Pivot table has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_GO_NES_Pivot.csv"))))
       
        combined_p_go <- df_common_go %>%
          group_by(ID) %>%
          summarise(combined_p = pchisq(-2 * sum(log(.data$pvalue)), df = 2 * n(), lower.tail = FALSE), .groups = 'drop')
        cat(sprintf("%s GO combined P-value table:\n", gene))
        print(combined_p_go)
       
        write.csv(combined_p_go, file.path(output_dir, paste0(gene, "_GO_Combined_P.csv")), row.names = FALSE)
        cat(sprintf("%s GO combined P-value table has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_GO_Combined_P.csv"))))
        # ==================== Original heatmap part → commented out ====================
        # nes_matrix_go <- as.matrix(nes_go_pivot[, intersect(colnames(nes_go_pivot), gse_ids), drop = FALSE])
        # rownames(nes_matrix_go) <- nes_go_pivot$ID
        # png(file.path(output_dir, paste0(gene, "_GO_NES_Heatmap.png")))
        # heatmap(nes_matrix_go, Colv = NA, Rowv = NA, col = colorRampPalette(c("blue", "white", "red"))(100),
        # main = paste0(gene, " GO NES Heatmap"), xlab = "Dataset", ylab = "GO ID")
        # dev.off()
        # cat(sprintf("%s GO NES heatmap has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_GO_NES_Heatmap.png"))))
      }
    }
  } else {
    warning_msg <- sprintf("Warning: %s GO files are all missing, generating empty NES and P-value files\n", gene)
    cat(warning_msg)
    sink(log_file, append = TRUE)
    cat(warning_msg)
    sink()
    write.csv(data.frame(), file.path(output_dir, paste0(gene, "_GO_NES_Pivot.csv")), row.names = FALSE)
    cat(sprintf("%s empty GO NES table has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_GO_NES_Pivot.csv"))))
    write.csv(data.frame(), file.path(output_dir, paste0(gene, "_GO_Combined_P.csv")), row.names = FALSE)
    cat(sprintf("%s empty GO combined P-value table has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_GO_Combined_P.csv"))))
  }
 
  # 3. Load GSEA Reactome data (using Code 1 naming rules)
  df_react_list <- list()
  for (gse_id in gse_ids) {
    file_name <- file.path(gse_paths[[gse_id]], gene, "gsea_Reactome.csv")
    if (file.exists(file_name)) {
      df <- read.csv(file_name, stringsAsFactors = FALSE)
      if (nrow(df) > 0) {
        df$数据集 <- gse_id
        df_react_list[[gse_id]] <- df
      } else {
        warning_msg <- sprintf("Warning: Reactome file %s is empty, skipping\n", file_name)
        cat(warning_msg)
        sink(log_file, append = TRUE)
        cat(warning_msg)
        sink()
      }
    } else {
      warning_msg <- sprintf("Warning: Reactome file does not exist: %s\n", file_name)
      cat(warning_msg)
      sink(log_file, append = TRUE)
      cat(warning_msg)
      sink()
    }
  }
  # All Reactome parts below are identical to Code 2 (omitted repetition here, copy corresponding part from Code 2 when running)
  if (length(df_react_list) > 0) {
    df_react <- bind_rows(df_react_list)
    cat(sprintf("Merged %s Reactome data:\n", gene))
    print(head(df_react))
   
    common_react <- df_react %>%
      group_by(ID) %>%
      summarise(Count = n()) %>%
      filter(Count >= 2)
    cat(sprintf("%s Common Reactome pathways table:\n", gene))
    print(common_react %>% left_join(df_react %>% select(ID, Description) %>% unique(), by = "ID") %>% select(ID, Description, Count))
   
    df_common_react <- df_react %>%
      filter(ID %in% common_react$ID)
    cat(sprintf("%s Common Reactome pathway data example:\n", gene))
    print(head(df_common_react))
   
    if (!"pvalue" %in% colnames(df_common_react)) {
      warning_msg <- sprintf("Warning: %s df_common_react is missing 'pvalue' column\n", gene)
      cat(warning_msg)
      sink(log_file, append = TRUE)
      cat(warning_msg)
      sink()
    } else {
      cat(sprintf("%s Reactome pvalue column exists, data type:\n", gene))
      print(str(df_common_react$pvalue))
      if (!is.numeric(df_common_react$pvalue) || any(is.na(df_common_react$pvalue))) {
        warning_msg <- sprintf("Warning: %s pvalue column contains non-numeric or NA values\n", gene)
        cat(warning_msg)
        sink(log_file, append = TRUE)
        cat(warning_msg)
        sink()
      } else {
        available_gse_ids <- unique(df_common_react$数据集)
        nes_react_pivot <- df_common_react %>%
          select(ID, NES, 数据集) %>%
          pivot_wider(names_from = 数据集, values_from = NES) %>%
          mutate(mean_NES = rowMeans(select(., all_of(intersect(colnames(.), available_gse_ids))), na.rm = TRUE)) %>%
          arrange(desc(mean_NES))
        cat(sprintf("%s Reactome NES Pivot table:\n", gene))
        print(nes_react_pivot)
       
        write.csv(nes_react_pivot, file.path(output_dir, paste0(gene, "_Reactome_NES_Pivot.csv")), row.names = FALSE)
        cat(sprintf("%s Reactome NES Pivot table has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_Reactome_NES_Pivot.csv"))))
       
        combined_p_react <- df_common_react %>%
          group_by(ID) %>%
          summarise(combined_p = pchisq(-2 * sum(log(.data$pvalue)), df = 2 * n(), lower.tail = FALSE), .groups = 'drop')
        cat(sprintf("%s Reactome combined P-value table:\n", gene))
        print(combined_p_react)
       
        write.csv(combined_p_react, file.path(output_dir, paste0(gene, "_Reactome_Combined_P.csv")), row.names = FALSE)
        cat(sprintf("%s Reactome combined P-value table has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_Reactome_Combined_P.csv"))))
        # ==================== Original heatmap part → commented out ====================
        # nes_matrix_react <- as.matrix(nes_react_pivot[, intersect(colnames(nes_react_pivot), available_gse_ids), drop = FALSE])
        # rownames(nes_matrix_react) <- nes_react_pivot$ID
        # png(file.path(output_dir, paste0(gene, "_Reactome_NES_Heatmap.png")))
        # heatmap(nes_matrix_react, Colv = NA, Rowv = NA, col = colorRampPalette(c("blue", "white", "red"))(100),
        # main = paste0(gene, " Reactome NES Heatmap"), xlab = "Dataset", ylab = "Reactome ID")
        # dev.off()
        # cat(sprintf("%s Reactome NES heatmap has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_Reactome_NES_Heatmap.png"))))
      }
    }
  } else {
    warning_msg <- sprintf("Warning: %s Reactome files are all missing, generating empty NES and P-value files\n", gene)
    cat(warning_msg)
    sink(log_file, append = TRUE)
    cat(warning_msg)
    sink()
    write.csv(data.frame(), file.path(output_dir, paste0(gene, "_Reactome_NES_Pivot.csv")), row.names = FALSE)
    cat(sprintf("%s empty Reactome NES table has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_Reactome_NES_Pivot.csv"))))
    write.csv(data.frame(), file.path(output_dir, paste0(gene, "_Reactome_Combined_P.csv")), row.names = FALSE)
    cat(sprintf("%s empty Reactome combined P-value table has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_Reactome_Combined_P.csv"))))
  }
 
  # 4. Load GSEA KEGG data (using Code 1 automatic matching rules)
  df_kegg_list <- list()
  for (gse_id in gse_ids) {
    kegg_dir <- file.path(gse_paths[[gse_id]], gene)
    kegg_files <- list.files(kegg_dir, pattern = "^gsea_KEGG.*\\.csv$", full.names = TRUE)
    if (length(kegg_files) > 0) {
      file_name <- kegg_files[1]
      df <- read.csv(file_name, stringsAsFactors = FALSE)
      if (nrow(df) > 0) {
        df$数据集 <- gse_id
        df_kegg_list[[gse_id]] <- df
      } else {
        warning_msg <- sprintf("Warning: KEGG file %s is empty, skipping\n", file_name)
        cat(warning_msg)
        sink(log_file, append = TRUE)
        cat(warning_msg)
        sink()
      }
    } else {
      warning_msg <- sprintf("Warning: KEGG file does not exist: %s\n", file_name)
      cat(warning_msg)
      sink(log_file, append = TRUE)
      cat(warning_msg)
      sink()
    }
  }
  # All KEGG parts below are identical to Code 2
  if (length(df_kegg_list) > 0) {
    df_kegg <- bind_rows(df_kegg_list)
    cat(sprintf("Merged %s KEGG data:\n", gene))
    print(head(df_kegg))
   
    if (!"pval" %in% colnames(df_kegg)) {
      warning_msg <- sprintf("Warning: %s df_kegg is missing 'pval' column\n", gene)
      cat(warning_msg)
      sink(log_file, append = TRUE)
      cat(warning_msg)
      sink()
    } else {
      cat(sprintf("%s KEGG pval column exists, data type:\n", gene))
      print(str(df_kegg$pval))
      if (!is.numeric(df_kegg$pval) || any(is.na(df_kegg$pval))) {
        warning_msg <- sprintf("Warning: %s pval column contains non-numeric or NA values\n", gene)
        cat(warning_msg)
        sink(log_file, append = TRUE)
        cat(warning_msg)
        sink()
      } else {
        available_gse_ids <- unique(df_kegg$数据集)
        nes_kegg_pivot <- df_kegg %>%
          select(pathway, NES, 数据集) %>%
          pivot_wider(names_from = 数据集, values_from = NES) %>%
          mutate(mean_NES = rowMeans(select(., all_of(intersect(colnames(.), available_gse_ids))), na.rm = TRUE))
        cat(sprintf("%s KEGG NES Pivot table:\n", gene))
        print(nes_kegg_pivot)
       
        write.csv(nes_kegg_pivot, file.path(output_dir, paste0(gene, "_KEGG_NES_Pivot.csv")), row.names = FALSE)
        cat(sprintf("%s KEGG NES table has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_KEGG_NES_Pivot.csv"))))
       
        combined_p_kegg <- df_kegg %>%
          group_by(pathway) %>%
          summarise(combined_p = pchisq(-2 * sum(log(.data$pval)), df = 2 * n(), lower.tail = FALSE), .groups = 'drop')
        cat(sprintf("%s KEGG combined P-value table:\n", gene))
        print(combined_p_kegg)
       
        write.csv(combined_p_kegg, file.path(output_dir, paste0(gene, "_KEGG_Combined_P.csv")), row.names = FALSE)
        cat(sprintf("%s KEGG combined P-value table has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_KEGG_Combined_P.csv"))))
      }
    }
  } else {
    warning_msg <- sprintf("Warning: %s KEGG files do not exist, generating empty files\n", gene)
    cat(warning_msg)
    sink(log_file, append = TRUE)
    cat(warning_msg)
    sink()
    write.csv(data.frame(), file.path(output_dir, paste0(gene, "_KEGG_NES_Pivot.csv")), row.names = FALSE)
    cat(sprintf("%s empty KEGG NES file has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_KEGG_NES_Pivot.csv"))))
    write.csv(data.frame(), file.path(output_dir, paste0(gene, "_KEGG_Combined_P.csv")), row.names = FALSE)
    cat(sprintf("%s empty KEGG combined P-value file has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_KEGG_Combined_P.csv"))))
  }
 
  # 5. Network edge analysis (using Code 1 naming rules + get_edge_modules)
  df_edges_list <- list()
  for (gse_id in gse_ids) {
    edge_modules <- get_edge_modules(gene, gse_id, gse_paths[[gse_id]])
    if (is.null(edge_modules)) {
      warning_msg <- sprintf("Warning: Gene %s has no edge files in %s\n", gene, gse_id)
      cat(warning_msg)
      sink(log_file, append = TRUE)
      cat(warning_msg)
      sink()
      next
    }
   
    for (edge_module in edge_modules) {
      file_name <- file.path(gse_paths[[gse_id]], gene, paste0("Cytoscape_edges_", edge_module, "_", gene, ".txt"))
      if (file.exists(file_name)) {
        df <- read.table(file_name, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        required_cols <- c("fromNode", "toNode", "weight", "direction")
        if (!all(required_cols %in% colnames(df))) {
          warning_msg <- sprintf("Warning: Edge file %s is missing required columns (%s), skipping\n", file_name, paste(required_cols, collapse = ", "))
          cat(warning_msg)
          sink(log_file, append = TRUE)
          cat(warning_msg)
          sink()
          next
        }
        if (nrow(df) > 0) {
          df$module <- edge_module
          df_edges_list[[paste(gse_id, edge_module, sep = "_")]] <- df
          cat(sprintf("Loaded edge file: %s, contains %d edges\n", file_name, nrow(df)))
        } else {
          warning_msg <- sprintf("Warning: Edge file %s is empty, skipping\n", file_name)
          cat(warning_msg)
          sink(log_file, append = TRUE)
          cat(warning_msg)
          sink()
        }
      } else {
        warning_msg <- sprintf("Warning: Edge file does not exist: %s\n", file_name)
        cat(warning_msg)
        sink(log_file, append = TRUE)
        cat(warning_msg)
        sink()
      }
    }
  }
  # All network analysis parts below are identical to Code 2 (merged_edges, nodes, Venn, network plots, empty files etc.)
  if (length(df_edges_list) > 0) {
    df_merged_edges <- bind_rows(df_edges_list)
    cat(sprintf("Initial merged %s edge data:\n", gene))
    print(head(df_merged_edges))
   
    df_merged_edges <- df_merged_edges %>%
      mutate(node1 = pmin(fromNode, toNode),
             node2 = pmax(fromNode, toNode))
    cat(sprintf("Standardized %s merged edge data example:\n", gene))
    print(head(df_merged_edges))
   
    df_unique_edges <- df_merged_edges %>%
      group_by(node1, node2) %>%
      summarise(weight = max(weight, na.rm = TRUE),
                direction = first(direction),
                module = paste(unique(module), collapse = ","),
                .groups = 'drop') %>%
      rename(fromNode = node1, toNode = node2)
    cat(sprintf("%s unique edge processed data:\n", gene))
    print(head(df_unique_edges))
   
    edge_file_csv <- file.path(output_dir, paste0(gene, "_Cytoscape_edges.csv"))
    write.csv(df_unique_edges, edge_file_csv, row.names = FALSE, quote = FALSE)
    cat(sprintf("%s merged edge file (CSV) has been saved as %s\n", gene, edge_file_csv))
   
    edge_file_txt <- file.path(output_dir, paste0(gene, "_Cytoscape_edges.txt"))
    write.table(df_unique_edges, edge_file_txt, sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("%s merged edge file (TXT) has been saved as %s\n", gene, edge_file_txt))
   
    all_nodes <- unique(c(df_unique_edges$fromNode, df_unique_edges$toNode))
    cat(sprintf("%s all unique nodes:\n", gene))
    print(all_nodes)
   
    g_merged <- graph_from_data_frame(df_unique_edges, directed = FALSE)
    node_degrees <- degree(g_merged)
    cat(sprintf("%s node degree example:\n", gene))
    print(head(node_degrees))
   
    node_modules <- df_merged_edges %>%
      pivot_longer(cols = c(fromNode, toNode), values_to = "node") %>%
      group_by(node) %>%
      summarise(modules = paste(unique(module), collapse = ","), .groups = 'drop')
    cat(sprintf("%s node module example:\n", gene))
    print(head(node_modules))
   
    df_nodes <- data.frame(node = all_nodes) %>%
      left_join(node_modules, by = "node") %>%
      mutate(degree = node_degrees[node])
    cat(sprintf("%s node data example:\n", gene))
    print(head(df_nodes))
   
    node_file_csv <- file.path(output_dir, paste0(gene, "_Cytoscape_nodes.csv"))
    write.csv(df_nodes, node_file_csv, row.names = FALSE, quote = FALSE)
    cat(sprintf("%s node file (CSV) has been saved as %s\n", gene, node_file_csv))
   
    node_file_txt <- file.path(output_dir, paste0(gene, "_Cytoscape_nodes.txt"))
    write.table(df_nodes, node_file_txt, sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("%s node file (TXT) has been saved as %s\n", gene, node_file_txt))
   
    partners <- lapply(df_edges_list, function(df) {
      unique(c(df$fromNode, df$toNode)[c(df$fromNode, df$toNode) != gene])
    })
   
    for (gse_id in names(partners)) {
      cat(sprintf("%s interaction partners:\n", gse_id))
      print(head(partners[[gse_id]]))
      cat(sprintf("%s interaction partner count: %d\n", gse_id, length(partners[[gse_id]])))
    }
   
    if (length(partners) >= 2) {
      common_partners_at_least2 <- unique(c(
        intersect(partners[[1]], partners[[2]]),
        intersect(partners[[1]], partners[[3]]),
        intersect(partners[[2]], partners[[3]])
      ))
      cat(sprintf("%s partners common to at least 2 datasets:\n", gene))
      print(common_partners_at_least2)
     
      venn_file <- file.path(output_dir, paste0(gene, "_edges_venn.png"))
      venn.diagram(x=partners, filename=venn_file, main=paste0(gene, " Interaction Partners Venn Diagram"))
      cat(sprintf("%s Venn diagram has been saved as %s\n", gene, venn_file))
    }
   
    if (length(df_edges_list) > 0) {
      last_gse_id <- names(df_edges_list)[length(df_edges_list)]
      g_last <- graph_from_data_frame(df_edges_list[[last_gse_id]], directed = FALSE)
      # [MODIFIED] Triple-format output for network plot
      save_base_fig(quote({
        plot(g_last, vertex.size = 5, vertex.label.cex = 0.8, main = paste0(gene, " ", last_gse_id, " Network"))
      }), file.path(output_dir, paste0(gene, "_", last_gse_id, "_network")), w_in = 8, h_in = 8)
      cat(sprintf("%s %s network plot has been saved as %s\n", gene, last_gse_id, file.path(output_dir, paste0(gene, "_", last_gse_id, "_network.png"))))
    }
  } else {
    warning_msg <- sprintf("Warning: %s all edge files do not exist or are empty, skipping network analysis\n", gene)
    cat(warning_msg)
    sink(log_file, append = TRUE)
    cat(warning_msg)
    sink()
   
    node_file_csv <- file.path(output_dir, paste0(gene, "_Cytoscape_nodes.csv"))
    write.csv(data.frame(node = character(), modules = character(), degree = numeric()),
              node_file_csv, row.names = FALSE, quote = FALSE)
    cat(sprintf("%s empty node file (CSV) has been saved as %s\n", gene, node_file_csv))
   
    node_file_txt <- file.path(output_dir, paste0(gene, "_Cytoscape_nodes.txt"))
    write.table(data.frame(node = character(), modules = character(), degree = numeric()),
                node_file_txt, sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("%s empty node file (TXT) has been saved as %s\n", gene, node_file_txt))
  }
 
  # 6. Pathway classification integration (using Code 1 naming rules)
  df_path_list <- list()
  for (gse_id in gse_ids) {
    file_name <- file.path(gse_paths[[gse_id]], gene, paste0(gene, "_pathway_classification.csv"))
    if (file.exists(file_name)) {
      df <- read.csv(file_name, stringsAsFactors = FALSE)
      if (nrow(df) > 0) {
        df$数据集 <- gse_id
        df_path_list[[gse_id]] <- df
      } else {
        warning_msg <- sprintf("Warning: Pathway classification file %s is empty, skipping\n", file_name)
        cat(warning_msg)
        sink(log_file, append = TRUE)
        cat(warning_msg)
        sink()
      }
    } else {
      warning_msg <- sprintf("Warning: Pathway classification file does not exist: %s\n", file_name)
      cat(warning_msg)
      sink(log_file, append = TRUE)
      cat(warning_msg)
      sink()
    }
  }
 
  if (length(df_path_list) > 0) {
    df_path <- bind_rows(df_path_list)
    cat(sprintf("Merged %s pathway classification data:\n", gene))
    print(head(df_path))
   
    common_path <- df_path %>%
      group_by(ID) %>%
      summarise(Count = n()) %>%
      filter(Count == 3)
    cat(sprintf("%s Common pathways table:\n", gene))
    print(common_path %>% left_join(df_path %>% select(ID, Term, Category) %>% unique(), by = "ID") %>%
            select(ID, Term, Category, Count))
   
    write.csv(common_path %>% left_join(df_path %>% select(ID, Term, Category) %>% unique(), by = "ID"),
              file.path(output_dir, paste0(gene, "_Common_Pathways.csv")), row.names = FALSE)
    cat(sprintf("%s common pathways table has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_Common_Pathways.csv"))))
   
    category_count <- df_path %>% group_by(Category) %>% summarise(Count = n_distinct(ID))
    cat(sprintf("%s pathway category distribution:\n", gene))
    print(category_count)
   
    p_category <- ggplot(category_count, aes(x = "", y = Count, fill = Category)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y") +
      theme_void() +
      labs(title = paste0(gene, " Pathway Category Distribution"))
    # [MODIFIED] Triple-format output via save_fig
    save_fig(p_category, file.path(output_dir, paste0(gene, "_pathway_category_pie")), w_in = 8, h_in = 6)
    cat(sprintf("%s pathway category pie chart has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_pathway_category_pie.png"))))
  } else {
    warning_msg <- sprintf("Warning: %s pathway classification files are all missing, skipping pathway analysis\n", gene)
    cat(warning_msg)
    sink(log_file, append = TRUE)
    cat(warning_msg)
    sink()
    write.csv(data.frame(), file.path(output_dir, paste0(gene, "_Common_Pathways.csv")), row.names = FALSE)
    cat(sprintf("%s empty common pathways table has been saved as %s\n", gene, file.path(output_dir, paste0(gene, "_Common_Pathways.csv"))))
  }
 
  # 7. Advanced NES bar charts (all identical to Code 2)
  if (exists("nes_go_pivot") && nrow(nes_go_pivot) > 0 && !all(is.na(nes_go_pivot$mean_NES))) {
    p_go <- ggplot(nes_go_pivot, aes(x = reorder(ID, -mean_NES), y = mean_NES)) +
      geom_bar(stat = "identity", fill = "skyblue") +
      theme_minimal() +
      labs(title = paste0(gene, " All GO Average NES"), x = "GO ID", y = "Average NES") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    # [MODIFIED] Triple-format output via save_fig
    save_fig(p_go, file.path(output_dir, paste0(gene, "_all_go_mean_nes_bar")), w_in = 8, h_in = 6)
  } else {
    warning_msg <- sprintf("Warning: %s nes_go_pivot data is missing or invalid, skipping GO NES bar chart\n", gene)
    cat(warning_msg)
    sink(log_file, append = TRUE)
    cat(warning_msg)
    sink()
  }
 
  if (exists("nes_react_pivot") && nrow(nes_react_pivot) > 0 && !all(is.na(nes_react_pivot$mean_NES))) {
    p_react <- ggplot(nes_react_pivot, aes(x = reorder(ID, -mean_NES), y = mean_NES)) +
      geom_bar(stat = "identity", fill = "lightgreen") +
      theme_minimal() +
      labs(title = paste0(gene, " All Reactome Average NES"), x = "Reactome ID", y = "Average NES") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    # [MODIFIED] Triple-format output via save_fig
    save_fig(p_react, file.path(output_dir, paste0(gene, "_all_reactome_mean_nes_bar")), w_in = 8, h_in = 6)
  } else {
    warning_msg <- sprintf("Warning: %s nes_react_pivot data is missing or invalid, skipping Reactome NES bar chart\n", gene)
    cat(warning_msg)
    sink(log_file, append = TRUE)
    cat(warning_msg)
    sink()
  }
 
  if (exists("nes_kegg_pivot") && nrow(nes_kegg_pivot) > 0 && !all(is.na(nes_kegg_pivot$mean_NES))) {
    p_kegg <- ggplot(nes_kegg_pivot, aes(x = pathway, y = mean_NES)) +
      geom_bar(stat = "identity", fill = "orange") +
      theme_minimal() +
      labs(title = paste0(gene, " KEGG Average NES"), x = "Pathway", y = "Average NES") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    # [MODIFIED] Triple-format output via save_fig
    save_fig(p_kegg, file.path(output_dir, paste0(gene, "_kegg_mean_nes_bar")), w_in = 8, h_in = 6)
  } else {
    warning_msg <- sprintf("Warning: %s nes_kegg_pivot data is missing or invalid, skipping KEGG NES bar chart\n", gene)
    cat(warning_msg)
    sink(log_file, append = TRUE)
    cat(warning_msg)
    sink()
  }
 
  cat(sprintf("Completed analysis for gene %s\n", gene))
  sink(log_file, append = TRUE)
  cat(sprintf("Completed analysis for gene %s Time: %s\n", gene, Sys.time()))
  sink()
}
cat("All gene analyses completed. Output files are saved in their respective gene_Integrated_Analysis folders.\n")