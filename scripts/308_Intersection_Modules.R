
library(rprojroot)
library(WGCNA)
library(dplyr)
library(pheatmap)
library(igraph)
library(RColorBrewer)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

project_root <- find_rstudio_root_file()
outputs_dir  <- file.path(project_root, "outputs")
wgcna_dir    <- file.path(outputs_dir, "03_WGCNA")   # Module 03 output root
out_dir      <- file.path(outputs_dir, "03_WGCNA", "Intersection_Modules")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(out_dir, "intersection_modules_log.txt")
writeLines(paste("Module 05 — Cross-Dataset Intersection\nStarted:", Sys.time()), log_path)
log_msg <- function(msg) {
  line <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_path, append = TRUE)
}

datasets      <- c("GSE55235", "GSE55457", "GSE55584")
R_THRESHOLD   <- 0.7    # eigengene-based module matching threshold (manuscript criterion)
KME_THRESHOLD <- 0.5    # minimum kME for gene inclusion in a CI module

# ==============================================================================
# STEP 1  Load expression matrix and module colors from Module 03 outputs
# ==============================================================================
# Module 03 Step 33 saves:
#   saveRDS(exprData, file.path(results_dir, "checkpoint_raw_data.rds"))
# exprData in Module 03 is samples × genes (rows = samples, cols = genes).
# Module 03 Step 8 saves:
#   write.csv(data.frame(Gene = colnames(exprData), ModuleColor = moduleColors), ...)

load_expr <- function(dataset) {
  path <- file.path(wgcna_dir, dataset,
                    "Results_WGCNA_GO-KEGG-Reactome",
                    "checkpoint_raw_data.rds")
  if (!file.exists(path))
    stop(paste(
      "checkpoint_raw_data.rds not found for", dataset, "\n",
      "Expected:", path, "\n",
      "Make sure Module 03 has completed successfully for this dataset."
    ))
  expr <- readRDS(path)
  # Module 03 saves exprData as samples × genes.
  # Defensive check (mirrors the preservation code in Module 03):
  # if somehow saved as genes × samples, transpose back.
  if (nrow(expr) > ncol(expr)) {
    log_msg(paste(dataset, ": expr appears to be genes×samples, transposing to samples×genes"))
    expr <- t(expr)
  }
  log_msg(paste(dataset, ": expr loaded —",
                nrow(expr), "samples ×", ncol(expr), "genes"))
  return(expr)   # returns: samples × genes matrix
}

load_module_colors <- function(dataset) {
  path <- file.path(wgcna_dir, dataset,
                    "Results_WGCNA_GO-KEGG-Reactome",
                    "moduleColors.csv")
  if (!file.exists(path))
    stop(paste("moduleColors.csv not found for", dataset, "\nExpected:", path))
  mc <- read.csv(path, stringsAsFactors = FALSE)
  # Module 03 writes: data.frame(Gene = colnames(exprData), ModuleColor = moduleColors)
  # so columns are exactly "Gene" and "ModuleColor"
  if (!all(c("Gene", "ModuleColor") %in% colnames(mc)))
    stop(paste("Unexpected columns in moduleColors.csv for", dataset,
               "— found:", paste(colnames(mc), collapse = ", ")))
  mc <- mc[mc$ModuleColor != "grey", ]   # exclude unassigned genes
  log_msg(paste(dataset, ": module colors loaded —",
                nrow(mc), "genes,", length(unique(mc$ModuleColor)), "non-grey modules"))
  return(mc)   # data.frame: Gene | ModuleColor
}

log_msg("=== STEP 1: Loading Module 03 outputs ===")
expr_list <- list()
mc_list   <- list()
for (ds in datasets) {
  expr_list[[ds]] <- load_expr(ds)
  mc_list[[ds]]   <- load_module_colors(ds)
}

# ==============================================================================
# STEP 2  Recompute module eigengenes (MEs) from saved expression + colors
# ==============================================================================
# moduleEigengenes(expr, colors) expects:
#   expr   : samples × genes  (matches Module 03 usage exactly)
#   colors : character vector of length = ncol(expr), names not required

compute_ME <- function(expr, mc_df, dataset) {
  # Align: keep only genes present in both expr and mc_df
  common_genes <- intersect(colnames(expr), mc_df$Gene)
  log_msg(paste(dataset, ": aligning genes — expr has", ncol(expr),
                "| mc_df has", nrow(mc_df),
                "| overlap:", length(common_genes)))
  if (length(common_genes) < 10)
    stop(paste(dataset, ": fewer than 10 genes overlap between expr and moduleColors"))
  
  expr_sub   <- expr[, common_genes]              # samples × common_genes
  colors_vec <- mc_df$ModuleColor[match(common_genes, mc_df$Gene)]  # aligned vector
  
  ME_result <- moduleEigengenes(expr_sub, colors = colors_vec, excludeGrey = TRUE)
  ME_mat    <- ME_result$eigengenes  # samples × modules, colnames = "ME{color}"
  log_msg(paste(dataset, ": eigengenes computed —",
                ncol(ME_mat), "modules, colnames:",
                paste(head(colnames(ME_mat), 5), collapse = ", ")))
  return(list(ME = ME_mat, expr_aligned = expr_sub, colors_aligned = colors_vec))
}

log_msg("=== STEP 2: Computing module eigengenes ===")
ME_data <- list()
for (ds in datasets) {
  ME_data[[ds]] <- compute_ME(expr_list[[ds]], mc_list[[ds]], ds)
}

# ==============================================================================
# STEP 3  Compute gene-module membership (kME) for cross-dataset comparison
# ==============================================================================
# kME[g, m] = cor(gene g expression across samples, module m eigengene across samples)
# Since both expr_aligned and ME are indexed by the same samples, this is simply:
#   cor(expr_aligned, ME)
# which R computes as: for each column pair (gene, module), Pearson r over rows (samples)
# Result shape: genes × modules  (colnames of expr × colnames of ME)

compute_kME <- function(expr_samples_x_genes, ME_samples_x_modules, dataset) {
  # Both inputs share the same row dimension (samples)
  stopifnot(nrow(expr_samples_x_genes) == nrow(ME_samples_x_modules))
  kME <- cor(expr_samples_x_genes, ME_samples_x_modules, use = "pairwise.complete.obs")
  # kME: rows = genes (colnames of expr), cols = modules (colnames of ME)
  log_msg(paste(dataset, ": kME computed —",
                nrow(kME), "genes ×", ncol(kME), "modules"))
  return(as.data.frame(kME))  # rownames = gene symbols, colnames = "ME{color}"
}

log_msg("=== STEP 3: Computing kME (gene-module membership) ===")
kME_list <- list()
for (ds in datasets) {
  kME_list[[ds]] <- compute_kME(
    ME_data[[ds]]$expr_aligned,
    ME_data[[ds]]$ME,
    ds
  )
  write.csv(kME_list[[ds]], file.path(out_dir, paste0(ds, "_kME.csv")))
}

# ==============================================================================
# STEP 4  Cross-dataset module-pair correlation via kME vectors
# ==============================================================================
# For each pair (module_i from ds_A, module_j from ds_B):
#   r = cor(kME_A[shared_genes, module_i], kME_B[shared_genes, module_j])
# This compares the gene-loading profiles of two modules across datasets.
# Shared genes = genes measured in both datasets (typically ~95% of probes
# on the same GPL570 platform, as used by all three GSE datasets here).

compute_cross_kME_cor <- function(kME_a, kME_b, ds_a, ds_b) {
  shared_genes <- intersect(rownames(kME_a), rownames(kME_b))
  log_msg(paste(" ", ds_a, "vs", ds_b, "— shared genes:", length(shared_genes)))
  if (length(shared_genes) < 20) {
    log_msg(paste("  WARNING: only", length(shared_genes),
                  "shared genes — correlation may be unreliable"))
  }
  sub_a <- as.matrix(kME_a[shared_genes, ])
  sub_b <- as.matrix(kME_b[shared_genes, ])
  # cor(sub_a, sub_b): rows = shared_genes (observations), cols = modules (variables)
  # Result: modules_a × modules_b matrix of Pearson r values
  cor_mat <- cor(sub_a, sub_b, use = "pairwise.complete.obs")
  return(cor_mat)
}

log_msg("=== STEP 4: Cross-dataset kME correlations ===")
ds_pairs <- list(
  c("GSE55235", "GSE55457"),
  c("GSE55235", "GSE55584"),
  c("GSE55457", "GSE55584")
)
cor_mats <- list()
for (p in ds_pairs) {
  key <- paste(p[1], p[2], sep = "_vs_")
  cor_mats[[key]] <- compute_cross_kME_cor(
    kME_list[[p[1]]], kME_list[[p[2]]], p[1], p[2])
  write.csv(cor_mats[[key]], file.path(out_dir, paste0("cormat_", key, ".csv")))
  log_msg(paste(" Saved:", key,
                "—", nrow(cor_mats[[key]]), "×", ncol(cor_mats[[key]]), "module pairs"))
}

# ==============================================================================
# STEP 5  Module matching at r > R_THRESHOLD
# ==============================================================================

find_matches <- function(cor_mat, ds_a, ds_b, threshold = R_THRESHOLD) {
  rows <- list()
  for (mod_a in rownames(cor_mat)) {
    for (mod_b in colnames(cor_mat)) {
      r <- cor_mat[mod_a, mod_b]
      if (!is.na(r) && abs(r) >= threshold) {
        rows[[length(rows) + 1]] <- data.frame(
          DS_A = ds_a, Module_A = mod_a,
          DS_B = ds_b, Module_B = mod_b,
          r    = round(r, 4),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

log_msg(paste("=== STEP 5: Finding module matches at |r| >", R_THRESHOLD, "==="))
match_list <- list()
for (key in names(cor_mats)) {
  parts <- strsplit(key, "_vs_")[[1]]
  mt    <- find_matches(cor_mats[[key]], parts[1], parts[2])
  if (nrow(mt) > 0) {
    match_list[[key]] <- mt
    log_msg(paste(" ", key, ":", nrow(mt), "matched pairs"))
    for (i in seq_len(min(8, nrow(mt))))
      log_msg(paste("   r =", mt$r[i], " | ", mt$Module_A[i], "<->", mt$Module_B[i]))
  } else {
    log_msg(paste(" ", key, ": 0 matched pairs at threshold", R_THRESHOLD))
  }
}

all_matches <- if (length(match_list) > 0) do.call(rbind, match_list) else data.frame()
write.csv(all_matches, file.path(out_dir, "module_matching_table.csv"), row.names = FALSE)

# If threshold gives no results, retry at 0.5 and report diagnostic
if (nrow(all_matches) == 0) {
  log_msg("WARNING: No matches at r > 0.7 — retrying at r > 0.5 for diagnostics")
  for (key in names(cor_mats)) {
    parts <- strsplit(key, "_vs_")[[1]]
    mt5   <- find_matches(cor_mats[[key]], parts[1], parts[2], threshold = 0.5)
    log_msg(paste(" ", key, "at r>0.5:", nrow(mt5), "pairs"))
    if (nrow(mt5) > 0)
      log_msg(paste("   Top pair: r =", mt5$r[1], mt5$Module_A[1], "<->", mt5$Module_B[1]))
  }
  log_msg("Check: (1) Was ComBat batch correction applied in Module 03?")
  log_msg("       (2) Are all three datasets on the same platform (GPL570)?")
  log_msg("       (3) Is the checkpoint_raw_data.rds the batch-corrected matrix?")
}

# ==============================================================================
# STEP 6  Define CI modules — cliques spanning all three datasets
# ==============================================================================

log_msg("=== STEP 6: Defining CI modules (3-dataset cliques) ===")
ci_modules <- list()

if (nrow(all_matches) > 0) {
  # Prefix nodes with dataset name to keep them unique across datasets
  edge_df <- all_matches %>%
    mutate(
      node_a = paste0(DS_A, "::", Module_A),
      node_b = paste0(DS_B, "::", Module_B)
    )
  
  g <- graph_from_data_frame(edge_df[, c("node_a", "node_b")], directed = FALSE)
  E(g)$weight <- edge_df$r
  
  # Find maximal cliques; keep those with one node per dataset
  all_cliques <- cliques(g, min = 3, max = 3)  # triangles = 3 datasets
  
  ci_idx <- 0
  for (cl in all_cliques) {
    node_names    <- V(g)$name[cl]
    ds_of_nodes   <- sub("::.*", "", node_names)
    if (length(unique(ds_of_nodes)) == 3 && all(datasets %in% ds_of_nodes)) {
      ci_idx <- ci_idx + 1
      ci_modules[[ci_idx]] <- node_names
      log_msg(paste("  CI-", ci_idx, ":", paste(node_names, collapse = " | ")))
    }
  }
  log_msg(paste("Total CI modules found:", length(ci_modules)))
} else {
  log_msg("No all-3-dataset cliques possible (no matches found).")
}

# ==============================================================================
# STEP 7  Extract gene membership for each CI module
# ==============================================================================
# A gene is included in a CI module if:
#   (a) it belongs to the matching Louvain module in that dataset (mc_list), AND
#   (b) its kME in that module exceeds KME_THRESHOLD (minimum membership strength)

extract_ci_genes <- function(ci_node_list, mc_list, kME_list, threshold = KME_THRESHOLD) {
  result_all <- data.frame()
  for (i in seq_along(ci_node_list)) {
    ci_name <- paste0("CI-", i)
    nodes   <- ci_node_list[[i]]
    per_ds  <- list()
    for (node in nodes) {
      parts   <- strsplit(node, "::")[[1]]
      ds      <- parts[1]
      me_col  <- parts[2]                  # e.g. "MEbrown"
      color   <- sub("^ME", "", me_col)    # e.g. "brown"
      
      mc  <- mc_list[[ds]]
      kME <- kME_list[[ds]]
      
      # Genes assigned to this module in Louvain results
      louvain_genes <- mc$Gene[mc$ModuleColor == color]
      
      # Filter by kME strength (me_col must exist as column name in kME)
      if (me_col %in% colnames(kME)) {
        kme_vals     <- kME[louvain_genes[louvain_genes %in% rownames(kME)], me_col]
        strong_genes <- names(kme_vals)[!is.na(kme_vals) & kme_vals >= threshold]
      } else {
        log_msg(paste("  WARNING:", me_col, "not in kME for", ds, "— using all Louvain genes"))
        strong_genes <- louvain_genes
      }
      
      per_ds[[ds]] <- data.frame(
        Gene    = strong_genes,
        Dataset = ds,
        Module  = color,
        kME     = if (me_col %in% colnames(kME))
          round(kME[strong_genes[strong_genes %in% rownames(kME)], me_col], 4)
        else NA_real_,
        CI      = ci_name,
        stringsAsFactors = FALSE
      )
    }
    result_all <- rbind(result_all, do.call(rbind, per_ds))
  }
  return(result_all)
}

log_msg("=== STEP 7: Extracting CI module gene lists ===")
if (length(ci_modules) > 0) {
  ci_genes_long <- extract_ci_genes(ci_modules, mc_list, kME_list)
  
  # Stability summary: per gene per CI module, how many datasets include it
  gene_stability <- ci_genes_long %>%
    group_by(Gene, CI) %>%
    summarise(
      N_Datasets      = n_distinct(Dataset),
      Datasets        = paste(sort(unique(Dataset)), collapse = "|"),
      Module_per_DS   = paste(paste0(Dataset, ":", Module), collapse = "|"),
      Mean_kME        = round(mean(kME, na.rm = TRUE), 4),
      Stability_Class = dplyr::case_when(
        N_Datasets == 3 ~ "Hard_core_3of3",
        N_Datasets == 2 ~ "Extended_2of3",
        TRUE            ~ "Dataset_specific_1of3"
      ),
      .groups = "drop"
    ) %>%
    arrange(CI, desc(N_Datasets), desc(Mean_kME))
  
  write.csv(gene_stability,
            file.path(out_dir, "intersection_gene_lists.csv"), row.names = FALSE)
  log_msg(paste("Gene stability table saved:",
                nrow(gene_stability), "gene-CI entries"))
  
  # Print hard-core summary per CI module
  hard_core <- gene_stability %>% filter(N_Datasets == 3)
  log_msg(paste("Hard core genes (3/3 datasets):", nrow(hard_core)))
  for (ci in unique(gene_stability$CI)) {
    hc <- hard_core$Gene[hard_core$CI == ci]
    log_msg(paste(" ", ci, "— hard core (n =", length(hc), "):",
                  paste(head(hc, 15), collapse = ", ")))
  }
} else {
  log_msg("No CI modules to extract genes for.")
  gene_stability <- data.frame()
}

# ==============================================================================
# STEP 8  Visualisation — cross-dataset kME correlation heatmaps
# ==============================================================================

plot_kME_heatmap <- function(cor_mat, title, file_path, threshold = R_THRESHOLD) {
  if (is.null(cor_mat) || nrow(cor_mat) == 0) return(invisible(NULL))
  # Clean "ME" prefix from labels
  rn <- sub("^ME", "", rownames(cor_mat))
  cn <- sub("^ME", "", colnames(cor_mat))
  rownames(cor_mat) <- rn
  colnames(cor_mat) <- cn
  
  breaks <- seq(-1, 1, length.out = 101)
  cols   <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  
  # Mark matched pairs
  display_mat <- round(cor_mat, 2)
  display_mat[abs(cor_mat) < threshold] <- NA  # grey out non-matches
  
  pdf(file_path, width = max(6, ncol(cor_mat) * 0.6 + 2),
      height = max(5, nrow(cor_mat) * 0.5 + 2))
  pheatmap(cor_mat,
           color           = cols,
           breaks          = breaks,
           display_numbers = round(cor_mat, 2),
           number_format   = "%.2f",
           fontsize_number = 7,
           main            = paste0(title, "\n(r > ", threshold, " = matched)"),
           cluster_rows    = TRUE,
           cluster_cols    = TRUE,
           border_color    = "grey90",
           na_col          = "grey90")
  dev.off()
  log_msg(paste("Heatmap saved:", basename(file_path)))
}

log_msg("=== STEP 8: Saving correlation heatmaps ===")
for (key in names(cor_mats)) {
  parts <- strsplit(key, "_vs_")[[1]]
  plot_kME_heatmap(
    cor_mats[[key]],
    title     = paste("kME-based module correlation:", parts[1], "vs", parts[2]),
    file_path = file.path(out_dir, paste0("cor_heatmap_", key, ".pdf"))
  )
}

# ==============================================================================
# STEP 9  Written summary report
# ==============================================================================

log_msg("=== STEP 9: Writing summary ===")
summary_lines <- c(
  "=================================================================",
  "MODULE 05: Cross-Dataset Intersection Module Summary",
  paste("Date:", Sys.time()),
  "=================================================================",
  "",
  paste("Datasets analysed:          ", paste(datasets, collapse = ", ")),
  paste("kME correlation threshold:  ", R_THRESHOLD, "(manuscript criterion: r > 0.7)"),
  paste("kME gene-inclusion cutoff:  ", KME_THRESHOLD),
  paste("Platform:                    GPL570 (HG-U133A) — all 3 datasets"),
  "",
  "--- Module matching results ---"
)

for (key in names(match_list)) {
  mt <- match_list[[key]]
  summary_lines <- c(summary_lines,
                     paste(" ", key, ":", nrow(mt), "matched pairs at |r| >", R_THRESHOLD))
  for (i in seq_len(min(5, nrow(mt))))
    summary_lines <- c(summary_lines,
                       paste0("    r = ", mt$r[i], "  |  ", mt$Module_A[i], "  <->  ", mt$Module_B[i]))
}

summary_lines <- c(summary_lines, "",
                   paste("--- CI modules defined:", length(ci_modules), "---"))

if (nrow(gene_stability) > 0) {
  for (ci in unique(gene_stability$CI)) {
    sub <- gene_stability[gene_stability$CI == ci, ]
    hc  <- sub$Gene[sub$N_Datasets == 3]
    ext <- sub$Gene[sub$N_Datasets == 2]
    summary_lines <- c(summary_lines,
                       paste(""),
                       paste(ci),
                       paste("  Total gene-dataset entries:", nrow(sub)),
                       paste("  Hard core (3/3):", length(hc), "genes:",
                             paste(head(hc, 20), collapse = ", ")),
                       paste("  Extended (2/3):", length(ext), "genes:",
                             paste(head(ext, 10), collapse = ", "))
    )
  }
}

summary_lines <- c(summary_lines, "",
                   "=================================================================",
                   "INTERPRETATION GUIDANCE FOR METHODS SECTION:",
                   "",
                   "Report as: 'Cross-dataset module unification was performed by computing",
                   "gene-level module membership (kME, Pearson r of each gene's expression",
                   "with each module eigengene) in each dataset independently. For each pair",
                   "of datasets, a module-pair correlation matrix was constructed by",
                   "correlating kME vectors across genes common to both datasets. Modules",
                   "from different datasets were declared equivalent when the Pearson",
                   "correlation of their kME vectors exceeded r > 0.7 across at least",
                   "2 of 3 pairwise comparisons, consistent with the eigengene-correlation",
                   "criterion used in the original WGCNA analysis. CI modules were defined",
                   "as sets of equivalent modules spanning all three datasets.'",
                   "=================================================================")

writeLines(summary_lines, file.path(out_dir, "intersection_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
log_msg("=== Module 05 complete ===")
log_msg(paste("Outputs saved to:", out_dir))

# ==============================================================================
# DIAGNOSTIC: If no CI modules found, print top correlations for inspection
# ==============================================================================
if (length(ci_modules) == 0 && length(cor_mats) > 0) {
  cat("\n=== DIAGNOSTIC: Top module-pair correlations (all pairs) ===\n")
  for (key in names(cor_mats)) {
    cm <- cor_mats[[key]]
    flat <- data.frame(
      Pair = paste0(rep(rownames(cm), ncol(cm)), "<->",
                    rep(colnames(cm), each = nrow(cm))),
      r    = as.vector(cm)
    )
    flat <- flat[!is.na(flat$r), ]
    flat <- flat[order(abs(flat$r), decreasing = TRUE), ]
    cat("\n", key, "— top 10 pairs:\n")
    print(head(flat, 10))
  }
  cat("\nIf top correlations are well below 0.7, possible causes:\n")
  cat("  1. checkpoint_raw_data.rds is the RAW (not batch-corrected) matrix.\n")
  cat("     Solution: In Module 03, batch correction is noted as 'skipped' for\n")
  cat("     individual datasets (each is normalised independently with RMA).\n")
  cat("     The three datasets share the same GPL570 platform, so inter-dataset\n")
  cat("     differences should be small after RMA. If correlations are still low,\n")
  cat("     consider applying ComBat across datasets before running Module 05.\n")
  cat("  2. The soft-power or min_cluster_size differed substantially across\n")
  cat("     datasets, producing incompatible module structures.\n")
  cat("  3. Grey (unassigned) genes were not filtered out — already handled above.\n")
}