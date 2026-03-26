# ========================================================
# 608_Module_BioValidation_and_SingleVsCombo.R  (v2 simplified)
#
# 精简说明：移除 Step1（统计筛选）和 Step2（ML排序）
#   - Step1/2 与 607 已有的筛选逻辑高度重叠，在 Top10 样本量下无实质价值
#   - 直接使用 607 Top10 输出作为分析对象
#
# 保留内容：
#   Step 3  GO-BP + KEGG 富集分析（功能通路支持）
#   Step 4  单基因(606) vs 最优组合(607) 因果效应对比
#   Step 5  整合 GO/KEGG 支持标记到最终结果表
#
# Input:
#   outputs/06_UVMR_Analysis/04_UVMR_Data_filtering/ComboWeights_Top10/
#     Top10_GZMB.csv / Top10_NKG7.csv
#   outputs/06_UVMR_Analysis/04_UVMR_Data_filtering/SingleGene_Top1/
#     SingleGene_Top1_All.csv
#   data/kegg_jp_link_hsa_pathway.txt
#   data/Module_{name}_All genes.txt
#
# Output:
#   outputs/06_UVMR_Analysis/05_Module_Screening_Results/
#     gene_match_report.csv
#     {module}/
#       module_{module}_data.csv
#       module_{module}_go_enrich.csv
#       module_{module}_kegg_enrich.csv
#       module_{module}_single_vs_combo.csv
#       module_{module}_final_top_combinations.csv
#       module_{module}_log.txt
#       Final_Results/
# ========================================================

suppressPackageStartupMessages({
  library(rprojroot)
  library(fs)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

project_root <- find_rstudio_root_file()

# ── Paths ──────────────────────────────────────────────────────────────────────
input_dir        <- file.path(project_root, "outputs", "06_UVMR_Analysis",
                              "04_UVMR_Data_filtering", "ComboWeights_Top10")
single_top1_file <- file.path(project_root, "outputs", "06_UVMR_Analysis",
                              "04_UVMR_Data_filtering", "SingleGene_Top1",
                              "SingleGene_Top1_All.csv")
kegg_local_file  <- file.path(project_root, "data",
                              "kegg_jp_link_hsa_pathway.txt")
output_dir       <- file.path(project_root, "outputs", "06_UVMR_Analysis",
                              "05_Module_Screening_Results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== 608 Biological Validation & Single-vs-Combo ===\n")
cat("Input dir        :", input_dir,        "\n")
cat("Single-gene file :", single_top1_file, "\n")
cat("Local KEGG file  :", kegg_local_file,  "\n")
cat("Output dir       :", output_dir,       "\n\n")

# ── Enrichment parameters ──────────────────────────────────────────────────────
enrich_pvalue_cutoff <- 0.05
enrich_qvalue_cutoff <- 0.20

# ── Load module gene definitions ───────────────────────────────────────────────
module_files <- c(
  "Module_GZMB_All genes.txt",
  "Module_NKG7_All genes.txt"
)

modules <- list()
for (f in module_files) {
  full_path <- file.path(project_root, "data", f)
  if (!file.exists(full_path)) stop("Module gene file not found: ", full_path)
  lines    <- trimws(readLines(full_path, warn = FALSE))
  lines    <- lines[nchar(lines) > 0]
  mod_name <- sub("^Module_(.*)_All genes\\.txt$", "\\1", basename(f))
  modules[[mod_name]] <- list(core = lines[1], extended = lines[-1])
}
cat("Modules loaded:", paste(names(modules), collapse = ", "), "\n\n")

# ── Load local KEGG pathway-gene mapping ──────────────────────────────────────
if (!file.exists(kegg_local_file))
  stop("Local KEGG file not found: ", kegg_local_file)

kegg_raw <- fread(kegg_local_file, sep = "\t", header = FALSE,
                  col.names = c("pathway_raw", "gene_raw"), encoding = "UTF-8")
kegg_db  <- kegg_raw[,
                     .(pathway_id = str_remove(pathway_raw, "^path:"),
                       entrez_id  = str_remove(gene_raw,    "^hsa:"))
]
cat("Local KEGG loaded:",
    nrow(kegg_db), "pairs |",
    uniqueN(kegg_db$pathway_id), "pathways |",
    uniqueN(kegg_db$entrez_id),  "genes\n\n")

# ── Helper: local KEGG enrichment ─────────────────────────────────────────────
run_kegg_local <- function(candidate_entrez, universe_entrez,
                           pval_cutoff, qval_cutoff, entrez2symbol) {
  if (length(candidate_entrez) == 0 || length(universe_entrez) == 0)
    return(data.table())
  
  kegg_univ <- kegg_db[entrez_id %in% universe_entrez]
  kegg_cand <- kegg_db[entrez_id %in% candidate_entrez]
  n_cand    <- length(unique(candidate_entrez))
  n_univ    <- length(unique(universe_entrez))
  if (nrow(kegg_univ) == 0) return(data.table())
  
  results <- rbindlist(lapply(unique(kegg_univ$pathway_id), function(pid) {
    pw_univ <- kegg_univ[pathway_id == pid, unique(entrez_id)]
    pw_cand <- kegg_cand[pathway_id == pid, unique(entrez_id)]
    k <- length(pw_cand); M <- length(pw_univ)
    if (k == 0 || M == 0) return(NULL)
    pval <- phyper(k - 1L, M, n_univ - M, n_cand, lower.tail = FALSE)
    syms <- entrez2symbol[pw_cand]; syms <- syms[!is.na(syms)]
    data.table(pathway_id = pid, k = k, M = M, n = n_cand, N = n_univ,
               GeneRatio  = paste0(k, "/", n_cand),
               BgRatio    = paste0(M, "/", n_univ),
               pvalue     = pval,
               geneID     = paste(sort(syms), collapse = "/"))
  }), fill = TRUE)
  
  if (nrow(results) == 0) return(data.table())
  results[, p.adjust := p.adjust(pvalue, method = "BH")]
  results[, qvalue   := p.adjust]
  results[pvalue <= pval_cutoff & p.adjust <= qval_cutoff] %>%
    as.data.frame() %>%
    dplyr::arrange(pvalue) %>%
    as.data.table()
}

# ── Load 607 Top10 files ───────────────────────────────────────────────────────
top10_files <- fs::dir_ls(input_dir, regexp = "Top10_.*\\.csv$")
if (length(top10_files) == 0)
  stop("No Top10_*.csv files found in: ", input_dir)

raw_data <- rbindlist(lapply(top10_files, fread), fill = TRUE)

required_cols <- c("module", "exposure", "outcome", "b", "se",
                   "pval_BH", "nsnp", "mean_F", "score")
missing_cols  <- setdiff(required_cols, names(raw_data))
if (length(missing_cols) > 0)
  stop("Missing columns in Top10 input: ", paste(missing_cols, collapse = ", "))

combo_data <- as.data.frame(raw_data) %>%
  dplyr::rename(gene            = exposure,
                composite_score = score,
                ivw_b           = b,
                ivw_se          = se,
                ivw_pval_BH     = pval_BH) %>%
  as.data.table()

cat("607 combo data loaded:", nrow(combo_data), "rows |",
    "modules:", paste(sort(unique(combo_data$module)), collapse = ", "), "\n\n")

# ── Gene-module match report ───────────────────────────────────────────────────
all_combo_genes   <- unique(unlist(str_split(combo_data$gene, "_")))
gene_match_report <- rbindlist(lapply(names(modules), function(mod_name) {
  mod_genes <- c(modules[[mod_name]]$core, modules[[mod_name]]$extended)
  data.table(Module  = mod_name,
             Gene    = mod_genes,
             Matched = mod_genes %in% all_combo_genes)
}))
fwrite(gene_match_report, file.path(output_dir, "gene_match_report.csv"))
cat("Gene match report saved.\n\n")

# ── Load 606 single-gene Top1 ─────────────────────────────────────────────────
single_data <- if (file.exists(single_top1_file)) {
  fread(single_top1_file)
} else {
  warning("SingleGene_Top1_All.csv not found; Step 4 will be skipped.")
  NULL
}


# ═══════════════════════════════════════════════════════════════════════════════
# process_module() — Step3 + Step4 + Step5
# ═══════════════════════════════════════════════════════════════════════════════
process_module <- function(module_df, module_name, output_dir,
                           module_def, single_data) {
  
  module_output_dir <- file.path(output_dir, module_name)
  final_results_dir <- file.path(module_output_dir, "Final_Results")
  dir.create(final_results_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ── Log ────────────────────────────────────────────────────────────────────
  log_path <- file.path(module_output_dir,
                        paste0("module_", module_name, "_log.txt"))
  con <- tryCatch(file(log_path, open = "wt", encoding = "UTF-8"),
                  error = function(e) NULL)
  on.exit({ if (!is.null(con)) close(con) }, add = TRUE)
  
  log_write <- function(...) {
    msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., sep = ""))
    cat(msg, "\n")
    if (!is.null(con)) writeLines(msg, con)
  }
  
  log_write("=== Module: ", module_name, " ===")
  log_write("Input rows (607 Top10): ", nrow(module_df))
  
  # Save raw 607 data for this module
  fwrite(module_df,
         file.path(module_output_dir,
                   paste0("module_", module_name, "_data.csv")))
  
  # Combinations sorted by composite_score (607 ranking preserved)
  top_combinations <- as.data.frame(module_df) %>%
    dplyr::arrange(desc(composite_score))
  
  # ── Step 3: GO-BP + KEGG enrichment ────────────────────────────────────────
  log_write("--- Step 3: GO + KEGG enrichment ---")
  
  candidate_genes <- unique(unlist(str_split(top_combinations$gene, "_")))
  universe_genes  <- unique(c(module_def$core, module_def$extended))
  log_write("Candidate genes (", length(candidate_genes), "): ",
            paste(candidate_genes, collapse = ", "))
  log_write("Universe genes  (", length(universe_genes),  "): ",
            paste(universe_genes,  collapse = ", "))
  
  g2e <- tryCatch(
    bitr(candidate_genes, fromType = "SYMBOL",
         toType = "ENTREZID", OrgDb = org.Hs.eg.db),
    error = function(e) {
      log_write("bitr (candidates) failed: ", conditionMessage(e))
      data.frame(SYMBOL = character(), ENTREZID = character())
    }
  )
  u2e <- tryCatch(
    bitr(universe_genes, fromType = "SYMBOL",
         toType = "ENTREZID", OrgDb = org.Hs.eg.db),
    error = function(e) {
      log_write("bitr (universe) failed: ", conditionMessage(e))
      data.frame(SYMBOL = character(), ENTREZID = character())
    }
  )
  
  unmapped <- setdiff(candidate_genes, g2e$SYMBOL)
  if (length(unmapped) > 0)
    log_write("Genes unmapped to Entrez: ", paste(unmapped, collapse = ", "))
  
  go_result   <- data.table()
  kegg_result <- data.table()
  
  if (nrow(g2e) >= 2) {
    
    # GO-BP
    go_res <- tryCatch(
      enrichGO(gene          = g2e$ENTREZID,
               universe      = u2e$ENTREZID,
               OrgDb         = org.Hs.eg.db,
               ont           = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff  = enrich_pvalue_cutoff,
               qvalueCutoff  = enrich_qvalue_cutoff,
               readable      = TRUE),
      error = function(e) {
        log_write("GO failed: ", conditionMessage(e)); NULL
      }
    )
    if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
      go_result <- as.data.table(go_res)
      log_write("GO BP: ", nrow(go_result), " significant terms")
    } else {
      log_write("GO BP: no significant terms")
    }
    fwrite(go_result,
           file.path(module_output_dir,
                     paste0("module_", module_name, "_go_enrich.csv")))
    
    # KEGG (local)
    entrez2symbol <- setNames(g2e$SYMBOL, g2e$ENTREZID)
    kegg_result <- tryCatch(
      run_kegg_local(g2e$ENTREZID, u2e$ENTREZID,
                     enrich_pvalue_cutoff, enrich_qvalue_cutoff,
                     entrez2symbol),
      error = function(e) {
        log_write("KEGG failed: ", conditionMessage(e)); data.table()
      }
    )
    if (nrow(kegg_result) > 0) {
      log_write("KEGG: ", nrow(kegg_result), " significant pathways")
    } else {
      log_write("KEGG: no significant pathways")
    }
    fwrite(kegg_result,
           file.path(module_output_dir,
                     paste0("module_", module_name, "_kegg_enrich.csv")))
    
  } else {
    log_write("Enrichment skipped: < 2 candidate genes mapped to Entrez ID")
    fwrite(data.table(),
           file.path(module_output_dir,
                     paste0("module_", module_name, "_go_enrich.csv")))
    fwrite(data.table(),
           file.path(module_output_dir,
                     paste0("module_", module_name, "_kegg_enrich.csv")))
  }
  
  # ── Step 4: Single-gene (606) vs Top-1 combination (607) ───────────────────
  log_write("--- Step 4: Single-gene vs combination comparison ---")
  
  single_vs_combo_path <- file.path(
    module_output_dir,
    paste0("module_", module_name, "_single_vs_combo.csv")
  )
  
  if (!is.null(single_data) && nrow(top_combinations) > 0) {
    
    all_module_genes <- unique(c(module_def$core, module_def$extended))
    combo_top        <- top_combinations[1, ]
    
    log_write("Top-1 combo : ", combo_top$gene,
              " | outcome: ",   combo_top$outcome,
              " | OR = ",       round(combo_top$OR, 4),
              " | pval_BH = ",  round(combo_top$ivw_pval_BH, 4))
    
    comp_rows <- rbindlist(lapply(all_module_genes, function(g) {
      
      sg <- as.data.frame(single_data) %>%
        dplyr::filter(exposure == g) %>%
        dplyr::slice_max(score, n = 1, with_ties = FALSE)
      
      if (nrow(sg) > 0) {
        score_imp  <- if (isTRUE(sg$score[1] > 0))
          combo_top$composite_score / sg$score[1] else NA_real_
        effect_imp <- if (isTRUE(abs(sg$b[1]) > 0))
          abs(combo_top$ivw_b) / abs(sg$b[1]) else NA_real_
        
        data.table(
          Module                 = module_name,
          Single_gene            = g,
          Single_passed_606      = TRUE,
          Single_outcome         = sg$outcome[1],
          Single_b               = sg$b[1],
          Single_OR              = sg$OR[1],
          Single_OR_lower95      = sg$OR_lower95[1],
          Single_OR_upper95      = sg$OR_upper95[1],
          Single_pval_BH         = sg$pval_BH[1],
          Single_nsnp            = sg$nsnp[1],
          Single_mean_F          = sg$mean_F[1],
          Single_het_Q_pval      = sg$het_Q_pval[1],
          Single_egger_flag      = sg$egger_flag[1],
          Single_presso_sig      = sg$presso_sig[1],
          Single_direction_ratio = sg$direction_ratio[1],
          Single_score           = sg$score[1],
          Best_combo             = combo_top$gene,
          Combo_outcome          = combo_top$outcome,
          Combo_b                = combo_top$ivw_b,
          Combo_OR               = combo_top$OR,
          Combo_OR_lower95       = combo_top$OR_lower95,
          Combo_OR_upper95       = combo_top$OR_upper95,
          Combo_pval_BH          = combo_top$ivw_pval_BH,
          Combo_nsnp             = combo_top$nsnp,
          Combo_mean_F           = combo_top$mean_F,
          Combo_het_Q_pval       = combo_top$het_Q_pval,
          Combo_egger_flag       = combo_top$egger_flag,
          Combo_presso_sig       = combo_top$presso_sig,
          Combo_direction_ratio  = combo_top$direction_ratio,
          Combo_score            = combo_top$composite_score,
          Score_Improvement      = score_imp,
          Effect_Improvement     = effect_imp
        )
        
      } else {
        
        data.table(
          Module                 = module_name,
          Single_gene            = g,
          Single_passed_606      = FALSE,
          Single_outcome         = NA_character_,
          Single_b               = NA_real_,
          Single_OR              = NA_real_,
          Single_OR_lower95      = NA_real_,
          Single_OR_upper95      = NA_real_,
          Single_pval_BH         = NA_real_,
          Single_nsnp            = NA_real_,
          Single_mean_F          = NA_real_,
          Single_het_Q_pval      = NA_real_,
          Single_egger_flag      = NA,
          Single_presso_sig      = NA,
          Single_direction_ratio = NA_real_,
          Single_score           = NA_real_,
          Best_combo             = combo_top$gene,
          Combo_outcome          = combo_top$outcome,
          Combo_b                = combo_top$ivw_b,
          Combo_OR               = combo_top$OR,
          Combo_OR_lower95       = combo_top$OR_lower95,
          Combo_OR_upper95       = combo_top$OR_upper95,
          Combo_pval_BH          = combo_top$ivw_pval_BH,
          Combo_nsnp             = combo_top$nsnp,
          Combo_mean_F           = combo_top$mean_F,
          Combo_het_Q_pval       = combo_top$het_Q_pval,
          Combo_egger_flag       = combo_top$egger_flag,
          Combo_presso_sig       = combo_top$presso_sig,
          Combo_direction_ratio  = combo_top$direction_ratio,
          Combo_score            = combo_top$composite_score,
          Score_Improvement      = NA_real_,
          Effect_Improvement     = NA_real_
        )
      }
    }), fill = TRUE)
    
    fwrite(comp_rows, single_vs_combo_path)
    
    n_passed <- sum(comp_rows$Single_passed_606)
    n_total  <- nrow(comp_rows)
    log_write("Module genes total:                    ", n_total)
    log_write("Passed 606 (single-gene significant):  ", n_passed)
    log_write("NOT passed 606 (combo-only effect):    ", n_total - n_passed)
    
  } else {
    fwrite(data.table(), single_vs_combo_path)
    log_write("Step 4 skipped (single_data missing or no combo results)")
  }
  
  # ── Step 5: Final table with GO/KEGG support flags ──────────────────────────
  log_write("--- Step 5: Final integration ---")
  
  go_genes <- if (nrow(go_result) > 0 && "geneID" %in% names(go_result))
    unique(unlist(str_split(go_result$geneID, "/"))) else character(0)
  kegg_genes <- if (nrow(kegg_result) > 0 && "geneID" %in% names(kegg_result))
    unique(unlist(str_split(kegg_result$geneID, "/"))) else character(0)
  
  final_top <- as.data.frame(top_combinations) %>%
    dplyr::mutate(
      has_go_support   = sapply(str_split(gene, "_"),
                                function(g) any(g %in% go_genes)),
      has_kegg_support = sapply(str_split(gene, "_"),
                                function(g) any(g %in% kegg_genes)),
      has_bio_support  = has_go_support | has_kegg_support
    ) %>%
    dplyr::arrange(desc(composite_score))
  
  fwrite(final_top,
         file.path(module_output_dir,
                   paste0("module_", module_name, "_final_top_combinations.csv")))
  
  n_bio <- sum(final_top$has_bio_support)
  log_write("Final table saved: ", nrow(final_top), " rows | ",
            "with bio support: ", n_bio)
  
  # ── Copy to Final_Results/ ──────────────────────────────────────────────────
  key_suffixes <- c("go_enrich", "kegg_enrich",
                    "single_vs_combo", "final_top_combinations")
  for (suf in key_suffixes) {
    src  <- file.path(module_output_dir,
                      paste0("module_", module_name, "_", suf, ".csv"))
    dest <- file.path(final_results_dir,
                      paste0("module_", module_name, "_", suf, ".csv"))
    if (file.exists(src)) {
      file.copy(src, dest, overwrite = TRUE)
    } else {
      log_write("Warning: not found for copy: ",
                paste0("module_", module_name, "_", suf, ".csv"))
    }
  }
  
  log_write("=== Module ", module_name, " completed ===")
  return(invisible(final_top))
}


# ═══════════════════════════════════════════════════════════════════════════════
# Main loop
# ═══════════════════════════════════════════════════════════════════════════════
cat("Starting module processing loop\n\n")

for (module_name in names(modules)) {
  
  module_df <- combo_data[module == module_name]
  
  cat("──────────────────────────────────────────\n")
  cat("Module:", module_name, "| rows:", nrow(module_df), "\n")
  
  if (nrow(module_df) == 0) {
    cat("  No 607 results for this module. Skipping.\n")
    next
  }
  
  tryCatch(
    process_module(
      module_df   = module_df,
      module_name = module_name,
      output_dir  = output_dir,
      module_def  = modules[[module_name]],
      single_data = single_data
    ),
    error = function(e) {
      cat("  ERROR in module", module_name, ":", conditionMessage(e), "\n")
    }
  )
}

cat("\n=== 608 completed ===\n")
cat("Results saved to:", output_dir, "\n")