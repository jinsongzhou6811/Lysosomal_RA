# ========================================================
# 702_Limma_TMT_Differential_Analysis.R
#
# Study:  PRIDE PXD027703 / P20180100010
# Paper:  Ren X et al., J Proteome Res. 2021, 20(10):4746-4757
# DOI:    10.1021/acs.jproteome.1c00399
# Title:  Quantitative Proteomic Analysis of Synovial Tissue Reveals That
#         Upregulated OLFM4 Aggravates Inflammation in Rheumatoid Arthritis
#
# Experimental design (confirmed from PRIDE metadata + PCA):
#   Tissue  : Human synovial tissue
#   Samples : RA patients (n=10) vs OA patients (n=12, control)
#   Design  : Every 3-4 patient samples pooled → 5 analytical channels
#   TMT 6-plex channel assignment (PCA-confirmed):
#     ch1, ch2, ch3  →  OA (control)  pools 1–3  (~12 patients)
#     ch4, ch6       →  RA (disease)  pools 1–2  (~10 patients)
#     ch5            →  EXCLUDED (outlier in PCA; likely pooled reference)
#   Note: OA is used as comparator (not healthy), so differential proteins
#         represent RA-specific changes relative to OA, not normal tissue.
#
# Inputs:
#   outputs/07_Proteomics_Analysis/01_Proteomics_Analysis/01c_tmt_normalized.csv
#     cols: ProteinID (int), Description (chr), TMT_ch1 … TMT_ch6 (num)
#           log2 + median-normalised intensities produced by script 610
#   outputs/07_Proteomics_Analysis/01_Proteomics_Analysis/03_protein_detected.csv
#     cols: Gene, ProteinIDs, UniProt_ACs, Module, In_top_combo, Passed_606
#   outputs/06_UVMR_Analysis/05_Module_Screening_Results/
#     {mod}/Final_Results/module_{mod}_final_top_combinations.csv
#     cols: gene, ivw_b, ivw_se, ivw_pval_BH, module, composite_score …
#   outputs/06_UVMR_Analysis/04_UVMR_Data_filtering/SingleGene_Top1/
#     SingleGene_Top1_All.csv
#     cols: exposure, b, se, pval_BH, score …
#
# Outputs:  outputs/07_Proteomics_Analysis/02_Differential_Proteomics/
#   01_limma_all_proteins.csv          Full limma results (all proteins)
#   02_target_genes_comparison.csv     Target gene TMT logFC vs MR β table
#   03_volcano.pdf / _300.png          Volcano plot (MR targets highlighted)
#   04_heatmap_targets.pdf             Heatmap of target genes across samples
#   05_MR_direction_consistency.csv    MR β vs TMT logFC direction summary
#   06_analysis_report.txt             Full analysis log
#
# Bugs fixed vs original draft:
#   [F1] rownames/ProteinID confusion  → ProteinID kept as column throughout
#   [F2] topTable rownames misused     → keep.rownames = "ProteinID_key"
#   [F3] res filtered to MR genes only → all proteins retained; annotation col added
#   [F4] heatmap match() logic error   → replaced with setNames lookup
#   [F5] NA values not handled         → group-aware filter + half-min imputation
#   [F6] sample grouping hard-coded    → confirmed from PRIDE + PCA; ch5 excluded
# ========================================================

suppressPackageStartupMessages({
  library(rprojroot)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

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

project_root <- find_rstudio_root_file()

# ══════════════════════════════════════════════════════════════════════════════
# SAMPLE CONFIGURATION  ← confirmed; modify only if design changes
# ══════════════════════════════════════════════════════════════════════════════

# ch5 excluded: PCA PC1 score = +495 (all others: -197 to +39).
# Likely pooled reference channel; excluded from differential analysis.
sample_info <- data.frame(
  Channel  = c("TMT_ch1", "TMT_ch2", "TMT_ch3", "TMT_ch4", "TMT_ch6"),
  Group    = c("OA",      "OA",      "OA",      "RA",      "RA"),
  SampleID = c("OA_pool1","OA_pool2","OA_pool3","RA_pool1","RA_pool2"),
  stringsAsFactors = FALSE
)

control_label  <- "OA"    # reference group
disease_label  <- "RA"    # comparison group

# Analysis thresholds
# Note: n=3 vs n=2 → low statistical power; exploratory thresholds applied
fc_threshold   <- 0.5     # |log2FC| cut-off
pval_threshold <- 0.1     # adj.P.Val cut-off (relaxed given small n)

# MR modules to query
modules <- c("GZMB", "LAMP3", "NKG7", "TRAF3IP3")

# ── Paths ─────────────────────────────────────────────────────────────────────
prot_dir <- file.path(project_root, "outputs", "07_Proteomics_Analysis",
                      "01_Proteomics_Analysis")
uvmr_dir <- file.path(project_root, "outputs", "06_UVMR_Analysis")
out_dir  <- file.path(project_root, "outputs", "07_Proteomics_Analysis",
                      "02_Differential_Proteomics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Log ───────────────────────────────────────────────────────────────────────
log_path <- file.path(out_dir, "06_analysis_report.txt")
log_con  <- file(log_path, open = "wt", encoding = "UTF-8")
on.exit(close(log_con), add = TRUE)

wlog <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., sep = ""))
  cat(msg, "\n")
  tryCatch(writeLines(msg, log_con), error = function(e) invisible())
}

safe_fread <- function(p, ...) {
  if (!file.exists(p)) return(data.table())
  tryCatch(fread(p, ...), error = function(e) data.table())
}

wlog("=== 611  Limma TMT Differential Analysis ===")
wlog("Study       : PRIDE PXD027703 / P20180100010")
wlog("Comparison  : ", disease_label, " vs ", control_label,
     " (synovial tissue, TMT 6-plex)")
wlog("Channels    : ", paste(sample_info$Channel, collapse = ", "))
wlog("Groups      : ", paste(sample_info$Group,   collapse = ", "))
wlog("Sample IDs  : ", paste(sample_info$SampleID,collapse = ", "))
wlog("ch5 excluded: outlier confirmed by PCA (PC1 = +495)")

# Validate grouping
n_ctrl <- sum(sample_info$Group == control_label)
n_case <- sum(sample_info$Group == disease_label)
wlog("n_", control_label, " = ", n_ctrl,
     " | n_", disease_label, " = ", n_case)

if (n_ctrl < 2 || n_case < 2)
  stop("Each group must have >= 2 samples. Check sample_info.")

wlog("WARNING: small n (", n_case, " vs ", n_ctrl,
     ") — results are exploratory; interpret with caution.")

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Load and prepare TMT normalised matrix
# File: 01c_tmt_normalized.csv
# Cols: ProteinID (int), Description (chr), TMT_ch1 … TMT_ch6 (num, log2)
# ══════════════════════════════════════════════════════════════════════════════
wlog("--- Step 1: Load TMT normalised matrix ---")

norm_file <- file.path(prot_dir, "01c_tmt_normalized.csv")
if (!file.exists(norm_file))
  stop("01c_tmt_normalized.csv not found — please run script 610 first.")

dat <- fread(norm_file)
wlog("Loaded: ", nrow(dat), " proteins x ", ncol(dat), " cols")

# Confirm expected columns exist
tmt_cols_all <- grep("^TMT_ch", names(dat), value = TRUE)
tmt_cols     <- intersect(sample_info$Channel, tmt_cols_all)
if (length(tmt_cols) < 2)
  stop("No matching TMT channel columns. Check sample_info$Channel vs file.")
wlog("Channels used: ", paste(tmt_cols, collapse = ", "))

# Extract gene symbol and UniProt AC from FASTA-style Description
# Format: ">sp|Q71U36|TBA1A_HUMAN ... GN=TUBA1A PE=1 SV=1"
dat[, Gene      := str_match(Description, "GN=([^\\s]+)")[, 2]]
dat[, UniProt_AC := str_match(Description,
                              "\\|([A-Z][0-9][A-Z0-9]{3}[0-9])\\|")[, 2]]
wlog("Gene symbol extracted: ",
     sum(!is.na(dat$Gene)), " / ", nrow(dat), " proteins")

# ── Step 1.5: Missing value handling [F5] ─────────────────────────────────────
wlog("--- Step 1.5: Missing value filtering and imputation ---")

expr_raw <- as.matrix(dat[, ..tmt_cols])
rownames(expr_raw) <- dat$ProteinID   # integer ProteinID as row key [F1]

n_before <- nrow(expr_raw)

# Remove proteins with all-NA in either group (cannot estimate group mean)
keep <- apply(expr_raw, 1, function(r) {
  ctrl_ok <- sum(!is.na(r[sample_info$Group == control_label])) >= 1
  case_ok <- sum(!is.na(r[sample_info$Group == disease_label])) >= 1
  ctrl_ok & case_ok
})
expr_raw  <- expr_raw[keep, , drop = FALSE]
dat_filt  <- dat[keep]
wlog("Proteins retained after group-NA filter: ",
     nrow(expr_raw), " / ", n_before)

# Half-minimum imputation (log2 space: subtract log2(2) = 1 from group min)
# Appropriate for missing-not-at-random (low-abundance proteins)
expr_imp <- expr_raw
for (i in seq_len(nrow(expr_imp))) {
  for (grp in c(control_label, disease_label)) {
    idx  <- which(sample_info$Group == grp)
    vals <- expr_imp[i, idx]
    na_i <- which(is.na(vals))
    if (length(na_i) > 0 && length(na_i) < length(idx)) {
      expr_imp[i, idx[na_i]] <- min(vals, na.rm = TRUE) - 1  # half-min in log2
    }
  }
}
n_imputed  <- sum(is.na(expr_raw))
pct_imputed <- round(100 * n_imputed / length(expr_raw), 2)
wlog("Imputed: ", n_imputed, " values (", pct_imputed, "% of matrix)")

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: limma differential analysis
# ══════════════════════════════════════════════════════════════════════════════
wlog("--- Step 2: limma differential analysis ---")

group_fac <- factor(sample_info$Group,
                    levels = c(control_label, disease_label))
design    <- model.matrix(~0 + group_fac)
colnames(design) <- levels(group_fac)

contrast_label <- paste0(disease_label, "_vs_", control_label)
cont_matrix    <- makeContrasts(
  contrasts = paste0(disease_label, " - ", control_label),
  levels    = design
)
colnames(cont_matrix) <- contrast_label

# trend = TRUE : intensity-dependent prior variance (appropriate for TMT)
# robust = TRUE: robust prior to downweight outlier proteins (critical for n=2)
fit  <- lmFit(expr_imp, design)
fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
wlog("eBayes completed | df.prior = ", round(fit2$df.prior, 2))

# Extract results — use keep.rownames to recover ProteinID [F2]
res_raw <- topTable(fit2, coef = contrast_label,
                    number = Inf, sort.by = "P")
res <- as.data.table(res_raw, keep.rownames = "ProteinID_key")
res[, ProteinID := as.integer(ProteinID_key)]

# Annotate with gene/description [F2: no rownames dependency]
anno_cols <- intersect(c("ProteinID", "Gene", "UniProt_AC", "Description"),
                       names(dat_filt))
res <- merge(res, dat_filt[, ..anno_cols], by = "ProteinID", all.x = TRUE)

# Significance flags
res[, Significant := abs(logFC) >= fc_threshold & adj.P.Val < pval_threshold]
res[, Direction   := fcase(
  logFC >  fc_threshold & adj.P.Val < pval_threshold, "Up",
  logFC < -fc_threshold & adj.P.Val < pval_threshold, "Down",
  default = "NS"
)]

n_up   <- sum(res$Direction == "Up",   na.rm = TRUE)
n_dn   <- sum(res$Direction == "Down", na.rm = TRUE)
wlog("Significant (|logFC|>=", fc_threshold,
     ", adj.P<", pval_threshold, "): ", n_up + n_dn,
     " proteins (Up=", n_up, ", Down=", n_dn, ")")

# ── Annotate MR target genes [F3: no filtering, only labelling] ──────────────
det <- safe_fread(file.path(prot_dir, "03_protein_detected.csv"))

target_genes_vec <- if (nrow(det) > 0 && "Gene" %in% names(det))
  det[Detected == TRUE, Gene] else character(0)
top_combo_genes  <- if (nrow(det) > 0 && "In_top_combo" %in% names(det))
  det[In_top_combo == TRUE, Gene] else character(0)

res[, Is_MR_target := Gene %in% target_genes_vec]
res[, Is_top_combo := Gene %in% top_combo_genes]

wlog("MR target genes present in results: ",
     sum(res$Is_MR_target, na.rm = TRUE))
wlog("MR target genes significant: ",
     sum(res$Is_MR_target & res$Significant, na.rm = TRUE))

fwrite(res, file.path(out_dir, "01_limma_all_proteins.csv"))
wlog("01_limma_all_proteins.csv saved: ", nrow(res), " proteins")

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: MR β vs TMT logFC — direction consistency analysis
# ══════════════════════════════════════════════════════════════════════════════
wlog("--- Step 3: MR direction consistency ---")

# Load MR combination results (script 608)
mr_combo <- rbindlist(lapply(modules, function(mod) {
  f <- file.path(uvmr_dir, "05_Module_Screening_Results", mod,
                 "Final_Results",
                 paste0("module_", mod, "_final_top_combinations.csv"))
  dt <- safe_fread(f)
  if (nrow(dt) == 0) return(data.table())
  dt[, Module_name := mod]; dt
}), fill = TRUE)

# Load single-gene MR results (script 606)
mr_single <- safe_fread(
  file.path(uvmr_dir, "04_UVMR_Data_filtering",
            "SingleGene_Top1", "SingleGene_Top1_All.csv"))

# Build comparison table for MR target genes
tmt_sub <- res[Is_MR_target == TRUE,
               .(Gene, UniProt_AC, logFC, AveExpr,
                 P.Value, adj.P.Val, Direction, Significant)]

if (nrow(tmt_sub) > 0) {
  
  # Merge single-gene MR results
  if (nrow(mr_single) > 0) {
    tmt_sub <- merge(
      tmt_sub,
      mr_single[, .(Gene        = exposure,
                    MR_single_b    = b,
                    MR_single_se   = se,
                    MR_single_pBH  = pval_BH,
                    MR_single_score = score)],
      by = "Gene", all.x = TRUE)
  }
  
  # Merge best combination MR result per gene
  if (nrow(mr_combo) > 0 && "gene" %in% names(mr_combo)) {
    combo_long <- rbindlist(lapply(seq_len(nrow(mr_combo)), function(i) {
      genes <- unlist(strsplit(mr_combo$gene[i], "_", fixed = TRUE))
      data.table(
        Gene          = genes,
        Combo_name    = mr_combo$gene[i],
        Module        = mr_combo$module[i],
        MR_combo_b    = mr_combo$ivw_b[i],
        MR_combo_se   = mr_combo$ivw_se[i],
        MR_combo_pBH  = mr_combo$ivw_pval_BH[i],
        MR_combo_score = mr_combo$composite_score[i])
    }), fill = TRUE)
    best_combo <- combo_long[order(-MR_combo_score)][!duplicated(Gene)]
    tmt_sub    <- merge(tmt_sub, best_combo, by = "Gene", all.x = TRUE)
  }
  
  # Direction labels
  # MR β > 0 → risk-increasing exposure → disease group over-expresses protein
  # MR β < 0 → protective exposure      → disease group under-expresses protein
  tmt_sub[, TMT_direction :=
            ifelse(logFC > 0, "Up_in_RA", "Down_in_RA")]
  tmt_sub[, MR_single_dir :=
            fcase(!is.na(MR_single_b) & MR_single_b > 0, "Risk_up",
                  !is.na(MR_single_b) & MR_single_b < 0, "Protective_down",
                  default = NA_character_)]
  tmt_sub[, MR_combo_dir :=
            fcase(!is.na(MR_combo_b)  & MR_combo_b  > 0, "Risk_up",
                  !is.na(MR_combo_b)  & MR_combo_b  < 0, "Protective_down",
                  default = NA_character_)]
  
  # Consistency: Risk_up + Up_in_RA  OR  Protective_down + Down_in_RA
  tmt_sub[, Single_consistent := fcase(
    MR_single_dir == "Risk_up"         & TMT_direction == "Up_in_RA",   TRUE,
    MR_single_dir == "Protective_down" & TMT_direction == "Down_in_RA", TRUE,
    !is.na(MR_single_dir),                                              FALSE,
    default = NA
  )]
  tmt_sub[, Combo_consistent := fcase(
    MR_combo_dir  == "Risk_up"         & TMT_direction == "Up_in_RA",   TRUE,
    MR_combo_dir  == "Protective_down" & TMT_direction == "Down_in_RA", TRUE,
    !is.na(MR_combo_dir),                                               FALSE,
    default = NA
  )]
  
  # Log per-gene results
  wlog("\nPer-gene direction consistency (MR beta vs TMT logFC):")
  wlog(sprintf("  %-10s | %8s | %10s | %-16s | %-8s | %-16s | %-8s | %s",
               "Gene","logFC","adj.P","MR_single_dir",
               "S_cons","MR_combo_dir","C_cons","Module"))
  for (i in seq_len(nrow(tmt_sub))) {
    r <- tmt_sub[i]
    wlog(sprintf("  %-10s | %+8.3f | %10.4f | %-16s | %-8s | %-16s | %-8s | %s",
                 r$Gene,
                 r$logFC, r$adj.P.Val,
                 ifelse(is.na(r$MR_single_dir), "—", r$MR_single_dir),
                 ifelse(is.na(r$Single_consistent), "NA",
                        as.character(r$Single_consistent)),
                 ifelse(is.na(r$MR_combo_dir),  "—", r$MR_combo_dir),
                 ifelse(is.na(r$Combo_consistent),  "NA",
                        as.character(r$Combo_consistent)),
                 ifelse(is.na(r$Module), "—", r$Module)))
  }
  
  n_s_cons <- sum(tmt_sub$Single_consistent == TRUE, na.rm = TRUE)
  n_c_cons <- sum(tmt_sub$Combo_consistent  == TRUE, na.rm = TRUE)
  n_s_eval <- sum(!is.na(tmt_sub$Single_consistent))
  n_c_eval <- sum(!is.na(tmt_sub$Combo_consistent))
  
  wlog("\nConsistency summary:")
  wlog("  Single-gene MR: ", n_s_cons, " / ", n_s_eval, " consistent")
  wlog("  Combo MR      : ", n_c_cons, " / ", n_c_eval, " consistent")
  
  fwrite(tmt_sub, file.path(out_dir, "02_target_genes_comparison.csv"))
  fwrite(
    tmt_sub[, .(Gene, TMT_direction, logFC, adj.P.Val, Significant,
                MR_single_dir, MR_single_b, Single_consistent,
                MR_combo_dir,  MR_combo_b,  Combo_consistent,
                Module, Combo_name)],
    file.path(out_dir, "05_MR_direction_consistency.csv"))
  wlog("02 and 05 saved.")
  
} else {
  wlog("No MR target genes found in differential results.")
  lapply(c("02_target_genes_comparison.csv",
           "05_MR_direction_consistency.csv"),
         function(f) fwrite(data.table(), file.path(out_dir, f)))
  tmt_sub <- data.table()
}

# ══════════════════════════════════════════════════════════════════════════════
# Step 4: Volcano plot
# All proteins as grey background; MR targets highlighted in red/orange
# ══════════════════════════════════════════════════════════════════════════════
wlog("--- Step 4: Volcano plot ---")

# Label: all MR target genes (significant ones get bold annotation)
res[, label_gene := ifelse(Is_MR_target, Gene, NA_character_)]

res[, point_cat := fcase(
  Is_MR_target & Significant, "MR target (sig)",
  Is_MR_target,               "MR target (ns)",
  Significant,                "Other sig",
  default = "NS"
)]

col_map   <- c("MR target (sig)" = "#C0392B",
               "MR target (ns)"  = "#E67E22",
               "Other sig"       = "#2980B9",
               "NS"              = "grey80")
size_map  <- c("MR target (sig)" = 3.5,
               "MR target (ns)"  = 3.0,
               "Other sig"       = 1.8,
               "NS"              = 0.8)
alpha_map <- c("MR target (sig)" = 1.0,
               "MR target (ns)"  = 0.9,
               "Other sig"       = 0.6,
               "NS"              = 0.25)

y_cap <- min(ceiling(-log10(min(res$adj.P.Val, na.rm = TRUE))) + 1, 30)

p_vol <- ggplot(res,
                aes(x     = logFC,
                    y     = -log10(adj.P.Val + 1e-300),
                    colour = point_cat,
                    size   = point_cat,
                    alpha  = point_cat)) +
  geom_point(shape = 16) +
  geom_text_repel(
    aes(label = label_gene),
    size          = 3.5,
    fontface      = "bold",
    max.overlaps  = 40,
    box.padding   = 0.4,
    point.padding = 0.3,
    segment.color = "grey40",
    segment.size  = 0.3,
    na.rm         = TRUE
  ) +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold),
             linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = -log10(pval_threshold),
             linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  scale_colour_manual(values = col_map,   name = NULL) +
  scale_size_manual(  values = size_map,  name = NULL) +
  scale_alpha_manual( values = alpha_map, name = NULL) +
  scale_y_continuous(limits = c(0, y_cap)) +
  labs(
    title    = paste0(disease_label, " vs ", control_label,
                      " — Synovial TMT proteomics (PXD027703)"),
    subtitle = paste0(
      "n = ", n_case, " RA pools vs ", n_ctrl, " OA pools",
      " | |log\u2082FC| \u2265 ", fc_threshold,
      ", adj.P < ", pval_threshold,
      " | Red = MR target + significant"),
    x        = paste0("log\u2082 fold-change (", disease_label,
                      " / ", control_label, ")"),
    y        = expression(-log[10](adj.P.Val))
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position    = "right",
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 9, colour = "grey40"),
    panel.grid.major.x = element_line(colour = "grey95", linewidth = 0.3),
    panel.grid.major.y = element_line(colour = "grey95", linewidth = 0.3)
  )

# [MODIFIED] Triple-format output via save_fig
save_fig(p_vol, file.path(out_dir, "03_volcano"), w_in = 10, h_in = 8)
wlog("03_volcano saved (.pdf/_300.png/_600.tiff)")

# ══════════════════════════════════════════════════════════════════════════════
# Step 5: Heatmap — MR target genes across all five samples
# [F4] Row labelling uses setNames lookup, not match(rownames, col)
# ══════════════════════════════════════════════════════════════════════════════
wlog("--- Step 5: Heatmap (MR target genes) ---")

tgt_res <- res[Is_MR_target == TRUE & !is.na(Gene)]

if (nrow(tgt_res) >= 2) {
  
  # [FIX] Deduplicate: keep only the best protein entry per gene (lowest adj.P.Val)
  # This prevents duplicate rownames when multiple ProteinIDs map to the same Gene
  tgt_res <- tgt_res[order(adj.P.Val)][!duplicated(Gene)]
  wlog("Heatmap target genes after dedup: ", nrow(tgt_res))
  
  # Extract rows from imputed matrix [F4]
  pid_chr  <- as.character(tgt_res$ProteinID)
  tgt_mat  <- expr_imp[rownames(expr_imp) %in% pid_chr, , drop = FALSE]
  
  # Replace ProteinID row names with gene symbols (now guaranteed unique)
  gene_lut       <- setNames(tgt_res$Gene, as.character(tgt_res$ProteinID))
  rownames(tgt_mat) <- gene_lut[rownames(tgt_mat)]
  
  # Replace channel column names with sample IDs
  colnames(tgt_mat) <- sample_info$SampleID
  
  # Column annotation (sample group)
  col_anno <- data.frame(Group    = sample_info$Group,
                         row.names = sample_info$SampleID)
  ann_col_colors <- list(
    Group = setNames(c("#2980B9", "#C0392B"),
                     c(control_label, disease_label))
  )
  
  # Row annotation (direction of TMT difference)
  dir_lut  <- setNames(tgt_res$Direction, gene_lut[as.character(tgt_res$ProteinID)])
  row_anno <- data.frame(TMT_direction = dir_lut[rownames(tgt_mat)],
                         row.names     = rownames(tgt_mat))
  ann_row_colors <- list(
    TMT_direction = c(Up = "#C0392B", Down = "#2980B9", NS = "grey75")
  )
  
  tryCatch({
    hm_stem <- file.path(out_dir, "04_heatmap_targets")
    hm_h    <- max(4, nrow(tgt_mat) * 0.5 + 2.5)
    hm_args <- list(
      mat                      = tgt_mat,
      scale                    = "row",
      clustering_distance_rows = "correlation",
      clustering_distance_cols = "euclidean",
      color                    = colorRampPalette(
        rev(brewer.pal(9, "RdBu")))(100),
      annotation_col           = col_anno,
      annotation_row           = row_anno,
      annotation_colors        = c(ann_col_colors, ann_row_colors),
      fontsize_row             = 9,
      fontsize_col             = 9,
      main                     = paste0("MR target proteins — ",
                                        disease_label, " vs ",
                                        control_label)
    )
    # [MODIFIED] Triple-format output for heatmap
    pdf(paste0(hm_stem, ".pdf"), width = 8, height = hm_h)
    do.call(pheatmap, hm_args); dev.off()
    png(paste0(hm_stem, "_300.png"), width = 8, height = hm_h, units = "in", res = 300)
    do.call(pheatmap, hm_args); dev.off()
    tiff(paste0(hm_stem, "_600.tiff"), width = 8, height = hm_h, units = "in", res = 600, compression = "lzw")
    do.call(pheatmap, hm_args); dev.off()
    wlog("04_heatmap_targets saved (", nrow(tgt_mat), " genes, .pdf/_300.png/_600.tiff)")
  }, error = function(e)
    wlog("Heatmap failed: ", conditionMessage(e)))
  
} else {
  wlog("Heatmap skipped: fewer than 2 MR target genes in results.")
}

# ══════════════════════════════════════════════════════════════════════════════
# Step 5.5: Combo-level TMT proteomics scoring
# Computes equal-weight mean of TMT intensities for top MR combinations
# Outputs: 06_combo_score_per_sample.csv, 07_combo_score_summary.csv
# ══════════════════════════════════════════════════════════════════════════════
wlog("--- Step 5.5: Combo-level TMT scoring ---")

if (nrow(mr_combo) > 0 && "gene" %in% names(mr_combo)) {
  
  # Get unique top combos (deduplicate by combo name, keep best score)
  top_combos <- mr_combo[order(-composite_score)][!duplicated(gene)]
  
  # Build gene-to-row lookup from imputed matrix
  gene_to_pid <- setNames(as.character(res$ProteinID), res$Gene)
  
  combo_per_sample <- list()
  combo_summary    <- list()
  
  for (i in seq_len(nrow(top_combos))) {
    combo_name  <- top_combos$gene[i]
    combo_genes <- unlist(strsplit(combo_name, "_", fixed = TRUE))
    mr_b        <- top_combos$ivw_b[i]
    mr_pBH      <- top_combos$ivw_pval_BH[i]
    mr_dir      <- ifelse(mr_b > 0, "Risk_up", "Protective_down")
    
    # Find which combo genes are detected in TMT data
    detected    <- combo_genes[combo_genes %in% res$Gene]
    if (length(detected) < 2) next
    
    # Get ProteinIDs for detected genes (best per gene, already deduped in res)
    res_dedup <- res[!duplicated(Gene)]
    pids <- as.character(res_dedup[Gene %in% detected, ProteinID])
    pids <- pids[pids %in% rownames(expr_imp)]
    if (length(pids) < 2) next
    
    # Equal-weight mean of TMT intensities per sample
    sub_mat    <- expr_imp[pids, , drop = FALSE]
    score_vec  <- colMeans(sub_mat, na.rm = TRUE)
    names(score_vec) <- sample_info$SampleID
    
    # Per-sample record
    combo_per_sample[[i]] <- data.table(
      Combo         = combo_name,
      Method        = "equal_weight_mean",
      N_genes       = length(pids),
      Genes_used    = paste(detected, collapse = ";"),
      MR_combo_b    = mr_b,
      MR_combo_pBH  = mr_pBH,
      MR_direction  = mr_dir,
      OA_pool1      = score_vec["OA_pool1"],
      OA_pool2      = score_vec["OA_pool2"],
      OA_pool3      = score_vec["OA_pool3"],
      RA_pool1      = score_vec["RA_pool1"],
      RA_pool2      = score_vec["RA_pool2"]
    )
    
    # Summary: mean per group + direction consistency
    oa_mean <- mean(score_vec[sample_info$Group == control_label], na.rm = TRUE)
    ra_mean <- mean(score_vec[sample_info$Group == disease_label], na.rm = TRUE)
    delta   <- round(ra_mean - oa_mean, 4)
    tmt_dir <- ifelse(delta > 0, "Up_in_RA", "Down_in_RA")
    consist <- (mr_dir == "Risk_up" & tmt_dir == "Up_in_RA") |
      (mr_dir == "Protective_down" & tmt_dir == "Down_in_RA")
    
    combo_summary[[i]] <- data.table(
      Combo          = combo_name,
      Method         = "equal_weight_mean",
      N_genes        = length(pids),
      Genes_used     = paste(detected, collapse = ";"),
      Score_OA_mean  = round(oa_mean, 4),
      Score_RA_mean  = round(ra_mean, 4),
      Delta_RA_OA    = delta,
      TMT_direction  = tmt_dir,
      MR_combo_b     = mr_b,
      MR_combo_pBH   = mr_pBH,
      MR_direction   = mr_dir,
      Consistent     = consist
    )
  }
  
  combo_per_sample_dt <- rbindlist(combo_per_sample[!sapply(combo_per_sample, is.null)])
  combo_summary_dt    <- rbindlist(combo_summary[!sapply(combo_summary, is.null)])
  
  if (nrow(combo_per_sample_dt) > 0) {
    fwrite(combo_per_sample_dt, file.path(out_dir, "06_combo_score_per_sample.csv"))
    fwrite(combo_summary_dt,    file.path(out_dir, "07_combo_score_summary.csv"))
    wlog("06_combo_score_per_sample.csv saved: ", nrow(combo_per_sample_dt), " combos")
    wlog("07_combo_score_summary.csv saved: ",    nrow(combo_summary_dt),    " combos")
    
    for (j in seq_len(nrow(combo_summary_dt))) {
      r <- combo_summary_dt[j]
      wlog(sprintf("  %s | OA=%.3f RA=%.3f | delta=%.4f | %s | MR=%s | %s",
                   r$Combo, r$Score_OA_mean, r$Score_RA_mean, r$Delta_RA_OA,
                   r$TMT_direction, r$MR_direction,
                   ifelse(r$Consistent, "CONSISTENT", "DISCORDANT")))
    }
  } else {
    wlog("No valid combos with >= 2 detected genes for scoring.")
  }
} else {
  wlog("No MR combo data available, skipping combo scoring.")
}

# ══════════════════════════════════════════════════════════════════════════════
# Final summary
# ══════════════════════════════════════════════════════════════════════════════
wlog("\n══════════ ANALYSIS SUMMARY ══════════")
wlog("Dataset          : PRIDE PXD027703 (P20180100010)")
wlog("Comparison       : ", disease_label, " vs ", control_label,
     " (synovial tissue)")
wlog("Channels used    : ", paste(sample_info$Channel, collapse = ", "),
     " (ch5 excluded)")
wlog("Proteins analysed: ", nrow(res))
wlog("Significant      : ", sum(res$Significant, na.rm = TRUE),
     " (|log2FC|>=", fc_threshold, ", adj.P<", pval_threshold, ")")
wlog("  Up in RA       : ", sum(res$Direction == "Up",   na.rm = TRUE))
wlog("  Down in RA     : ", sum(res$Direction == "Down", na.rm = TRUE))
wlog("MR targets found : ", sum(res$Is_MR_target, na.rm = TRUE))
wlog("MR targets sig   : ", sum(res$Is_MR_target & res$Significant, na.rm = TRUE))

if (nrow(tmt_sub) > 0) {
  wlog("Direction consistent (single MR): ", n_s_cons, " / ", n_s_eval)
  wlog("Direction consistent (combo MR) : ", n_c_cons, " / ", n_c_eval)
}

wlog("\n--- Key MR target gene results ---")
if (nrow(res[Is_MR_target == TRUE]) > 0) {
  key <- intersect(c("Gene","logFC","AveExpr","P.Value",
                     "adj.P.Val","Direction","Significant"),
                   names(res))
  tbl <- res[Is_MR_target == TRUE][order(adj.P.Val), ..key]
  for (i in seq_len(nrow(tbl))) {
    r <- tbl[i]
    wlog(sprintf("  %-10s | logFC=%+.3f | adj.P=%.4f | %-4s %s",
                 r$Gene, r$logFC, r$adj.P.Val,
                 r$Direction,
                 ifelse(isTRUE(r$Significant), "*** SIG", "")))
  }
}

wlog("\n--- Interpretation note ---")
wlog("  Control = OA (osteoarthritis), not healthy tissue.")
wlog("  Significant proteins represent RA-specific changes vs OA,")
wlog("  not RA vs normal. Interpret accordingly in the manuscript.")
wlog("  Recommended phrasing: 'differentially expressed between RA")
wlog("  and OA synovial tissue'.")

wlog("\n=== 611 completed ===")
wlog("Outputs → ", out_dir)

# Console summary
cat("\n╔══════════════ SUMMARY ══════════════╗\n")
cat(sprintf("  Proteins analysed : %d\n",  nrow(res)))
cat(sprintf("  Significant       : %d\n",  sum(res$Significant, na.rm=TRUE)))
cat(sprintf("    Up in RA        : %d\n",  n_up))
cat(sprintf("    Down in RA      : %d\n",  n_dn))
cat(sprintf("  MR targets        : %d\n",  sum(res$Is_MR_target, na.rm=TRUE)))
cat(sprintf("  MR targets sig    : %d\n",
            sum(res$Is_MR_target & res$Significant, na.rm=TRUE)))
cat("╚═════════════════════════════════════╝\n\n")
cat("Key output files:\n")
cat("  01_limma_all_proteins.csv       — full differential results\n")
cat("  02_target_genes_comparison.csv  — MR target gene details\n")
cat("  05_MR_direction_consistency.csv — MR beta vs TMT logFC direction\n")
cat("  03_volcano.pdf / _300.png\n")
cat("  04_heatmap_targets.pdf\n")
cat("\nOutputs →", out_dir, "\n")