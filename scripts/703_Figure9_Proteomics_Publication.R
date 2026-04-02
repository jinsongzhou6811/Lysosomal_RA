# ============================================================
# 703_Figure9_Proteomics_Publication.R  (v3 — data-aware fixes, Figure 9)
# ============================================================

suppressPackageStartupMessages({
  library(rprojroot)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
})

# ── Font fallback ─────────────────────────────────────────────────────────────
BASE_FONT <- tryCatch({
  if (requireNamespace("extrafont", quietly=TRUE)) {
    extrafont::loadfonts(device="pdf",        quiet=TRUE)
    extrafont::loadfonts(device="postscript", quiet=TRUE)
    if ("Arial" %in% extrafont::fonts()) "Arial" else "sans"
  } else "sans"
}, error=function(e) "sans")
message("Font: ", BASE_FONT)

project_root <- find_rstudio_root_file()

out_dir <- file.path(project_root, "outputs", "07_Proteomics_Analysis",
                     "03_Figure9")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

log_con <- file(file.path(out_dir, "figure9_log.txt"),
                open="wt", encoding="UTF-8")
on.exit(close(log_con), add=TRUE)
wlog <- function(...) {
  msg <- paste0(format(Sys.time(),"%H:%M:%S")," | ", paste(...,sep=""))
  cat(msg,"\n"); writeLines(msg, log_con)
}
wlog("=== 703 v3 started === font=", BASE_FONT)

# ── Helpers ───────────────────────────────────────────────────────────────────
safe_fread <- function(p) {
  if (!file.exists(p)) { wlog("MISSING: ",p); return(data.table()) }
  tryCatch(fread(p), error=function(e){ wlog("ERR ",basename(p)); data.table() })
}

save_fig <- function(p, stem, w_in, h_in) {
  ggsave(paste0(stem,".pdf"),      plot=p, width=w_in, height=h_in,
         units="in", device="pdf", useDingbats=FALSE)
  ggsave(paste0(stem,"_300.png"),  plot=p, width=w_in, height=h_in,
         units="in", device="png",  dpi=300)
  ggsave(paste0(stem,"_600.tiff"), plot=p, width=w_in, height=h_in,
         units="in", device="tiff", dpi=600, compression="lzw")
  wlog("Saved: ", basename(stem))
}

theme_pub <- function(bs=11) {
  theme_classic(base_size=bs) +
    theme(
      text               = element_text(family=BASE_FONT),
      plot.title         = element_text(face="bold",size=bs+1,family=BASE_FONT),
      plot.subtitle      = element_text(size=bs-2,colour="grey40",family=BASE_FONT),
      axis.text          = element_text(size=bs-2,family=BASE_FONT),
      strip.background   = element_rect(fill="grey92",colour=NA),
      strip.text         = element_text(face="bold",size=bs-1,family=BASE_FONT),
      legend.text        = element_text(size=bs-2,family=BASE_FONT),
      legend.title       = element_text(size=bs-1,family=BASE_FONT),
      panel.grid.major.x = element_line(colour="grey92",linewidth=0.3),
      panel.grid.major.y = element_line(colour="grey92",linewidth=0.3)
    )
}

COL_MR_SIG  <- "#C0392B"
COL_MR_NS   <- "#E67E22"
COL_SIG     <- "#2980B9"
COL_NS      <- "grey78"
COL_CONSIST <- "#27AE60"
COL_INCONS  <- "#C0392B"
COL_NOMB    <- "#BDC3C7"

FC_THR   <- 0.5
PVAL_THR <- 0.1

sample_info <- data.frame(
  Channel  = c("TMT_ch1","TMT_ch2","TMT_ch3","TMT_ch4","TMT_ch6"),
  Group    = c("OA","OA","OA","RA","RA"),
  SampleID = c("OA_pool1","OA_pool2","OA_pool3","RA_pool1","RA_pool2"),
  stringsAsFactors=FALSE
)

# ── Paths ─────────────────────────────────────────────────────────────────────
prot_dir <- file.path(project_root,"outputs","07_Proteomics_Analysis",
                      "01_Proteomics_Analysis")
diff_dir <- file.path(project_root,"outputs","07_Proteomics_Analysis",
                      "02_Differential_Proteomics")

limma_all <- safe_fread(file.path(diff_dir,"01_limma_all_proteins.csv"))
tgt_cmp   <- safe_fread(file.path(diff_dir,"02_target_genes_comparison.csv"))
dir_cons  <- safe_fread(file.path(diff_dir,"05_MR_direction_consistency.csv"))
tmt_norm  <- safe_fread(file.path(prot_dir, "01c_tmt_normalized.csv"))
det       <- safe_fread(file.path(prot_dir, "03_protein_detected.csv"))

wlog("Raw rows — limma:",nrow(limma_all)," tgt_cmp:",nrow(tgt_cmp),
     " dir_cons:",nrow(dir_cons)," tmt_norm:",nrow(tmt_norm)," det:",nrow(det))

# ── [D1] Add Gene column to tmt_norm via Description parsing ─────────────────
if (nrow(tmt_norm) > 0 && !"Gene" %in% names(tmt_norm)) {
  if ("Description" %in% names(tmt_norm)) {
    tmt_norm[, Gene := str_match(Description, "GN=([^\\s;]+)")[,2]]
    wlog("tmt_norm: parsed Gene from Description, non-NA=",
         sum(!is.na(tmt_norm$Gene)))
  } else {
    wlog("tmt_norm: no Gene or Description column — heatmap will be skipped")
  }
}

# ── [D2] Deduplicate dir_cons: one row per gene ───────────────────────────────
# Strategy: lowest adj.P.Val; for consistency cols take majority vote
dedup_gene <- function(dt) {
  if (nrow(dt) == 0 || !"Gene" %in% names(dt)) return(dt)
  
  # For logical consistency columns, majority vote (TRUE if >= 50% TRUE)
  cons_cols <- intersect(c("Single_consistent","Combo_consistent"), names(dt))
  
  dt[, .SD[which.min(adj.P.Val)], by=Gene] -> best_row
  
  if (length(cons_cols) > 0) {
    vote <- dt[, lapply(.SD, function(x) {
      nT <- sum(x == TRUE,  na.rm=TRUE)
      nF <- sum(x == FALSE, na.rm=TRUE)
      if (nT + nF == 0) return(NA)
      nT >= nF
    }), .SDcols=cons_cols, by=Gene]
    setnames(vote, cons_cols, paste0(cons_cols,"_vote"))
    best_row <- merge(best_row, vote, by="Gene", all.x=TRUE)
    for (cc in cons_cols) {
      vc <- paste0(cc,"_vote")
      if (vc %in% names(best_row)) {
        best_row[[cc]] <- best_row[[vc]]
        best_row[, (vc) := NULL]
      }
    }
  }
  best_row
}

dir_cons_d <- dedup_gene(dir_cons)
tgt_cmp_d  <- dedup_gene(tgt_cmp)
wlog("After dedup — dir_cons:",nrow(dir_cons_d)," tgt_cmp:",nrow(tgt_cmp_d))

# ── [D3] Deduplicate limma_all: one row per Gene ──────────────────────────────
limma_d <- if (nrow(limma_all) > 0 && "Gene" %in% names(limma_all)) {
  tmp <- as.data.table(limma_all)
  tmp[!is.na(Gene)][order(adj.P.Val)][!duplicated(Gene)]
} else as.data.table(limma_all)
wlog("limma_all after dedup:",nrow(limma_d)," unique genes")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 9A — Volcano Plot
# ══════════════════════════════════════════════════════════════════════════════
wlog("--- 9A ---")

make_9A <- function(limma_d, tgt_cmp_d, fc_thr, pval_thr) {
  if (nrow(limma_d) == 0) {
    return(ggplot()+theme_void()+
             annotate("text",x=0.5,y=0.5,label="No limma data",
                      family=BASE_FONT,size=5))
  }
  
  vol <- copy(limma_d)
  # Flexible adj.P.Val column name
  if (!"adj.P.Val" %in% names(vol)) {
    ac <- grep("adj.*[Pp]|FDR|fdr", names(vol), value=TRUE)[1]
    if (!is.na(ac)) setnames(vol, ac, "adj.P.Val")
  }
  
  mr_genes  <- if (nrow(tgt_cmp_d)>0) unique(tgt_cmp_d$Gene) else character(0)
  sig_genes <- if (nrow(tgt_cmp_d)>0 && "Significant" %in% names(tgt_cmp_d))
    tgt_cmp_d[Significant==TRUE, Gene] else character(0)
  
  vol[, Is_MR := Gene %in% mr_genes]
  vol[, cat := fcase(
    Gene %in% sig_genes,                               "MR target (sig)",
    Is_MR & !Gene %in% sig_genes,                     "MR target (ns)",
    !Is_MR & abs(logFC)>=fc_thr & adj.P.Val<pval_thr, "Other sig",
    default = "NS"
  )]
  # [D4] Label each MR target gene ONCE
  vol[, label := ifelse(Is_MR & !duplicated(Gene), Gene, NA_character_)]
  
  col_map   <- c("MR target (sig)"=COL_MR_SIG,"MR target (ns)"=COL_MR_NS,
                 "Other sig"=COL_SIG,"NS"=COL_NS)
  size_map  <- c("MR target (sig)"=3.2,"MR target (ns)"=2.8,
                 "Other sig"=1.2,"NS"=0.6)
  alpha_map <- c("MR target (sig)"=1.0,"MR target (ns)"=0.85,
                 "Other sig"=0.50,"NS"=0.18)
  
  n_up <- sum(vol$logFC >= fc_thr & vol$adj.P.Val < pval_thr, na.rm=TRUE)
  n_dn <- sum(vol$logFC <= -fc_thr & vol$adj.P.Val < pval_thr, na.rm=TRUE)
  y_cap <- pmin(ceiling(-log10(min(vol$adj.P.Val+1e-300,na.rm=TRUE)))+1, 20)
  x_max <- max(abs(vol$logFC), na.rm=TRUE)
  
  ggplot(vol[order(cat)],
         aes(x=logFC, y=-log10(adj.P.Val+1e-300),
             colour=cat, size=cat, alpha=cat)) +
    geom_point(shape=16) +
    geom_text_repel(
      aes(label=label),
      size=3.2, fontface="bold", family=BASE_FONT,
      max.overlaps=50,
      box.padding=0.5, point.padding=0.3,
      segment.colour="grey45", segment.size=0.3,
      na.rm=TRUE
    ) +
    geom_vline(xintercept=c(-fc_thr,fc_thr),
               linetype="dashed",colour="grey50",linewidth=0.4) +
    geom_hline(yintercept=-log10(pval_thr),
               linetype="dashed",colour="grey50",linewidth=0.4) +
    annotate("text",x=x_max*0.85,y=y_cap*0.96,
             label=paste0("Up: ",n_up),
             colour=COL_MR_SIG,size=3.0,family=BASE_FONT,hjust=1) +
    annotate("text",x=-x_max*0.85,y=y_cap*0.96,
             label=paste0("Down: ",n_dn),
             colour=COL_SIG,size=3.0,family=BASE_FONT,hjust=0) +
    scale_colour_manual(values=col_map,name=NULL) +
    scale_size_manual(  values=size_map,name=NULL) +
    scale_alpha_manual( values=alpha_map,name=NULL) +
    scale_y_continuous(limits=c(0,y_cap)) +
    labs(
      title    = "(A)  TMT 6-plex proteomics: RA vs OA synovial tissue",
      subtitle = sprintf(
        "PRIDE PXD027703 | n = 2 RA pools vs 3 OA pools (ch5 excluded, pooled ref)\n|log\u2082FC| \u2265 %.1f, adj.P < %.1f (exploratory) | Labelled = MR candidate hub genes",
        fc_thr, pval_thr),
      x = expression(log[2]~"fold-change  (RA / OA)"),
      y = expression(-log[10](adj.P))
    ) +
    theme_pub() +
    theme(legend.position="right", legend.key.size=unit(0.35,"cm"))
}

p9A <- make_9A(limma_d, tgt_cmp_d, FC_THR, PVAL_THR)
save_fig(p9A, file.path(out_dir,"Figure9A_volcano"), 10, 7)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 9B — Direction Consistency (deduplicated, NO pivot_longer)
# ══════════════════════════════════════════════════════════════════════════════
wlog("--- 9B ---")

make_9B <- function(dir_cons_d, tgt_cmp_d) {
  
  dc <- if (nrow(dir_cons_d)>0 && "Gene" %in% names(dir_cons_d)) {
    as.data.table(dir_cons_d)
  } else if (nrow(tgt_cmp_d)>0 && "Gene" %in% names(tgt_cmp_d)) {
    wlog("9B: using tgt_cmp_d as fallback")
    as.data.table(tgt_cmp_d)
  } else {
    wlog("9B: no data")
    return(ggplot()+theme_void()+
             annotate("text",x=0.5,y=0.5,label="No direction data",
                      family=BASE_FONT,size=5))
  }
  
  # Standardise logFC
  lfc_col <- intersect(c("logFC","log2FC","LogFC"),names(dc))[1]
  if (!is.na(lfc_col) && lfc_col!="logFC") { dc<-copy(dc); setnames(dc,lfc_col,"logFC") }
  
  # Ensure consistency columns exist
  if (!"Single_consistent" %in% names(dc)) dc[, Single_consistent := NA]
  if (!"Combo_consistent"  %in% names(dc)) dc[, Combo_consistent  := NA]
  if (!"logFC"             %in% names(dc)) dc[, logFC             := NA_real_]
  
  keep <- intersect(c("Gene","logFC","adj.P.Val","Single_consistent",
                      "Combo_consistent","Module"), names(dc))
  rows <- unique(dc[, ..keep])
  
  # Final dedup: one row per Gene
  if (any(duplicated(rows$Gene))) rows <- rows[!duplicated(Gene)]
  
  wlog("9B: ",nrow(rows)," unique genes | S_cons non-NA=",
       sum(!is.na(rows$Single_consistent)),
       " | C_cons non-NA=",sum(!is.na(rows$Combo_consistent)))
  
  status_lbl <- function(x)
    ifelse(is.na(x), "No MR \u03b2",
           ifelse(x==TRUE, "Consistent", "Inconsistent"))
  
  rows[, single_status := status_lbl(Single_consistent)]
  rows[, combo_status  := status_lbl(Combo_consistent)]
  rows[, logFC_val     := ifelse(is.na(logFC), 0, as.numeric(logFC))]
  
  # Gene order
  g_ord <- c("SLC39A8","LAMP3","PRF1","GZMB",
             "TUBA1A","TUBA1B","TUBA1C","TUBA3D",
             "TUBA3E","TUBA4A","TUBA8","TUBAL3",
             "NKG7","CLEC10A","HLA-F")
  all_g  <- unique(rows$Gene)
  pres   <- intersect(g_ord, all_g)
  extra  <- setdiff(all_g, g_ord)
  g_lvls <- rev(c(pres, extra))
  
  # Long format via rbind (NO pivot_longer)
  single_df <- data.table(Gene=rows$Gene, MR_type="Single-gene MR",
                          status=rows$single_status, logFC_val=rows$logFC_val)
  combo_df  <- data.table(Gene=rows$Gene, MR_type="Combination MR",
                          status=rows$combo_status,  logFC_val=rows$logFC_val)
  long_df <- rbind(single_df, combo_df)
  long_df[, Gene    := factor(Gene,    levels=g_lvls)]
  long_df[, MR_type := factor(MR_type, levels=c("Single-gene MR","Combination MR"))]
  long_df[, status  := factor(status,
                              levels=c("Consistent","Inconsistent","No MR \u03b2"))]
  
  col_st <- c("Consistent"=COL_CONSIST,"Inconsistent"=COL_INCONS,
              "No MR \u03b2"=COL_NOMB)
  
  lfc_ann <- unique(long_df[,.(Gene,logFC_val)])
  
  ggplot(long_df, aes(x=Gene,y=1,fill=status)) +
    geom_col(colour="white",linewidth=0.3,width=0.72) +
    coord_flip() +
    facet_wrap(~MR_type, ncol=2) +
    geom_text(data=lfc_ann,
              aes(x=Gene,y=0.5,label=sprintf("%+.2f",logFC_val),fill=NULL),
              size=2.8,family=BASE_FONT,colour="white",fontface="bold") +
    scale_fill_manual(values=col_st,name="Direction\nconsistency",drop=FALSE) +
    scale_y_continuous(breaks=NULL,labels=NULL,expand=expansion(mult=c(0,0))) +
    labs(
      title    = "(B)  MR \u03b2 vs TMT log\u2082FC direction consistency",
      subtitle = "Bar text = TMT log\u2082FC (RA / OA) | Green = MR and proteomics directions agree\nControl = OA synovial tissue",
      x=NULL, y=NULL
    ) +
    theme_pub() +
    theme(legend.position="right",
          panel.grid=element_blank(),
          axis.line=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text=element_text(face="bold"))
}

p9B <- make_9B(dir_cons_d, tgt_cmp_d)
save_fig(p9B, file.path(out_dir,"Figure9B_direction"),
         10, max(5, nrow(dir_cons_d)*0.42+3))

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 9C — Heatmap
# [D1] tmt_norm: ProteinID + Description + TMT channels (no Gene column)
#      → parse Gene from "GN=XXXX" in Description
#      → merge with det for Gene→ProteinID link
#      → one row per gene: pick row with highest mean intensity
# ══════════════════════════════════════════════════════════════════════════════
wlog("--- 9C ---")

make_9C <- function(tmt_norm, tgt_cmp_d, det, sample_info,
                    pdf_path, png_path, tiff_path) {
  
  if (nrow(tmt_norm)==0) { wlog("9C: no tmt_norm"); return(invisible(NULL)) }
  
  HUB <- c("GZMB","TUBA1A","TUBA1B","TUBA1C","TUBA3D","TUBA3E",
           "TUBA4A","TUBA8","TUBAL3","LAMP3","NKG7","PRF1",
           "SLC39A8","CLEC10A","HLA-F")
  
  # ── Channel columns (exclude ch5 = pooled reference outlier) ────────────────
  tmt_cols <- intersect(sample_info$Channel, names(tmt_norm))
  if (length(tmt_cols)==0)
    tmt_cols <- grep("^TMT_ch[^5]|^TMT_ch[1-4]$|^TMT_ch6$",
                     names(tmt_norm), value=TRUE)
  si_use   <- sample_info[sample_info$Channel %in% tmt_cols, ]
  tmt_cols <- si_use$Channel
  wlog("9C: tmt_cols=", paste(tmt_cols, collapse=","))
  
  # ── [H1] Parse Gene from Description ────────────────────────────────────────
  tm <- copy(tmt_norm)
  if (!"Gene" %in% names(tm)) {
    if (!"Description" %in% names(tm)) {
      wlog("9C: no Gene or Description column — skipped")
      return(invisible(NULL))
    }
    # Primary: GN=XXXX pattern
    tm[, Gene := str_match(Description, "GN=([^\\s;|>]+)")[, 2]]
    # Fallback: match hub gene names directly in Description string
    no_gene <- is.na(tm$Gene)
    if (any(no_gene)) {
      for (g in HUB) {
        pattern <- paste0("(?i)(^|\\s|\\|)", g, "(\\s|_|\\||$)")
        hits <- no_gene & grepl(pattern, tm$Description, perl=TRUE)
        if (any(hits)) tm[hits, Gene := g]
      }
    }
    wlog("9C: Gene parsed — non-NA=", sum(!is.na(tm$Gene)),
         " / ", nrow(tm))
  }
  
  # ── Filter to hub genes ──────────────────────────────────────────────────────
  tm_hub <- tm[!is.na(Gene) & Gene %in% HUB]
  wlog("9C: hub rows=", nrow(tm_hub),
       " | unique genes=", length(unique(tm_hub$Gene)))
  
  if (nrow(tm_hub)==0) {
    wlog("9C: no hub gene rows found — skipped")
    return(invisible(NULL))
  }
  
  # ── [H2] Impute NA: half-minimum per row ────────────────────────────────────
  # [H2] Impute NA/NaN/Inf using matrix operations (avoids data.table index recycling bug)
  mat_imp <- as.matrix(tm_hub[, ..tmt_cols])
  storage.mode(mat_imp) <- "numeric"
  mat_imp[!is.finite(mat_imp)] <- NA_real_   # unify all bad values to NA
  
  # Pass 1: half-row-minimum imputation
  row_min_vec <- apply(mat_imp, 1, function(x) {
    v <- x[is.finite(x) & !is.na(x)]
    if (length(v) > 0) min(v) else NA_real_
  })
  for (i in seq_len(nrow(mat_imp))) {
    bad <- is.na(mat_imp[i, ])
    if (any(bad) && is.finite(row_min_vec[i]))
      mat_imp[i, bad] <- row_min_vec[i] / 2
  }
  
  # Pass 2: row-mean imputation for any remaining NA
  for (i in seq_len(nrow(mat_imp))) {
    bad <- is.na(mat_imp[i, ])
    if (any(bad)) {
      good_mean <- mean(mat_imp[i, !bad], na.rm = TRUE)
      mat_imp[i, bad] <- if (is.finite(good_mean)) good_mean else 0
    }
  }
  
  # Write imputed values back to tm_hub columns
  for (j in seq_along(tmt_cols))
    set(tm_hub, j = tmt_cols[j], value = mat_imp[, j])
  
  # ── One row per gene: highest mean intensity ─────────────────────────────────
  tm_hub[, mean_int := rowMeans(.SD, na.rm=TRUE), .SDcols=tmt_cols]
  tm_hub <- tm_hub[order(-mean_int)][!duplicated(Gene)]
  wlog("9C: after dedup=", nrow(tm_hub), " genes")
  
  mat <- as.matrix(tm_hub[, ..tmt_cols])
  rownames(mat) <- tm_hub$Gene
  storage.mode(mat) <- "numeric"
  
  # ── [H3] Remove constant rows (sd=0 → NaN in correlation distance) ──────────
  row_sd <- apply(mat, 1, sd, na.rm=TRUE)
  const_rows <- which(row_sd < 1e-10 | !is.finite(row_sd))
  if (length(const_rows) > 0) {
    wlog("9C: removing ", length(const_rows),
         " constant rows: ", paste(rownames(mat)[const_rows], collapse=","))
    mat <- mat[-const_rows, , drop=FALSE]
  }
  
  # Rename cols to SampleID
  colnames(mat) <- si_use$SampleID[match(colnames(mat), si_use$Channel)]
  
  # ── [H4] Final NA/Inf safety sweep ──────────────────────────────────────────
  mat[!is.finite(mat)] <- NA_real_
  if (any(is.na(mat))) {
    # Column mean imputation as last resort
    col_means <- colMeans(mat, na.rm=TRUE)
    for (j in seq_len(ncol(mat)))
      mat[is.na(mat[,j]), j] <- col_means[j]
  }
  
  wlog("9C: final matrix dim=", nrow(mat), "×", ncol(mat),
       " | NA=", sum(is.na(mat)),
       " | Inf=", sum(!is.finite(mat)))
  
  if (nrow(mat)<2 || ncol(mat)<2) {
    wlog("9C: matrix too small after cleaning — skipped")
    return(invisible(NULL))
  }
  
  # ── Row annotation ────────────────────────────────────────────────────────────
  row_anno <- data.frame(row.names=rownames(mat))
  
  if (nrow(det)>0 && all(c("Gene","Module") %in% names(det))) {
    mod_lut        <- setNames(det$Module, det$Gene)
    row_anno$Module <- mod_lut[rownames(mat)]
    row_anno$Module[is.na(row_anno$Module)] <- "External"
  }
  
  if (nrow(tgt_cmp_d)>0 && "Gene" %in% names(tgt_cmp_d)) {
    dir_lut         <- setNames(tgt_cmp_d$TMT_direction, tgt_cmp_d$Gene)
    sig_lut         <- setNames(as.logical(tgt_cmp_d$Significant), tgt_cmp_d$Gene)
    row_anno$TMT_Dir <- dir_lut[rownames(mat)]
    row_anno$MR_Sig  <- ifelse(isTRUE(sig_lut[rownames(mat)]), "Yes", "No")
    row_anno$MR_Sig[is.na(row_anno$MR_Sig)] <- "No"
  }
  
  # Row labels with * for adj.P < 0.1
  if (nrow(tgt_cmp_d)>0 && "adj.P.Val" %in% names(tgt_cmp_d)) {
    star_lut <- setNames(
      ifelse(!is.na(tgt_cmp_d$adj.P.Val) & tgt_cmp_d$adj.P.Val < 0.1, " *", ""),
      tgt_cmp_d$Gene)
    row_labels <- paste0(rownames(mat),
                         ifelse(rownames(mat) %in% names(star_lut),
                                star_lut[rownames(mat)], ""))
  } else {
    row_labels <- rownames(mat)
  }
  
  # Column annotation
  col_anno <- data.frame(
    Group=si_use$Group[match(colnames(mat), si_use$SampleID)],
    row.names=colnames(mat)
  )
  
  # Annotation colours
  ann_colors <- list(Group=c(RA=COL_MR_SIG, OA=COL_SIG))
  
  if ("Module" %in% names(row_anno)) {
    mod_pal <- c(GZMB="#4A7C2F", LAMP3="#1A5276", NKG7="#784212",
                 TRAF3IP3="#7B241C", External="#707B7C")
    mods_here  <- unique(na.omit(row_anno$Module))
    extra_mods <- setdiff(mods_here, names(mod_pal))
    if (length(extra_mods)>0)
      mod_pal <- c(mod_pal,
                   setNames(colorRampPalette(c("#A9CCE3","#85C1E9"))(length(extra_mods)),
                            extra_mods))
    ann_colors$Module <- mod_pal[mods_here]
  }
  if ("TMT_Dir" %in% names(row_anno))
    ann_colors$TMT_Dir <- c(Up_in_RA=COL_MR_SIG, Down_in_RA=COL_SIG, NS="grey70")
  if ("MR_Sig" %in% names(row_anno))
    ann_colors$MR_Sig  <- c(Yes="#F39C12", No="grey85")
  
  h_in <- max(5, nrow(mat)*0.50 + 3)
  w_in <- 9
  
  # Use euclidean for both row and col when correlation may be unstable
  # (small n=5 samples → euclidean safer for col; correlation OK for rows
  #  after cleaning)
  row_dist <- tryCatch(
    "correlation",
    error=function(e) "euclidean"
  )
  
  draw_hmap <- function() {
    pheatmap(
      mat,
      scale                    = "row",
      clustering_distance_rows = "euclidean",   # safer than correlation for n=15
      clustering_distance_cols = "euclidean",
      clustering_method        = "ward.D2",
      color    = colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
      annotation_col    = col_anno,
      annotation_row    = row_anno,
      annotation_colors = ann_colors,
      labels_row        = row_labels,
      fontsize_row      = 10,
      fontsize_col      = 10,
      fontfamily        = BASE_FONT,
      border_color      = NA,
      main = paste0(
        "(C)  Normalised TMT intensities: 15 hub genes (RA vs OA)\n",
        "Row-scaled log\u2082 intensity | * adj.P < 0.1 (exploratory)"
      )
    )
  }
  
  for (fmt in c("pdf","png","tiff")) {
    fpath <- switch(fmt, pdf=pdf_path, png=png_path, tiff=tiff_path)
    tryCatch({
      if (fmt=="pdf")  pdf( fpath, width=w_in, height=h_in)
      if (fmt=="png")  png( fpath, width=w_in, height=h_in, units="in", res=300)
      if (fmt=="tiff") tiff(fpath, width=w_in, height=h_in, units="in",
                            res=600, compression="lzw")
      draw_hmap()
      dev.off()
      wlog("9C ", toupper(fmt), " saved: ", basename(fpath))
    }, error=function(e) {
      tryCatch(dev.off(), error=function(x) NULL)
      wlog("9C ", toupper(fmt), " FAILED: ", conditionMessage(e))
    })
  }
}

make_9C(
  tmt_norm, tgt_cmp_d, det, sample_info,
  pdf_path  = file.path(out_dir,"Figure9C_heatmap.pdf"),
  png_path  = file.path(out_dir,"Figure9C_heatmap_300.png"),
  tiff_path = file.path(out_dir,"Figure9C_heatmap_600.tiff")
)

# ══════════════════════════════════════════════════════════════════════════════
# Combined Figure 9
# ══════════════════════════════════════════════════════════════════════════════
wlog("--- Combined ---")
tryCatch({
  p9C_wrap <- if (file.exists(file.path(out_dir,"Figure9C_heatmap_300.png")) &&
                  requireNamespace("png",quietly=TRUE)) {
    img <- png::readPNG(file.path(out_dir,"Figure9C_heatmap_300.png"))
    ggplot() +
      annotation_raster(img,xmin=0,xmax=1,ymin=0,ymax=1) +
      coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand=FALSE) +
      theme_void()
  } else {
    ggplot()+theme_void()+
      annotate("text",x=0.5,y=0.5,family=BASE_FONT,size=4,
               label="Panel C: see Figure9C_heatmap_300.png")
  }
  
  fig9 <- (p9A / p9B) | p9C_wrap
  save_fig(fig9, file.path(out_dir,"Figure9_combined"), w_in=20, h_in=14)
  wlog("Combined saved")
}, error=function(e) wlog("Combined failed: ",conditionMessage(e)))

wlog("=== 703 v3 completed === outputs: ", out_dir)
cat("\n=== 703 v3 completed ===\nOutputs:", out_dir, "\n")