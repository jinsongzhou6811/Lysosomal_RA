# ============================================================
# 609_Figure8_Table2_Publication.R (v2 — font + layout fixes)
#
# Figure 8:
# (A) Top-10 GZMB weighted-combination MR forest plot
# (B) Score improvement bar chart (combo vs single-gene)
# (C) Single-gene forest: LAMP3, TUBA4A, GZMB
#
# Table 2: publication MR summary → CSV + DOCX
#
# Inputs (project root relative):
# outputs/06_UVMR_Analysis/05_Module_Screening_Results/
# GZMB/Final_Results/module_GZMB_final_top_combinations.csv
# GZMB/Final_Results/module_GZMB_single_vs_combo.csv
# outputs/06_UVMR_Analysis/04_UVMR_Data_filtering/
# SingleGene_AllPassing.csv
#
# Outputs → outputs/06_UVMR_Analysis/07_Figure8_Table2/
# ============================================================
suppressPackageStartupMessages({
  library(rprojroot)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork) # panel assembly — no cowplot height/width issue
  library(scales)
})
# ── [F1] Font: try Arial, fall back to "sans" ─────────────────────────────────
BASE_FONT <- tryCatch({
  if (requireNamespace("extrafont", quietly = TRUE)) {
    extrafont::loadfonts(device = "pdf", quiet = TRUE)
    extrafont::loadfonts(device = "postscript", quiet = TRUE)
    if ("Arial" %in% extrafont::fonts()) { message("Using Arial"); "Arial" }
    else { message("Arial not found, using sans"); "sans" }
  } else { message("extrafont not installed, using sans"); "sans" }
}, error = function(e) "sans")

project_root <- find_rstudio_root_file()

out_dir <- file.path(project_root, "outputs", "06_UVMR_Analysis",
                     "07_Figure8_Table2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "figure8_table2_log.txt")
log_con <- file(log_file, open = "wt", encoding = "UTF-8")
on.exit(close(log_con), add = TRUE)

wlog <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., sep=""))
  cat(msg, "\n"); writeLines(msg, log_con)
}

wlog("=== 610 Figure8 v2 started === font=", BASE_FONT)

safe_fread <- function(p) {
  if (!file.exists(p)) { wlog("MISSING: ", p); return(data.table()) }
  tryCatch(fread(p), error=function(e){ wlog("READ ERR: ",p); data.table() })
}

# [F3] save helper with explicit units
save_fig <- function(p, stem, w_in, h_in) {
  ggsave(paste0(stem, ".pdf"), plot=p, width=w_in, height=h_in,
         units="in", device="pdf", useDingbats=FALSE)
  ggsave(paste0(stem, "_300.png"), plot=p, width=w_in, height=h_in,
         units="in", device="png", dpi=300)
  ggsave(paste0(stem, "_600.tiff"), plot=p, width=w_in, height=h_in,
         units="in", device="tiff", dpi=600, compression="lzw")
  wlog("Saved: ", basename(stem))
}

# [F1] theme without base_family argument
theme_pub <- function(bs=11) {
  theme_classic(base_size=bs) +
    theme(
      text = element_text(family=BASE_FONT),
      plot.title = element_text(face="bold", size=bs+1, family=BASE_FONT),
      plot.subtitle = element_text(size=bs-2, colour="grey40", family=BASE_FONT),
      axis.text = element_text(size=bs-2, family=BASE_FONT),
      strip.background = element_rect(fill="grey92", colour=NA),
      strip.text = element_text(face="bold", size=bs-1, family=BASE_FONT),
      legend.text = element_text(size=bs-2, family=BASE_FONT),
      legend.title = element_text(size=bs-1, family=BASE_FONT),
      panel.grid.major.x = element_line(colour="grey92", linewidth=0.3)
    )
}

COL_RISK <- "#C0392B"
COL_PROT <- "#2980B9"
COL_COMBO <- "#E67E22"

# ── Input files ───────────────────────────────────────────────────────────────
res_dir <- file.path(project_root,"outputs","06_UVMR_Analysis",
                     "05_Module_Screening_Results")
filter_dir <- file.path(project_root,"outputs","06_UVMR_Analysis",
                        "04_UVMR_Data_filtering")

combo_raw <- safe_fread(file.path(res_dir,"GZMB","Final_Results",
                                  "module_GZMB_final_top_combinations.csv"))
svc_raw <- safe_fread(file.path(res_dir,"GZMB","Final_Results",
                                "module_GZMB_single_vs_combo.csv"))
sg_all <- safe_fread(file.path(filter_dir,"SingleGene_AllPassing.csv"))

wlog("combo=",nrow(combo_raw)," svc=",nrow(svc_raw)," sg_all=",nrow(sg_all))

# ══════════════════════════════════════════════════════════════════════════════
# 8A — Top-10 Combination Forest Plot
# ══════════════════════════════════════════════════════════════════════════════
wlog("--- 8A ---")
if (nrow(combo_raw) > 0) {
  c10 <- combo_raw %>%
    arrange(desc(composite_score)) %>% slice_head(n=10) %>%
    mutate(
      rank = row_number(),
      lab = factor(
        sprintf("#%d %s", rank, str_replace_all(gene,"_","+")),
        levels = rev(sprintf("#%d %s", rank,
                             str_replace_all(gene,"_","+")))),
      OR = exp(ivw_b),
      OR_lo = exp(ivw_b - 1.96*ivw_se),
      OR_hi = exp(ivw_b + 1.96*ivw_se),
      sig = case_when(ivw_pval_BH<0.001~"***",
                      ivw_pval_BH<0.01 ~"**",
                      ivw_pval_BH<0.05 ~"*", TRUE~""),
      annot = sprintf("nSNP=%d F=%.0f p=%.1e%s",
                      nsnp, mean_F, ivw_pval_BH, sig),
      direction = ifelse(ivw_b>0,"Risk","Protective")
    )
  
  xlo <- max(floor(min(c10$OR_lo)*500)/500, 0.990)
  xhi <- min(ceiling(max(c10$OR_hi)*500)/500, 1.040)
  
  p8A <- ggplot(c10, aes(x=OR, y=lab, colour=direction)) +
    geom_vline(xintercept=1, linetype="dashed", colour="grey40", linewidth=0.5) +
    geom_errorbarh(aes(xmin=OR_lo, xmax=OR_hi), height=0.3, linewidth=0.6) +
    geom_point(size=4, shape=18) +
    geom_text(aes(x=xhi+0.002, label=annot), hjust=0, size=2.7,
              colour="grey25", family=BASE_FONT) +
    scale_colour_manual(values=c(Risk=COL_RISK, Protective=COL_PROT),
                        name="Direction") +
    scale_x_continuous(limits=c(xlo-0.002, xhi+0.028),
                       breaks=pretty(c(xlo,xhi),n=6),
                       labels=function(x) sprintf("%.3f",x)) +
    labs(title="(A) GZMB module: Top-10 weighted-combination MR",
         subtitle="Outcome: ieu-a-832 (Okada 2014, European RA)\negger_flag=FALSE, direction_ratio=1.0 for all combinations",
         x="Odds Ratio (95% CI)", y=NULL) +
    theme_pub() +
    theme(legend.position="none", plot.margin=margin(5,90,5,5,"pt"))
  
  save_fig(p8A, file.path(out_dir,"Figure8A_combo_forest"), 12, 6)
} else {
  p8A <- ggplot()+theme_void()+
    annotate("text",x=0.5,y=0.5,label="No combo data",family=BASE_FONT,size=5)
  wlog("8A: no data")
}

# ══════════════════════════════════════════════════════════════════════════════
# 8B — Score Improvement
# ══════════════════════════════════════════════════════════════════════════════
wlog("--- 8B ---")
if (nrow(svc_raw) > 0) {
  sb <- svc_raw %>%
    filter(!is.na(Score_Improvement)) %>%
    mutate(
      combo_dir = ifelse(Combo_b>0,"Risk","Protective"),
      passed_606 = as.character(Single_passed_606),
      gene = factor(Single_gene,
                    levels=Single_gene[order(Score_Improvement)])
    )
  
  p8B <- ggplot(sb, aes(x=Score_Improvement, y=gene,
                        fill=combo_dir, alpha=passed_606)) +
    geom_col(colour="white", linewidth=0.3, width=0.70) +
    geom_vline(xintercept=0, linewidth=0.5, colour="grey30") +
    geom_text(aes(label=sprintf("%+.1f", Score_Improvement)),
              hjust=ifelse(sb$Score_Improvement>=0,-0.2,1.2),
              size=3, family=BASE_FONT, colour="grey20") +
    scale_fill_manual(name="Combo direction",
                      values=c(Risk=COL_RISK,Protective=COL_PROT)) +
    scale_alpha_manual(name="Single-gene sig.",
                       values=c("TRUE"=1.0,"FALSE"=0.55),
                       labels=c("TRUE"="Yes","FALSE"="No")) +
    scale_x_continuous(expand=expansion(mult=c(0.2,0.22))) +
    labs(title="(B) Score improvement: combination vs single-gene MR",
         subtitle="\u0394Score = combo composite score \u2212 single-gene composite score",
         x="\u0394 Composite score", y=NULL) +
    theme_pub() +
    theme(legend.position="right")
  
  save_fig(p8B, file.path(out_dir,"Figure8B_score_improvement"),
           10, max(4.5, nrow(sb)*0.55+2.5))
} else {
  p8B <- ggplot()+theme_void()+
    annotate("text",x=0.5,y=0.5,label="No svc data",family=BASE_FONT,size=5)
  wlog("8B: no data")
}

# ══════════════════════════════════════════════════════════════════════════════
# 8C — Single-Gene Forest (hardcoded from MR report)
# ══════════════════════════════════════════════════════════════════════════════
wlog("--- 8C ---")
sg <- tribble(
  ~gene, ~outcome, ~pop_cat, ~OR, ~OR_lo, ~OR_hi, ~pBH, ~soft,
  "LAMP3", "ieu-a-831", "European", 1.103, 1.057, 1.151, 5.9e-5, FALSE,
  "LAMP3", "ieu-a-833", "European", 1.057, 1.029, 1.086, 2.3e-4, FALSE,
  "LAMP3", "bbj-a-151", "East Asian", 1.049, 1.013, 1.086, 1.6e-2, FALSE,
  "LAMP3", "bbj-a-74", "East Asian", 1.039, 1.010, 1.069, 1.6e-2, FALSE,
  "TUBA4A", "bbj-a-151", "East Asian", 0.884, 0.844, 0.926, 1.9e-6, FALSE,
  "TUBA4A", "ieu-a-833", "European", 0.910, 0.876, 0.944, 2.8e-6, FALSE,
  "TUBA4A", "bbj-a-74", "East Asian", 0.913, 0.881, 0.947, 2.8e-6, FALSE,
  "TUBA4A", "ieu-a-831", "European", 0.887, 0.841, 0.936, 3.1e-5, FALSE,
  "TUBA4A", "bbj-a-72", "East Asian", 0.891, 0.845, 0.939, 3.4e-5, FALSE,
  "TUBA4A", "ieu-a-832", "European", 0.931, 0.885, 0.979, 7.5e-3, FALSE,
  "TUBA4A", "bbj-a-73", "East Asian", 0.933, 0.888, 0.979, 7.5e-3, FALSE,
  "GZMB", "ukb-d-M06", "European", 0.9998, 0.9996, 0.9999, 8.1e-3, TRUE
) %>%
  mutate(
    sig = case_when(pBH<0.001~"***", pBH<0.01~"**",
                    pBH<0.05~"*", TRUE~""),
    direction = ifelse(OR>=1, "Risk", "Protective"),
    flag_tag = ifelse(soft, " \u2020", ""),
    y_lab = paste0(outcome, flag_tag, " ", sig),
    gene_fac = factor(
      case_when(
        gene=="LAMP3" ~ "LAMP3\n(Risk \u2191, 4 outcomes, EU+EA)",
        gene=="TUBA4A" ~ "TUBA4A\n(Protective \u2193, 7 outcomes, EU+EA)",
        gene=="GZMB" ~ "GZMB\n(Minimal, 1 outcome, EU ICD-10 \u2020)"
      ),
      levels = c("LAMP3\n(Risk \u2191, 4 outcomes, EU+EA)",
                 "TUBA4A\n(Protective \u2193, 7 outcomes, EU+EA)",
                 "GZMB\n(Minimal, 1 outcome, EU ICD-10 \u2020)")
    ),
    y_ord = reorder(paste0(gene,"__",y_lab), -pBH)
  )

or_lo_g <- min(sg$OR_lo) - 0.01
or_hi_g <- max(sg$OR_hi) + 0.01

p8C <- ggplot(sg, aes(x=OR, y=y_ord, colour=direction, shape=pop_cat)) +
  geom_vline(xintercept=1, linetype="dashed", colour="grey40", linewidth=0.5) +
  geom_errorbarh(aes(xmin=OR_lo, xmax=OR_hi), height=0.30, linewidth=0.55) +
  geom_point(size=3.2) +
  facet_wrap(~gene_fac, ncol=1, scales="free_y") +
  scale_colour_manual(name="Direction",
                      values=c(Risk=COL_RISK, Protective=COL_PROT)) +
  scale_shape_manual(name="Population",
                     values=c("European"=16, "East Asian"=17)) +
  scale_y_discrete(labels=function(x) sub("^[^_]+__","",x)) +
  scale_x_continuous(limits=c(or_lo_g, or_hi_g+0.005),
                     breaks=pretty(c(or_lo_g, or_hi_g), n=6),
                     labels=function(x) sprintf("%.3f",x)) +
  labs(title="(C) Single-gene UVMR: three significant genes",
       subtitle="\u2020 GZMB: PRESSO significant, direction_ratio=0.75 (soft QC flags only)\nAll other entries: egger_flag=FALSE, direction_ratio=1.0",
       x="Odds Ratio (95% CI)", y=NULL) +
  theme_pub() +
  theme(legend.position="right", panel.spacing=unit(0.5,"cm"))

save_fig(p8C, file.path(out_dir,"Figure8C_singlegene_forest"), 10, 10)
wlog("8C saved")

# ══════════════════════════════════════════════════════════════════════════════
# Combined Figure 8 [F2] patchwork
# ══════════════════════════════════════════════════════════════════════════════
wlog("--- Combined Figure 8 ---")
fig8 <- (p8A / p8B) | p8C
save_fig(fig8, file.path(out_dir,"Figure8_combined"), w_in=22, h_in=14)
wlog("Combined Figure8 saved")

# ══════════════════════════════════════════════════════════════════════════════
# Table 2 — CSV + DOCX
# ══════════════════════════════════════════════════════════════════════════════
wlog("--- Table 2 ---")
t2 <- tribble(
  ~Section, ~Gene, ~Outcome, ~Population, ~`OR (95% CI)`,
  ~`pBH`, ~nSNP, ~`Mean F`, ~Direction, ~`QC notes`,
  
  "Single gene","LAMP3","ieu-a-831","European (Okada)",
  "1.103 (1.057-1.151)","5.9e-5",13,67.9,"Risk","All pass",
  
  "Single gene","LAMP3","ieu-a-833","European (Okada)",
  "1.057 (1.029-1.086)","2.3e-4",13,67.9,"Risk","All pass",
  
  "Single gene","LAMP3","bbj-a-151","East Asian (BBJ)",
  "1.049 (1.013-1.086)","1.6e-2",13,67.9,"Risk","All pass",
  
  "Single gene","LAMP3","bbj-a-74","East Asian (BBJ)",
  "1.039 (1.010-1.069)","1.6e-2",13,67.9,"Risk","All pass",
  
  "Single gene","TUBA4A","bbj-a-151","East Asian (BBJ)",
  "0.884 (0.844-0.926)","1.9e-6",11,27.8,"Protective","All pass",
  
  "Single gene","TUBA4A","ieu-a-833","European (Okada)",
  "0.910 (0.876-0.944)","2.8e-6",11,27.8,"Protective","All pass",
  
  "Single gene","TUBA4A","bbj-a-74","East Asian (BBJ)",
  "0.913 (0.881-0.947)","2.8e-6",11,27.8,"Protective","All pass",
  
  "Single gene","TUBA4A","ieu-a-831","European (Okada)",
  "0.887 (0.841-0.936)","3.1e-5",11,27.8,"Protective","All pass",
  
  "Single gene","TUBA4A","bbj-a-72","East Asian (BBJ)",
  "0.891 (0.845-0.939)","3.4e-5",11,27.8,"Protective","All pass",
  
  "Single gene","TUBA4A","ieu-a-832","European (Okada)",
  "0.931 (0.885-0.979)","7.5e-3",11,27.8,"Protective","All pass",
  
  "Single gene","TUBA4A","bbj-a-73","East Asian (BBJ)",
  "0.933 (0.888-0.979)","7.5e-3",11,27.8,"Protective","All pass",
  
  "Single gene","GZMB","ukb-d-M06 (ICD-10)","European (UKB)",
  "0.9998 (0.9996-0.9999)","8.1e-3",24,60.8,"Minimal",
  "PRESSO sig+; dir=0.75+",
  
  "Weighted combination",
  "GZMB+TUBA1A+TUBA1B+TUBA3D+TUBA3E+TUBAL3",
  "ieu-a-832","European (Okada)",
  "1.017 (1.009-1.025)","1.8e-4",82,162.0,"Risk",
  "All pass; egger=FALSE; dir=1.0"
)

fwrite(t2, file.path(out_dir,"Table2_MR_Results.csv"))
wlog("Table 2 CSV saved")

# ── DOCX (optional) ───────────────────────────────────────────────────────────
tryCatch({
  if (!requireNamespace("flextable",quietly=TRUE) ||
      !requireNamespace("officer", quietly=TRUE))
    stop("flextable/officer not installed")
  library(flextable); library(officer)
  
  ft <- flextable(t2) %>%
    bold(part="header") %>%
    bg(i= ~Gene=="LAMP3", bg="#FFF3F3", part="body") %>%
    bg(i= ~Gene=="TUBA4A", bg="#F0F4FF", part="body") %>%
    bg(i= ~Section=="Weighted combination", bg="#FFF0E6", part="body") %>%
    hline(i=sum(t2$Section=="Single gene"),
          border=fp_border_default(color="#333",width=1.5), part="body") %>%
    font(fontname="Arial", part="all") %>%
    fontsize(size=9, part="body") %>%
    fontsize(size=10, part="header") %>%
    set_table_properties(layout="autofit") %>%
    add_header_lines(
      "Table 2. Statistically significant MR results for lysosomal-cytoskeletal hub genes") %>%
    add_footer_lines(paste0(
      "OR=odds ratio (IVW); pBH=Benjamini-Hochberg adjusted P; F=mean instrument F-statistic. ",
      "+ GZMB soft QC flags: PRESSO significant, direction_ratio=0.75. ",
      "EU=European; EA=East Asian.")) %>%
    italic(part="footer") %>% fontsize(size=8, part="footer")
  
  doc <- read_docx() %>%
    body_add_par("Table 2", style="heading 1") %>%
    body_add_flextable(ft)
  print(doc, target=file.path(out_dir,"Table2_MR_Results.docx"))
  wlog("Table 2 DOCX saved")
}, error=function(e) {
  wlog("DOCX skipped (flextable/officer missing or error): ", conditionMessage(e))
  wlog("Install: install.packages(c('flextable','officer'))")
})

wlog("=== 610 Figure8 v2 completed === outputs in: ", out_dir)
cat("\n=== 610 Figure8 v2 completed ===\nOutputs:", out_dir, "\n")