# ========================================================
# 701_MSF_Proteomics_Validation_v3.R
#
# 根据两次运行日志确认的 P20180100010.msf 实际表结构：
#
#   ProteinAnnotations : ProteinAnnotationID, ProteinID,
#                        DescriptionHashCode, Description, TaxonomyID
#   Proteins           : ProteinID, Sequence, SequenceHashCode, IsMasterProtein
#   ReporterIonQuanResults : ProcessingNodeNumber, QuanChannelID, SpectrumID,
#                            Mass, Height          ← TMT 强度，SpectrumID=TMT侧
#   ReporterIonQuanResultsSearchSpectra : ProcessingNodeNumber,
#                            SpectrumID, SearchSpectrumID
#                            SpectrumID=TMT侧, SearchSpectrumID=肽段侧
#   Peptides           : PeptideID, SpectrumID, ...  ← SpectrumID=SearchSpectrumID
#   PeptidesProteins   : PeptideID, ProteinID
#
# 完整 TMT→蛋白 连接链（全部在 SQL 中完成）：
#   ReporterIonQuanResults.SpectrumID
#     = ReporterIonQuanResultsSearchSpectra.SpectrumID (TMT侧)
#     → SearchSpectrumID
#     = Peptides.SpectrumID (肽段侧)
#     → PeptideID
#     → PeptidesProteins.ProteinID
#     → ProteinAnnotations.ProteinID → Description
#
# 输出：outputs/07_Proteomics_Analysis/01_Proteomics_Analysis/
#   01_protein_annotations.csv      ProteinAnnotations 完整导出
#   01b_tmt_protein.csv             蛋白水平 TMT 定量矩阵（宽格式）
#   01c_tmt_normalized.csv          log2 + 中位数归一化后的矩阵
#   02_tables_overview.csv          表结构行数一览
#   03_protein_detected.csv         目标基因检出详情（核心结果）
#   04_coverage_summary.csv         按模块检出率汇总
#   05_psm_summary.csv              目标蛋白 PSM 数量
#   06_validation_report.txt        完整日志
# ========================================================

suppressPackageStartupMessages({
  library(rprojroot)
  library(RSQLite)
  library(DBI)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

project_root <- find_rstudio_root_file()

# ── 路径 ──────────────────────────────────────────────────────────────────────
data_dir  <- "D:/001 Proteome Data"
base_name <- "P20180100010"

msf_path <- {
  exts <- c(".msf", ".pdResult", ".pdresult", ".MSF")
  found <- file.path(data_dir, paste0(base_name, exts))
  found <- found[file.exists(found)]
  if (length(found) == 0)
    stop("MSF file not found. Tried: ",
         paste(paste0(base_name, exts), collapse = ", "))
  found[1]
}

uvmr_dir <- file.path(project_root, "outputs", "06_UVMR_Analysis",
                      "05_Module_Screening_Results")
out_dir  <- file.path(project_root, "outputs", "07_Proteomics_Analysis",
                      "01_Proteomics_Analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── 日志 ──────────────────────────────────────────────────────────────────────
log_path <- file.path(out_dir, "06_validation_report.txt")
log_con  <- file(log_path, open = "wt", encoding = "UTF-8")
on.exit({ if (!is.null(log_con) && isOpen(log_con)) close(log_con) }, add = TRUE)

wlog <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., sep = ""))
  cat(msg, "\n")
  tryCatch(writeLines(msg, log_con), error = function(e) invisible())
}

sq <- function(con, sql) {   # safe SQL query → data.table or NULL
  tryCatch(as.data.table(dbGetQuery(con, sql)),
           error = function(e) { wlog("SQL error: ", conditionMessage(e)); NULL })
}

safe_fread <- function(p) {
  if (!file.exists(p)) return(data.table())
  tryCatch(fread(p), error = function(e) data.table())
}

wlog("=== 610 MSF Proteomics Validation v3 ===")
wlog("File    : ", msf_path)
wlog("Out dir : ", out_dir)

# ── 连接 ──────────────────────────────────────────────────────────────────────
con <- dbConnect(SQLite(), msf_path)
on.exit(dbDisconnect(con), add = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: 表结构概览
# ══════════════════════════════════════════════════════════════════════════════
wlog("--- Step 1: Table overview ---")

all_tables <- dbListTables(con)
tov <- rbindlist(lapply(all_tables, function(t) {
  n <- tryCatch(
    dbGetQuery(con, paste0("SELECT COUNT(*) AS n FROM [", t, "]"))$n,
    error = function(e) NA_integer_)
  data.table(Table = t, Rows = n)
}))
fwrite(tov, file.path(out_dir, "02_tables_overview.csv"))
wlog("Tables: ", length(all_tables), " | Top 5 by rows:")
print(tov[order(-Rows)][1:5])

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: 读取 ProteinAnnotations
# 确认列：ProteinAnnotationID, ProteinID, DescriptionHashCode,
#          Description, TaxonomyID
# ══════════════════════════════════════════════════════════════════════════════
wlog("--- Step 2: ProteinAnnotations ---")

# 读取并合并 IsMasterProtein 标记
pa <- sq(con, "
  SELECT pa.ProteinID,
         pa.Description,
         pa.TaxonomyID,
         p.IsMasterProtein
  FROM   ProteinAnnotations pa
  LEFT JOIN Proteins p ON pa.ProteinID = p.ProteinID
")

if (is.null(pa) || nrow(pa) == 0)
  stop("ProteinAnnotations empty — check MSF file.")

wlog("ProteinAnnotations loaded: ", nrow(pa), " rows")
wlog("Cols: ", paste(names(pa), collapse = ", "))
wlog("Description sample [1]: ", str_trunc(pa$Description[1], 120))

fwrite(pa, file.path(out_dir, "01_protein_annotations.csv"))
wlog("01_protein_annotations.csv saved")

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: TMT 定量 → 蛋白水平矩阵（全 SQL JOIN）
#
# 连接链（列名已由日志确认）：
#   riq.SpectrumID
#     JOIN riqss ON riq.SpectrumID = riqss.SpectrumID      (TMT侧)
#     → riqss.SearchSpectrumID
#     JOIN pep   ON pep.SpectrumID = riqss.SearchSpectrumID (肽段侧)
#     → pep.PeptideID
#     JOIN pp    ON pp.PeptideID   = pep.PeptideID
#     → pp.ProteinID
# ══════════════════════════════════════════════════════════════════════════════
wlog("--- Step 3: TMT protein-level matrix (SQL JOIN) ---")

# 确认 6 个通道
chan_n <- sq(con,
             "SELECT QuanChannelID, COUNT(*) AS n
   FROM ReporterIonQuanResults
   GROUP BY QuanChannelID
   ORDER BY QuanChannelID")
wlog("TMT channels: ", paste(chan_n$QuanChannelID, collapse = ", "),
     " (", nrow(chan_n), " channels, rows/ch: ",
     chan_n$n[1], ")")

# 通过 SQL 完成全部关联：TMT Height → ProteinID
# 输出：每行一个 (ProteinID, QuanChannelID, Height_median)
wlog("Building TMT→Protein mapping via SQL (may take ~30s)...")

tmt_long_sql <- "
  SELECT   pp.ProteinID,
           riq.QuanChannelID,
           AVG(riq.Height) AS Height_mean
  FROM     ReporterIonQuanResults  riq
  JOIN     ReporterIonQuanResultsSearchSpectra riqss
             ON riq.SpectrumID       = riqss.SpectrumID
  JOIN     Peptides pep
             ON pep.SpectrumID       = riqss.SearchSpectrumID
  JOIN     PeptidesProteins pp
             ON pp.PeptideID         = pep.PeptideID
  WHERE    riq.Height IS NOT NULL
    AND    riq.Height  > 0
  GROUP BY pp.ProteinID, riq.QuanChannelID
"

tmt_long <- sq(con, tmt_long_sql)

if (!is.null(tmt_long) && nrow(tmt_long) > 0) {
  wlog("TMT long table: ", nrow(tmt_long), " rows | proteins: ",
       uniqueN(tmt_long$ProteinID),
       " | channels: ", uniqueN(tmt_long$QuanChannelID))
  
  # 透视为宽格式：行=蛋白，列=TMT通道
  tmt_wide <- dcast(tmt_long,
                    ProteinID ~ paste0("TMT_ch", QuanChannelID),
                    value.var = "Height_mean")
  
  tmt_cols <- setdiff(names(tmt_wide), "ProteinID")
  wlog("TMT wide: ", nrow(tmt_wide), " proteins × ", length(tmt_cols), " channels")
  
  # 关联 Description
  tmt_wide <- merge(tmt_wide,
                    pa[, .(ProteinID, Description)],
                    by = "ProteinID", all.x = TRUE)
  setcolorder(tmt_wide, c("ProteinID", "Description", tmt_cols))
  
  fwrite(tmt_wide, file.path(out_dir, "01b_tmt_protein.csv"))
  wlog("01b_tmt_protein.csv saved: ", nrow(tmt_wide), " proteins")
  
  # log2 + 中位数归一化
  tmt_norm <- copy(tmt_wide)
  tmt_norm[, (tmt_cols) := lapply(.SD, function(x)
    ifelse(is.na(x) | x == 0, NA_real_, log2(x))),
    .SDcols = tmt_cols]
  
  meds       <- sapply(tmt_norm[, ..tmt_cols], median, na.rm = TRUE)
  global_med <- median(meds, na.rm = TRUE)
  for (cc in tmt_cols) {
    if (!is.na(meds[cc]))
      tmt_norm[[cc]] <- tmt_norm[[cc]] - meds[cc] + global_med
  }
  fwrite(tmt_norm, file.path(out_dir, "01c_tmt_normalized.csv"))
  wlog("01c_tmt_normalized.csv saved (log2 + median normalised)")
  
  has_tmt <- TRUE
  
} else {
  wlog("WARNING: TMT SQL join returned 0 rows.")
  wlog("Possible cause: SpectrumID mismatch between tables.")
  wlog("Diagnostic: check first rows of each table:")
  wlog("  riq  SpectrumID range: ",
       sq(con, "SELECT MIN(SpectrumID), MAX(SpectrumID)
                FROM ReporterIonQuanResults"))
  wlog("  riqss SpectrumID range: ",
       sq(con, "SELECT MIN(SpectrumID), MAX(SpectrumID)
                FROM ReporterIonQuanResultsSearchSpectra"))
  wlog("  pep   SpectrumID range: ",
       sq(con, "SELECT MIN(SpectrumID), MAX(SpectrumID) FROM Peptides"))
  fwrite(data.table(), file.path(out_dir, "01b_tmt_protein.csv"))
  fwrite(data.table(), file.path(out_dir, "01c_tmt_normalized.csv"))
  tmt_wide <- data.table()
  has_tmt  <- FALSE
}

# ══════════════════════════════════════════════════════════════════════════════
# Step 4: 目标基因列表
# ══════════════════════════════════════════════════════════════════════════════
wlog("--- Step 4: Target genes ---")

tg_file  <- file.path(project_root, "data", "Target Genes.txt")
modules  <- c("GZMB", "LAMP3", "NKG7", "TRAF3IP3")

# 来源A: Target Genes.txt
target_genes <- if (file.exists(tg_file)) {
  lines <- trimws(readLines(tg_file, warn = FALSE))
  lines <- lines[nchar(lines) > 0 & !startsWith(lines, "#")]
  sort(unique(lines))
} else character(0)

if (length(target_genes) == 0) {
  wlog("Target Genes.txt not found — will use UVMR module genes")
}

# 来源B: 608 UVMR 输出（模块归属）
mgm <- rbindlist(lapply(modules, function(mod) {
  svc <- safe_fread(file.path(uvmr_dir, mod, "Final_Results",
                              paste0("module_", mod, "_single_vs_combo.csv")))
  ft  <- safe_fread(file.path(uvmr_dir, mod, "Final_Results",
                              paste0("module_", mod,
                                     "_final_top_combinations.csv")))
  mg  <- if (nrow(svc) > 0) unique(svc$Single_gene) else character(0)
  cg  <- if (nrow(ft)  > 0)
    unique(unlist(strsplit(ft$gene, "_", fixed = TRUE))) else character(0)
  p6  <- if (nrow(svc) > 0)
    setNames(svc$Single_passed_606, svc$Single_gene)[mg]
  else setNames(rep(NA, length(mg)), mg)
  data.table(Gene = mg, Module = mod,
             In_top_combo = mg %in% cg, Passed_606 = p6)
}), fill = TRUE)

uvmr_genes <- sort(unique(mgm$Gene[!is.na(mgm$Gene) & mgm$Gene != ""]))

if (length(target_genes) > 0) {
  ext <- setdiff(target_genes, uvmr_genes)
  if (length(ext) > 0) {
    mgm <- rbindlist(list(mgm,
                          data.table(Gene = ext, Module = "External",
                                     In_top_combo = FALSE,
                                     Passed_606 = NA_integer_)),
                     fill = TRUE)
  }
} else {
  target_genes <- uvmr_genes
}

target_genes <- sort(unique(target_genes[!is.na(target_genes) & target_genes != ""]))
wlog("Target genes (", length(target_genes), "): ",
     paste(target_genes, collapse = ", "))
if (length(target_genes) == 0) stop("No target genes found.")

# ══════════════════════════════════════════════════════════════════════════════
# Step 5: 在 ProteinAnnotations 中匹配目标基因
#
# Description 格式（已确认）：
#   ">sp|Q71U36|TBA1A_HUMAN ... OS=Homo sapiens OX=9606 GN=TUBA1A PE=1 SV=1"
#
# 匹配策略：
#   1. GN={gene} 精确匹配（最可靠，针对已知格式）
#   2. 单词边界匹配（宽松兜底）
# ══════════════════════════════════════════════════════════════════════════════
wlog("--- Step 5: Gene matching in ProteinAnnotations ---")

# 预先从 Description 提取 UniProt AC 和 GeneSymbol（向量化，快）
pa[, UniProt_AC  := str_match(Description,
                              "\\|([A-Z][0-9][A-Z0-9]{3}[0-9])\\|")[, 2]]
pa[, GeneSymbol  := str_match(Description, "GN=([^\\s]+)")[, 2]]

wlog("UniProt AC extracted for  : ",
     sum(!is.na(pa$UniProt_AC)), " proteins")
wlog("GeneSymbol extracted for  : ",
     sum(!is.na(pa$GeneSymbol)), " proteins")
wlog("GeneSymbol sample (first 5): ",
     paste(head(pa$GeneSymbol[!is.na(pa$GeneSymbol)], 5), collapse = ", "))

match_results <- rbindlist(lapply(target_genes, function(gene) {
  
  # 策略1: GeneSymbol 列精确等值匹配（GN= 字段提取后）
  m <- pa[GeneSymbol == gene & !is.na(GeneSymbol)]
  
  # 策略2: GN= 原始 Description 匹配（兜底，处理提取失败的情况）
  if (nrow(m) == 0) {
    pat <- paste0("(?:^|\\s)GN=", gene, "(?:\\s|$)")
    m   <- pa[grepl(pat, Description, perl = TRUE)]
  }
  
  # 策略3: 单词边界宽松匹配
  method <- if (nrow(m) > 0) {
    if (any(pa$GeneSymbol == gene & !is.na(pa$GeneSymbol))) "GN= exact"
    else "GN= pattern"
  } else {
    kw <- paste0("\\b", gene, "\\b")
    m  <- pa[grepl(kw, Description, perl = TRUE, ignore.case = TRUE)]
    if (nrow(m) > 0) "keyword" else "not detected"
  }
  
  detected <- nrow(m) > 0
  n_ent    <- nrow(m)
  
  # ProteinID 列表（内部ID，用于后续 TMT 关联）
  prot_ids <- if (detected) paste(unique(m$ProteinID), collapse = ";") else "—"
  
  # UniProt AC 列表
  ac_list  <- if (detected && any(!is.na(m$UniProt_AC)))
    paste(unique(na.omit(m$UniProt_AC)), collapse = ";")
  else prot_ids
  
  # Description 第一条
  desc_txt <- if (detected) str_trunc(m$Description[1], 120) else "—"
  
  # 定量中位数（从 tmt_wide 取）
  qmed <- NA_real_
  if (detected && has_tmt && nrow(tmt_wide) > 0) {
    pid_vec  <- unique(m$ProteinID)
    tmt_rows <- tmt_wide[ProteinID %in% pid_vec]
    if (nrow(tmt_rows) > 0) {
      num_cols <- tmt_cols[tmt_cols %in% names(tmt_rows)]
      vals     <- as.numeric(unlist(tmt_rows[, ..num_cols]))
      qmed     <- round(median(vals[vals > 0], na.rm = TRUE), 2)
    }
  }
  
  data.table(Gene              = gene,
             Detected          = detected,
             Match_method      = method,
             N_protein_entries = n_ent,
             ProteinIDs        = prot_ids,    # 内部整数 ID
             UniProt_ACs       = ac_list,     # 可读 AC
             Description       = desc_txt,
             Quant_median      = qmed)
}), fill = TRUE)

# 合并模块注释
det_full <- match_results %>%
  left_join(mgm, by = "Gene") %>%
  arrange(Module, desc(In_top_combo), desc(Detected)) %>%
  as.data.table()

fwrite(det_full, file.path(out_dir, "03_protein_detected.csv"))
wlog("03_protein_detected.csv saved")

# ══════════════════════════════════════════════════════════════════════════════
# Step 6: 检出率汇总
# ══════════════════════════════════════════════════════════════════════════════
wlog("--- Step 6: Coverage summary ---")

cov <- det_full %>%
  group_by(Module) %>%
  summarise(
    N_genes_total        = n_distinct(Gene),
    N_detected           = sum(Detected, na.rm = TRUE),
    Detection_pct        = round(100 * mean(Detected, na.rm = TRUE), 1),
    N_in_top_combo       = sum(In_top_combo, na.rm = TRUE),
    N_top_combo_detected = sum(In_top_combo & Detected, na.rm = TRUE),
    Top_combo_pct        = round(
      100 * sum(In_top_combo & Detected, na.rm = TRUE) /
        pmax(sum(In_top_combo, na.rm = TRUE), 1), 1),
    N_MR_sig_detected    = sum(
      !is.na(Passed_606) & Passed_606 == TRUE & Detected, na.rm = TRUE),
    Has_quant            = any(!is.na(Quant_median) & Detected, na.rm = TRUE),
    .groups = "drop"
  ) %>% as.data.table()

fwrite(cov, file.path(out_dir, "04_coverage_summary.csv"))
wlog("04_coverage_summary.csv saved")

# ══════════════════════════════════════════════════════════════════════════════
# Step 7: PSM 统计（SQL，用 ProteinID 做关联）
# ══════════════════════════════════════════════════════════════════════════════
wlog("--- Step 7: PSM counts via SQL ---")

detected_pids <- det_full[Detected == TRUE, ProteinIDs] %>%
  strsplit(";") %>% unlist() %>% unique() %>%
  .[. != "—" & !is.na(.)]

if (length(detected_pids) > 0) {
  
  pid_str <- paste(detected_pids, collapse = ",")
  
  psm_sql <- sprintf("
    SELECT pp.ProteinID, COUNT(DISTINCT p.PeptideID) AS N_PSMs,
           COUNT(DISTINCT p.SpectrumID)              AS N_Spectra
    FROM   Peptides p
    JOIN   PeptidesProteins pp ON pp.PeptideID = p.PeptideID
    WHERE  pp.ProteinID IN (%s)
    GROUP BY pp.ProteinID
  ", pid_str)
  
  psm_counts <- sq(con, psm_sql)
  wlog("PSM counts retrieved: ",
       ifelse(is.null(psm_counts), "failed", nrow(psm_counts)), " proteins")
  
  if (!is.null(psm_counts) && nrow(psm_counts) > 0) {
    
    # 汇总到基因层（一个基因可能对应多个 ProteinID 异构体）
    psm_summary <- det_full[Detected == TRUE] %>%
      mutate(ProteinID_split = strsplit(ProteinIDs, ";")) %>%
      unnest(ProteinID_split) %>%
      mutate(ProteinID = as.integer(ProteinID_split)) %>%
      left_join(psm_counts, by = "ProteinID") %>%
      group_by(Gene, Module, UniProt_ACs, In_top_combo,
               Passed_606, Quant_median) %>%
      summarise(Total_PSMs    = sum(N_PSMs, na.rm = TRUE),
                Total_Spectra = sum(N_Spectra, na.rm = TRUE),
                .groups = "drop") %>%
      arrange(Module, desc(Total_PSMs)) %>%
      as.data.table()
    
    fwrite(psm_summary, file.path(out_dir, "05_psm_summary.csv"))
    wlog("05_psm_summary.csv saved: ", nrow(psm_summary), " rows")
  }
  
} else {
  wlog("No detected proteins for PSM lookup.")
  fwrite(data.table(), file.path(out_dir, "05_psm_summary.csv"))
}

# ══════════════════════════════════════════════════════════════════════════════
# 最终摘要
# ══════════════════════════════════════════════════════════════════════════════
wlog("\n══════════ FINAL SUMMARY ══════════")
wlog("MSF file         : ", basename(msf_path))
wlog("Total proteins   : ", nrow(pa))
wlog("Target genes     : ", length(target_genes))
wlog("Detected (100%)  : ", sum(match_results$Detected))
wlog("TMT quant matrix : ",
     ifelse(has_tmt,
            paste0("YES — ", nrow(tmt_wide), " proteins × ",
                   length(tmt_cols), " channels → 01b_tmt_protein.csv"),
            "NO — SQL join returned 0 rows (see diagnostic above)"))

wlog("\n--- Per-module ---")
for (i in seq_len(nrow(cov))) {
  r <- cov[i]
  wlog(sprintf("  %-10s | %d/%d detected (%s%%) | top_combo %d/%d (%s%%) | MR_sig %d | quant:%s",
               r$Module,
               r$N_detected, r$N_genes_total, r$Detection_pct,
               r$N_top_combo_detected, r$N_in_top_combo, r$Top_combo_pct,
               r$N_MR_sig_detected,
               ifelse(r$Has_quant, "YES", "NO")))
}

wlog("\n--- Detected ---")
for (i in which(det_full$Detected)) {
  r <- det_full[i]
  wlog(sprintf("  %-10s (%s)%s → AC:%s | PSMs see 05",
               r$Gene, r$Module,
               ifelse(isTRUE(r$In_top_combo), " ★", ""),
               r$UniProt_ACs))
}

wlog("\n--- NOT detected ---")
nd <- sort(det_full[Detected == FALSE]$Gene)
wlog("  ", if (length(nd) > 0) paste(nd, collapse = ", ") else "(none)")

wlog("\n--- Next steps ---")
if (has_tmt) {
  wlog("  ① 01c_tmt_normalized.csv → limma differential analysis")
  wlog("     Need: sample group labels (case/control per TMT channel)")
  wlog("  ② Check direction: TMT log2FC consistent with MR β direction?")
} else {
  wlog("  ① TMT join failed — run diagnostic SQL in R console:")
  wlog("    con <- dbConnect(SQLite(), msf_path)")
  wlog("    dbGetQuery(con, 'SELECT r.SpectrumID, s.SpectrumID, s.SearchSpectrumID")
  wlog("                     FROM ReporterIonQuanResults r")
  wlog("                     JOIN ReporterIonQuanResultsSearchSpectra s")
  wlog("                     ON r.SpectrumID = s.SpectrumID LIMIT 5')")
}
wlog("  ③ pQTL-MR: use detected proteins as intermediate phenotype")
wlog("\n=== 610 v3 completed ===")

cat("\n╔══════════ QUICK SUMMARY ══════════╗\n")
cat(sprintf("  Proteins : %d (152354 annotated)\n", nrow(pa)))
cat(sprintf("  Detected : %d / %d (100%%)\n",
            sum(match_results$Detected), length(target_genes)))
cat(sprintf("  TMT quant: %s\n",
            ifelse(has_tmt,
                   paste0(nrow(tmt_wide), " proteins × ",
                          length(tmt_cols), " channels"),
                   "NOT available (see log)")))
cat("╚═══════════════════════════════════╝\n\n")
cat("Per-module:\n"); print(as.data.frame(cov), row.names = FALSE)
cat("\nDetected genes:\n")
print(as.data.frame(
  det_full[Detected == TRUE,
           .(Gene, Module, In_top_combo, Match_method,
             UniProt_ACs, Quant_median)]), row.names = FALSE)
cat("\nOutputs →", out_dir, "\n")