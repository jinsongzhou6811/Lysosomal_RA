# ============================================================
# Figure3B_tile_diagram.R
# 自动化脚本：读取3个 functional_modules.csv → 分类基因
# → 生成出版级瓦片图 (600 dpi) + Cross-dataset module assignment of the eight seed lysosomal DEGs.csv
# ============================================================

library(rprojroot)
library(dplyr)
library(grid)

project_root <- find_rstudio_root_file()
outputs_dir  <- file.path(project_root, "outputs")

# ============================================================
# ██  用户可配置参数（所有调整只需修改此区域）  ██
# ============================================================

# -------- 数据路径 --------
fm_paths <- list(
  GSE55235 = file.path(outputs_dir, "04_Lysosomal_Pathway_Expansion", "GSE55235", "functional_modules.csv"),
  GSE55457 = file.path(outputs_dir, "04_Lysosomal_Pathway_Expansion", "GSE55457", "functional_modules.csv"),
  GSE55584 = file.path(outputs_dir, "04_Lysosomal_Pathway_Expansion", "GSE55584", "functional_modules.csv")
)
out_dir <- file.path(outputs_dir, "04_Lysosomal_Pathway_Expansion", "Core_Modules")

# -------- 画布尺寸 --------
FIG_W   <- 16      # 画布宽度（英寸）
FIG_H   <- 16      # 画布高度（英寸）——缩小此值可收紧上下排间距
FIG_DPI <- 600     # 输出分辨率

# -------- 字体大小（单位：pt） --------
FONT_TITLE       <- 19   # 模块标题，如 "Module 1 · GZMB–TUBA"
FONT_SUBTITLE    <- 16   # 副标题斜体，如 "Cytotoxic / microtubule axis"
FONT_STARS       <- 19   # ★★★ 保守星号
FONT_NOTE        <- 15   # 模块备注斜体
FONT_SECTION     <- 17   # 小节标题，如 "CORE GENES · 3/3"
FONT_GENE_TAG    <- 16   # 基因标签内文字
FONT_DS_LABEL    <- 14   # 标签旁数据集编号，如 "457·584"
FONT_BOTTOM_BAR  <- 15   # 面板底部 "GSE55235 · GSE55457 · GSE55584"
FONT_LEGEND_NUM  <- 16   # 图例栏数据集编号 "235"
FONT_LEGEND_NAME <- 15   # 图例栏 "= GSE55235"
FONT_LEGEND_TAG  <- 14   # 图例栏标签说明 "3/3 datasets ..."
FONT_LEGEND_STAR <- 15   # 图例栏 "★★★ = conserved 3/3"

# -------- 字体族与字形 --------
FONT_FAMILY       <- "sans"   # 字体族（"sans"=系统无衬线体）
FONT_FACE_GENE    <- "bold"   # 基因标签字形
FONT_FACE_DS      <- "bold"   # 数据集编号字形
FONT_FACE_TITLE   <- "bold"   # 标题字形
FONT_FACE_SECTION <- "bold"   # 小节标题字形

# -------- 四模块面板布局 --------
# c(x中心, y中心)：x控制左右(0=最左,1=最右)，y控制上下(0=最底,1=最顶)
# 上排y值 > 下排y值；缩小两排y差值 → 模块间距变小
PANEL_POS <- list(
  M1 = c(0.26, 0.78),   # 左上
  M2 = c(0.74, 0.78),   # 右上
  M3 = c(0.26, 0.45),   # 左下
  M4 = c(0.74, 0.45)    # 右下
)
PANEL_W <- 0.47          # 面板宽度（占画布比例）
# 每个模块面板高度（可独立调整，基因多的模块设大值）
PANEL_H <- list(
  M1 = 0.30,
  M2 = 0.30,
  M3 = 0.30,
  M4 = 0.30
)

# -------- 面板内部元素位置 --------
BOTTOM_BAR_Y <- 0.10    # 面板底部数据集文字的y位置（面板内坐标，0=底边，1=顶边）

# -------- 基因标签尺寸（面板视口内的比例） --------
TAG_W_CORE  <- 0.22     # 核心基因(3/3)标签宽度
TAG_W_EXT   <- 0.19     # 扩展/特异性基因标签宽度
TAG_H       <- 0.04     # 标签高度
COLS_CORE   <- 4        # 核心基因每行列数
COLS_EXT    <- 3        # 扩展/特异性基因每行列数
ROW_GAP     <- 0.08     # 同一小节内标签行间距
SECTION_GAP <- 0.10     # 小节之间的间距（CORE→EXTENDED→DS之间）

# -------- 图例栏整体 --------
LEGEND_H        <- 0.03     # 图例栏高度（占画布比例）
LEGEND_Y        <- 0.975    # 图例栏y中心位置
LEGEND_BG       <- "#2C3E50"   # 图例栏背景色
LEGEND_TEXT_COL <- "#BDC3C7"   # 图例栏浅色文字
LEGEND_STAR_COL <- "#F1C40F"   # 星号颜色

# -------- 图例栏内部布局（x坐标为图例视口内比例 0~1） --------
# 左侧：三个数据集的编号与全称位置
LEG_DS_X <- list(
  n235 = 0.02,    # "235" x位置
  t235 = 0.05,    # "= GSE55235" x位置
  n457 = 0.15,    # "457" x位置
  t457 = 0.18,    # "= GSE55457" x位置
  n584 = 0.28,    # "584" x位置
  t584 = 0.31     # "= GSE55584" x位置
)
# 右侧：保守性tag标签
LEG_TAG_X0  <- 0.43    # 第一个tag起始x位置（向右移可避免与左侧重叠）
LEG_TAG_GAP <- 0.14    # tag之间的间距（加大可避免tag之间重叠）
# 星号标记
LEG_STAR_X  <- 0.97    # "★★★ = conserved 3/3" x位置
LEG_STAR_Y  <- 0.5     # y位置（0.5=图例栏内垂直居中）

# -------- 模块颜色（与 Figure 3A 一致） --------
# M1=绿色, M2=蓝色, M3=橙棕色, M4=红褐色
MOD_COLORS <- list(
  M1 = list(bg="#EAF5E6", border="#5A8A50", header="#3A6A32", title_bg="#C8DFC0"),
  M2 = list(bg="#E8F0F8", border="#4A7FA5", header="#2C5F7C", title_bg="#BDD4E8"),
  M3 = list(bg="#FDF2E4", border="#C17F3E", header="#8B5E2B", title_bg="#EBCFA8"),
  M4 = list(bg="#F8ECE4", border="#B0624D", header="#7D3E2E", title_bg="#DDB8A8")
)

# -------- 基因标签颜色（按保守程度分三级） --------
TAG_COLORS <- list(
  `3` = list(fc="#5B7553", ec="#3D5A3A", tc="white"),     # 3/3 深色
  `2` = list(fc="#8FA888", ec="#6B8F6B", tc="#1C2E1A"),   # 2/3 中色
  `1` = list(fc="#C8D5C0", ec="#A0B098", tc="#3A4A34")    # 1/3 浅色
)

# -------- 模块元数据（标题、副标题、星号、备注） --------
MODULE_INFO <- list(
  M1 = list(title  = "Module 1 \u00b7 GZMB\u2013TUBA",
            sub    = "Cytotoxic / microtubule axis",
            stars  = "\u2605\u2605\u2605",
            note   = NULL),
  M2 = list(title  = "Module 2 \u00b7 LAMP3",
            sub    = "Lysosomal transport hub",
            stars  = "\u2605\u2605\u2605",
            note   = "Always separate from GZMB\u2013TUBA axis"),
  M3 = list(title  = "Module 3 \u00b7 NKG7\u2013PRF1",
            sub    = "Granule-mediated cytotoxicity",
            stars  = "\u2605\u2605\u2605",
            note   = NULL),
  M4 = list(title  = "Module 4 \u00b7 SLC39A8",
            sub    = "Zinc transport \u00b7 isolated hub",
            stars  = "\u2605\u2605\u2605",
            note   = "Always isolated from other seed genes")
)

# -------- 锚定基因（自动将Louvain模块映射到功能模块） --------
ANCHOR_GENES <- list(
  M1 = c("GZMB", "TUBA1A", "TUBA1B"),
  M2 = c("LAMP3"),
  M3 = c("NKG7", "PRF1"),
  M4 = c("SLC39A8")
)
# M3备用标记基因（当Louvain模块不含锚定基因时用于识别M3子模块）
M3_MARKERS <- c("TRAF3IP3", "LAG3", "PPP2CA", "CHMP7", "SIDT1",
                "FASLG", "CST7", "CTSW", "CXCR4")

# ============================================================
# ██  参数配置结束 — 以下代码无需修改  ██
# ============================================================

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
datasets <- names(fm_paths)

# ============================================================
# STEP 1: 读取三个数据集的 functional_modules.csv
# ============================================================
fm_list <- list()
for (ds in datasets) {
  if (!file.exists(fm_paths[[ds]])) stop(paste("文件未找到:", fm_paths[[ds]]))
  fm_list[[ds]] <- read.csv(fm_paths[[ds]], stringsAsFactors = FALSE)
  cat("已加载:", ds, "\u2014", nrow(fm_list[[ds]]), "个模块\n")
}

# ============================================================
# STEP 2: 将 Louvain 模块编号 → 功能模块（M1/M2/M3/M4）
# 逻辑：每个Louvain模块中若含有锚定基因则归入对应功能模块；
#        若不含锚定基因但含有M3标记基因则归入M3。
# ============================================================
map_modules <- function(fm_df) {
  mapping <- list()
  for (r in seq_len(nrow(fm_df))) {
    genes <- trimws(unlist(strsplit(fm_df$Genes[r], ",")))
    assigned <- FALSE
    for (fm in names(ANCHOR_GENES)) {
      if (any(ANCHOR_GENES[[fm]] %in% genes)) {
        mapping[[fm]] <- c(mapping[[fm]], list(genes))
        assigned <- TRUE
        break
      }
    }
    if (!assigned && any(M3_MARKERS %in% genes)) {
      mapping[["M3"]] <- c(mapping[["M3"]], list(genes))
    }
  }
  lapply(mapping, function(x) unique(unlist(x)))
}

all_maps <- lapply(fm_list, map_modules)

# ============================================================
# STEP 3: 跨数据集交集分类
# 3/3 = CORE（三个数据集都有）
# 2/3 = EXTENDED（两个数据集有）
# 1/3 = DATASET-SPECIFIC（仅一个数据集有）
# ============================================================
classify_genes <- function(fm_name) {
  gene_ds <- list()
  for (ds in datasets) {
    for (g in all_maps[[ds]][[fm_name]]) {
      gene_ds[[g]] <- unique(c(gene_ds[[g]], ds))
    }
  }
  core <- ext <- ds_spec <- data.frame(gene = character(), ds_label = character(),
                                       stringsAsFactors = FALSE)
  for (g in sort(names(gene_ds))) {
    n <- length(gene_ds[[g]])
    ds_lab <- paste(gsub("GSE55", "", sort(gene_ds[[g]])), collapse = "\u00b7")
    row <- data.frame(gene = g, ds_label = ds_lab, stringsAsFactors = FALSE)
    if (n == 3)      core    <- rbind(core, row)
    else if (n == 2) ext     <- rbind(ext, row)
    else             ds_spec <- rbind(ds_spec, row)
  }
  list(core = core, ext = ext, ds_spec = ds_spec)
}

module_data <- lapply(c("M1", "M2", "M3", "M4"), function(fm) {
  md <- classify_genes(fm)
  cat(fm, ": 核心=", nrow(md$core),
      " 扩展=", nrow(md$ext),
      " 特异=", nrow(md$ds_spec), "\n")
  md
})
names(module_data) <- c("M1", "M2", "M3", "M4")

# ============================================================
# STEP 4: 绘图辅助函数
# ============================================================

# --- 绘制单个基因标签 ---
draw_gene_tag <- function(vp_x, vp_y, gene_text, ds_text, level, tag_w) {
  tc <- TAG_COLORS[[as.character(level)]]
  # 圆角矩形背景
  grid.roundrect(x = unit(vp_x, "npc"), y = unit(vp_y, "npc"),
                 width = unit(tag_w, "npc"), height = unit(TAG_H, "npc"),
                 r = unit(0.008, "npc"),
                 gp = gpar(fill = tc$fc, col = tc$ec, lwd = 1.2))
  # 基因名文字
  grid.text(gene_text, x = unit(vp_x, "npc"), y = unit(vp_y, "npc"),
            gp = gpar(fontsize = FONT_GENE_TAG, fontface = FONT_FACE_GENE,
                      col = tc$tc, fontfamily = FONT_FAMILY))
  # 数据集标注（仅2/3和1/3显示）
  if (level < 3 && nchar(ds_text) > 0) {
    grid.text(ds_text, x = unit(vp_x + tag_w/2 + 0.015, "npc"),
              y = unit(vp_y, "npc"), just = "left",
              gp = gpar(fontsize = FONT_DS_LABEL, fontface = FONT_FACE_DS,
                        col = "#7F8C8D", fontfamily = FONT_FAMILY))
  }
}

# --- 按行排列一组基因标签，返回当前y位置 ---
draw_tag_rows <- function(df, y_start, level, tag_w, cols) {
  y <- y_start
  if (nrow(df) == 0) return(y)
  for (i in seq_len(nrow(df))) {
    col_idx <- (i - 1) %% cols
    if (col_idx == 0 && i > 1) y <- y - ROW_GAP   # 换行
    col_spacing <- if (level == 3) (tag_w + 0.005) else (tag_w + 0.12)
    x <- 0.06 + tag_w/2 + col_idx * col_spacing
    draw_gene_tag(x, y, df$gene[i], df$ds_label[i], level, tag_w)
  }
  return(y - SECTION_GAP)
}

# --- 绘制单个模块面板 ---
draw_module_panel <- function(fm_key, gd, info) {
  mc <- MOD_COLORS[[fm_key]]
  
  # 面板背景
  grid.roundrect(x = 0.5, y = 0.5, width = 0.96, height = 0.96,
                 r = unit(0.015, "npc"),
                 gp = gpar(fill = mc$bg, col = mc$border, lwd = 3))
  # 标题栏背景
  grid.roundrect(x = 0.5, y = 0.935, width = 0.96, height = 0.10,
                 r = unit(0.01, "npc"),
                 gp = gpar(fill = mc$title_bg, col = NA))
  # 模块标题
  grid.text(info$title, x = 0.06, y = 0.955, just = "left",
            gp = gpar(fontsize = FONT_TITLE, fontface = FONT_FACE_TITLE,
                      col = mc$header, fontfamily = FONT_FAMILY))
  # 保守星号
  grid.text(info$stars, x = 0.94, y = 0.955, just = "right",
            gp = gpar(fontsize = FONT_STARS, col = "#D4AC0D",
                      fontfamily = FONT_FAMILY))
  # 副标题
  grid.text(info$sub, x = 0.06, y = 0.905, just = "left",
            gp = gpar(fontsize = FONT_SUBTITLE, fontface = "italic",
                      col = "#6C7A89", fontfamily = FONT_FAMILY))
  
  y <- 0.84
  
  # 备注（如有）
  if (!is.null(info$note)) {
    grid.text(info$note, x = 0.06, y = y, just = "left",
              gp = gpar(fontsize = FONT_NOTE, fontface = "bold.italic",
                        col = mc$header, fontfamily = FONT_FAMILY))
    y <- y - 0.04
  }
  
  # ── CORE GENES · 3/3 ──
  grid.text("CORE GENES \u00b7 3/3", x = 0.06, y = y, just = "left",
            gp = gpar(fontsize = FONT_SECTION, fontface = FONT_FACE_SECTION,
                      col = mc$header, fontfamily = FONT_FAMILY))
  y <- y - 0.015
  grid.lines(x = c(0.06, 0.94), y = c(y, y),
             gp = gpar(col = mc$border, lwd = 0.8, alpha = 0.5))
  y <- y - 0.035
  y <- draw_tag_rows(gd$core, y, 3, TAG_W_CORE, COLS_CORE)
  
  # ── EXTENDED · 2/3 ──
  if (nrow(gd$ext) > 0) {
    grid.text("EXTENDED \u00b7 2/3", x = 0.06, y = y, just = "left",
              gp = gpar(fontsize = FONT_SECTION, fontface = FONT_FACE_SECTION,
                        col = "#5D6D7E", fontfamily = FONT_FAMILY))
    y <- y - 0.015
    grid.lines(x = c(0.06, 0.94), y = c(y, y),
               gp = gpar(col = "#BDC3C7", lwd = 0.8))
    y <- y - 0.035
    y <- draw_tag_rows(gd$ext, y, 2, TAG_W_EXT, COLS_EXT)
  }
  
  # ── DATASET-SPECIFIC · 1/3 ──
  if (nrow(gd$ds_spec) > 0) {
    grid.text("DATASET-SPECIFIC \u00b7 1/3", x = 0.06, y = y, just = "left",
              gp = gpar(fontsize = FONT_SECTION, fontface = FONT_FACE_SECTION,
                        col = "#95A5A6", fontfamily = FONT_FAMILY))
    y <- y - 0.015
    grid.lines(x = c(0.06, 0.94), y = c(y, y),
               gp = gpar(col = "#D5DBDB", lwd = 0.8))
    y <- y - 0.035
    draw_tag_rows(gd$ds_spec, y, 1, TAG_W_EXT, COLS_EXT)
  }
  
  # 面板底部数据集标签
  grid.roundrect(x = 0.5, y = BOTTOM_BAR_Y, width = 0.55, height = 0.05,
                 r = unit(0.008, "npc"),
                 gp = gpar(fill = "white", col = "#D5DBDB", lwd = 1))
  grid.text("GSE55235 \u00b7 GSE55457 \u00b7 GSE55584",
            x = 0.5, y = BOTTOM_BAR_Y,
            gp = gpar(fontsize = FONT_BOTTOM_BAR, fontface = "bold",
                      col = "#7F8C8D", fontfamily = FONT_FAMILY))
}

# ============================================================
# STEP 5: 组装完整图形（图例栏 + 四面板）
# ============================================================
draw_full_figure <- function() {
  grid.newpage()
  
  # ── 顶部图例栏 ──
  pushViewport(viewport(x = 0.5, y = LEGEND_Y, width = 0.96, height = LEGEND_H))
  grid.roundrect(gp = gpar(fill = LEGEND_BG, col = NA))
  
  # 数据集编号与全称（位置由 LEG_DS_X 控制）
  grid.text(expression(bold("235")), x = LEG_DS_X$n235, y = 0.5, just = "left",
            gp = gpar(fontsize = FONT_LEGEND_NUM, col = "white", fontfamily = FONT_FAMILY))
  grid.text("= GSE55235", x = LEG_DS_X$t235, y = 0.5, just = "left",
            gp = gpar(fontsize = FONT_LEGEND_NAME, fontface = "bold",
                      col = LEGEND_TEXT_COL, fontfamily = FONT_FAMILY))
  grid.text(expression(bold("457")), x = LEG_DS_X$n457, y = 0.5, just = "left",
            gp = gpar(fontsize = FONT_LEGEND_NUM, col = "white", fontfamily = FONT_FAMILY))
  grid.text("= GSE55457", x = LEG_DS_X$t457, y = 0.5, just = "left",
            gp = gpar(fontsize = FONT_LEGEND_NAME, fontface = "bold",
                      col = LEGEND_TEXT_COL, fontfamily = FONT_FAMILY))
  grid.text(expression(bold("584")), x = LEG_DS_X$n584, y = 0.5, just = "left",
            gp = gpar(fontsize = FONT_LEGEND_NUM, col = "white", fontfamily = FONT_FAMILY))
  grid.text("= GSE55584", x = LEG_DS_X$t584, y = 0.5, just = "left",
            gp = gpar(fontsize = FONT_LEGEND_NAME, fontface = "bold",
                      col = LEGEND_TEXT_COL, fontfamily = FONT_FAMILY))
  
  # 保守性tag标签（位置由 LEG_TAG_X0 和 LEG_TAG_GAP 控制）
  tag_labs <- c("3/3 datasets ", "2/3 datasets", "1/3 dataset")
  tag_fcs  <- c("#5B7553", "#8FA888", "#C8D5C0")
  tag_ecs  <- c("#3D5A3A", "#6B8F6B", "#A0B098")
  for (j in 1:3) {
    xo <- LEG_TAG_X0 + (j - 1) * LEG_TAG_GAP
    grid.roundrect(x = xo, y = 0.5, width = 0.018, height = 0.5,
                   r = unit(0.1, "npc"),
                   gp = gpar(fill = tag_fcs[j], col = tag_ecs[j], lwd = 1))
    grid.text(tag_labs[j], x = xo + 0.015, y = 0.5, just = "left",
              gp = gpar(fontsize = FONT_LEGEND_TAG, fontface = "bold",
                        col = "#ECF0F1", fontfamily = FONT_FAMILY))
  }
  
  # 星号标记（位置由 LEG_STAR_X 和 LEG_STAR_Y 控制）
  grid.text("\u2605\u2605\u2605 = conserved 3/3",
            x = LEG_STAR_X, y = LEG_STAR_Y, just = "right",
            gp = gpar(fontsize = FONT_LEGEND_STAR, fontface = "bold",
                      col = LEGEND_STAR_COL, fontfamily = FONT_FAMILY))
  popViewport()
  
  # ── 四个模块面板 ──
  for (fm in c("M1", "M2", "M3", "M4")) {
    pushViewport(viewport(x = PANEL_POS[[fm]][1], y = PANEL_POS[[fm]][2],
                          width = PANEL_W, height = PANEL_H[[fm]]))
    draw_module_panel(fm, module_data[[fm]], MODULE_INFO[[fm]])
    popViewport()
  }
}

# ============================================================
# STEP 6: 保存图像
# ============================================================
png(file.path(out_dir, "Figure3B_600dpi.png"),
    width = FIG_W, height = FIG_H, units = "in", res = FIG_DPI)
draw_full_figure()
dev.off()

tiff(file.path(out_dir, "Figure3B_600dpi.tiff"),
     width = FIG_W, height = FIG_H, units = "in", res = FIG_DPI,
     compression = "lzw")
draw_full_figure()
dev.off()

cat("\n\u2713 Figure 3B 已保存至:", out_dir, "\n")
cat("  PNG:  Figure3B_600dpi.png\n")
cat("  TIFF: Figure3B_600dpi.tiff\n")

# ============================================================
# STEP 7: 导出 Cross-dataset module assignment of the eight seed lysosomal DEGs.csv（含全部 3/3 + 2/3 + 1/3 基因）
#         每个基因显示在各数据集中的 Louvain 模块编号
# ============================================================
module_names <- c(M1 = "GZMB-TUBA", M2 = "LAMP3", M3 = "NKG7-PRF1", M4 = "SLC39A8")
module_descs <- c(M1 = "Cytotoxic / microtubule axis",
                  M2 = "Lysosomal transport hub",
                  M3 = "Granule-mediated cytotoxicity",
                  M4 = "Zinc transport / isolated hub")

# --- 辅助函数：查找基因在某个数据集中的 Louvain 模块编号及功能标签 ---
get_louvain_module <- function(gene, ds) {
  fm_df <- fm_list[[ds]]
  for (r in seq_len(nrow(fm_df))) {
    genes_in_mod <- trimws(unlist(strsplit(fm_df$Genes[r], ",")))
    if (gene %in% genes_in_mod) {
      # 获取Louvain模块编号
      mod_num <- fm_df$Module[r]
      # 判断该Louvain模块属于哪个功能模块
      func_label <- ""
      for (fm_key in names(ANCHOR_GENES)) {
        if (any(ANCHOR_GENES[[fm_key]] %in% genes_in_mod)) {
          func_label <- paste0(" (", module_names[fm_key], ")")
          break
        }
      }
      return(paste0(mod_num, func_label))
    }
  }
  return("Not detected")
}

# --- 构建完整的基因表（遍历四个功能模块的全部基因） ---
all_rows <- list()
for (fm in c("M1", "M2", "M3", "M4")) {
  md <- module_data[[fm]]
  
  # 合并三个层级的基因
  levels_list <- list(
    list(df = md$core,    conserved = "Core 3/3")
  )
  if (nrow(md$ext) > 0) {
    levels_list[[length(levels_list) + 1]] <- list(df = md$ext, conserved = "Extended 2/3")
  }
  if (nrow(md$ds_spec) > 0) {
    levels_list[[length(levels_list) + 1]] <- list(df = md$ds_spec, conserved = "Dataset-specific 1/3")
  }
  
  for (lvl in levels_list) {
    if (nrow(lvl$df) == 0) next
    for (i in seq_len(nrow(lvl$df))) {
      gene <- lvl$df$gene[i]
      all_rows[[length(all_rows) + 1]] <- data.frame(
        Module               = fm,
        Module_Name          = module_names[fm],
        Description          = module_descs[fm],
        Gene                 = gene,
        GSE55235_Module      = get_louvain_module(gene, "GSE55235"),
        GSE55457_Module      = get_louvain_module(gene, "GSE55457"),
        GSE55584_Module      = get_louvain_module(gene, "GSE55584"),
        Final_Classification = lvl$conserved,
        Assigned_Module      = paste0(fm, ": ", module_names[fm]),
        stringsAsFactors     = FALSE
      )
    }
  }
}

core_df <- do.call(rbind, all_rows)
core_csv_path <- file.path(out_dir, "Cross-dataset module assignment of the eight seed lysosomal DEGs.csv")
write.csv(core_df, core_csv_path, row.names = FALSE)

cat("\n\u2713 Cross-dataset module assignment of the eight seed lysosomal DEGs.csv 已保存:", core_csv_path, "\n")
cat("  基因总数:", nrow(core_df),
    " (Core 3/3:", sum(core_df$Final_Classification == "Core 3/3"),
    " | Extended 2/3:", sum(core_df$Final_Classification == "Extended 2/3"),
    " | Dataset-specific 1/3:", sum(core_df$Final_Classification == "Dataset-specific 1/3"), ")\n")

# ============================================================
# STEP 8: 打印分类结果摘要
# ============================================================
cat("\n=== 基因分类结果摘要 ===\n")
for (fm in c("M1", "M2", "M3", "M4")) {
  md <- module_data[[fm]]
  cat(sprintf("\n%s: %d 核心 | %d 扩展 | %d 特异\n",
              fm, nrow(md$core), nrow(md$ext), nrow(md$ds_spec)))
  if (nrow(md$core) > 0)
    cat("  Core:", paste(md$core$gene, collapse = ", "), "\n")
  if (nrow(md$ext) > 0)
    cat("  Ext:",  paste0(md$ext$gene, "(", md$ext$ds_label, ")", collapse = ", "), "\n")
  if (nrow(md$ds_spec) > 0)
    cat("  DS:",   paste0(md$ds_spec$gene, "(", md$ds_spec$ds_label, ")", collapse = ", "), "\n")
}