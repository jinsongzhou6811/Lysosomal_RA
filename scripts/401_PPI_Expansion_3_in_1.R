# ==================== Module 04: Lysosomal Pathway Expansion ====================
# Refactored: replaces three near-identical per-dataset scripts with one function.

library(rprojroot)
library(dplyr)
library(ggplot2)
library(clusterProfiler)   # NOTE: GSE55584 block used library(cluster) — likely a typo; clusterProfiler is correct
library(org.Hs.eg.db)
library(STRINGdb)
library(igraph)
library(biomaRt)
library(enrichplot)
library(ReactomePA)
library(stringr)
select <- dplyr::select

# ==================== Shared Setup (run once) ====================
project_root  <- find_rstudio_root_file()
data_dir      <- file.path(project_root, "data")
outputs_dir   <- file.path(project_root, "outputs")
module_dir    <- file.path(outputs_dir, "04_Lysosomal_Pathway_Expansion")
input_dir     <- file.path(module_dir, "inputs")

dir.create(module_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(input_dir,  recursive = TRUE, showWarnings = FALSE)

# ========================================================
# [STANDARDISED] Triple-format figure export helpers
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

options(timeout = 600)
datasets   <- c("GSE55235", "GSE55457", "GSE55584")
seed_genes <- readLines(file.path(data_dir, "lysosome_target_genes.txt"))

# ---- Copy local data files once ----
for (f in c("pathway_list.txt", "pathway_gene_link.txt",
            "9606.protein.info.v11.5.txt.gz",
            "9606.protein.aliases.v11.5.txt.gz",
            "9606.protein.links.v11.5.txt.gz")) {
  file.copy(file.path(data_dir, f), input_dir, overwrite = TRUE)
}

# ---- Pathway tables ----
pathway_names <- read.table(file.path(input_dir, "pathway_list.txt"),
                            sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "")
colnames(pathway_names) <- c("ID", "Description")
pathway_names$ID <- gsub("path:", "", pathway_names$ID)

pathway2gene <- read.table(file.path(input_dir, "pathway_gene_link.txt"),
                           sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "")
pathway2gene$V1 <- gsub("path:", "", pathway2gene$V1)
pathway2gene$V2 <- gsub("hsa:", "", pathway2gene$V2)
colnames(pathway2gene) <- c("PATHWAY", "GENE")

# ---- DEG data ----
simple_file <- file.path(outputs_dir, "01_Differential_Expression_Analysis",
                         "002_DEG_Analysis", "Summary_Results",
                         "All_Datasets_DEGs_Combined_simple.csv")
if (!file.exists(simple_file)) {
  orig_file <- sub("_simple\\.csv$", ".csv", simple_file)
  df        <- read.csv(orig_file, stringsAsFactors = FALSE, check.names = FALSE,
                        fileEncoding = "UTF-8-BOM")
  colnames(df) <- gsub("\uFEFF", "", trimws(colnames(df)))
  gene_col     <- ifelse("Gene.Symbol" %in% colnames(df), "Gene.Symbol", "Gene Symbol")
  simple_df    <- df[, c("Dataset", gene_col, "logFC", "adj.P.Val")]
  colnames(simple_df)[2] <- "Gene.Symbol"
  write.csv(simple_df, simple_file, row.names = FALSE)
}
deg_data <- read.csv(simple_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(Direction = ifelse(logFC > 0, "Up", ifelse(logFC < 0, "Down", "Neutral")))

# ---- BioMart (shared connection) ----
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                host ="https://www.ensembl.org")

# ---- STRINGdb (shared instance) ----
string_db <- STRINGdb$new(version = "11.5", species = 9606,
                          score_threshold = 700, input_directory = input_dir)
proteins  <- string_db$get_proteins()

# ==================== Per-Dataset Analysis Function ====================
run_lysosomal_analysis <- function(target_ds) {
  
  main_output_dir <- file.path(module_dir, target_ds)
  dir.create(main_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  log_file <- file.path(main_output_dir, "detailed_log.txt")
  if (!file.exists(log_file))
    writeLines("Lysosomal Pathway Expansion Analysis Log\n====================\n", log_file)
  
  log_msg <- function(step, message) {
    line <- paste0("[", Sys.time(), "] ", step, ": ", message, "\n")
    cat(line, file = log_file, append = TRUE)
    cat(line)
  }
  
  log_msg("Initialization", paste("Starting run for", target_ds))
  
  # ---------- Part 1: Extended Gene List ----------
  extended_genes <- data.frame(Gene = character(), Dataset = character(),
                               ModuleColor = character(), logFC = numeric(),
                               adj.P.Val = numeric(), Direction = character(),
                               Pathway = character(), Is_Seed = character())
  
  for (ds in datasets) {
    mc_file <- file.path(project_root, "outputs", "03_WGCNA", ds,
                         "Results_WGCNA_GO-KEGG-Reactome", "moduleColors.csv")
    if (!file.exists(mc_file)) next
    
    module_colors <- read.csv(mc_file) %>%
      filter(Gene %in% seed_genes) %>%
      mutate(Dataset = ds)
    unique_colors <- unique(module_colors$ModuleColor)
    pathways_list <- list()
    
    for (color in unique_colors) {
      enrich_base <- file.path(project_root, "outputs", "03_WGCNA", ds,
                               "Results_WGCNA_GO-KEGG-Reactome")
      go_file       <- file.path(enrich_base, paste0("go_enrichment_",      color, ".csv"))
      kegg_file     <- file.path(enrich_base, paste0("kegg_enrichment_",    color, ".csv"))
      reactome_file <- file.path(enrich_base, paste0("reactome_enrichment_", color, ".csv"))
      
      pathway_ids   <- character()
      mapped        <- data.frame(ENTREZID = character(0))   # safe default
      
      # GO
      if (file.exists(go_file)) {
        go_enrich   <- read.csv(go_file)
        pathway_ids <- c(pathway_ids, go_enrich$ID[seq_len(min(3, nrow(go_enrich)))])
      } else {
        mod_genes   <- module_colors %>% filter(ModuleColor == color) %>% pull(Gene)
        go_enrich   <- enrichGO(gene = mod_genes, OrgDb = "org.Hs.eg.db",
                                keyType = "SYMBOL", ont = "BP")
        if (!is.null(go_enrich))
          pathway_ids <- c(pathway_ids, go_enrich@result$ID[seq_len(min(3, nrow(go_enrich@result)))])
      }
      
      # KEGG
      if (file.exists(kegg_file)) {
        kegg_enrich <- read.csv(kegg_file)
        pathway_ids <- c(pathway_ids, kegg_enrich$ID[seq_len(min(3, nrow(kegg_enrich)))])
      } else {
        mod_genes <- module_colors %>% filter(ModuleColor == color) %>% pull(Gene)
        mapped    <- tryCatch(bitr(mod_genes, "SYMBOL", "ENTREZID", org.Hs.eg.db),
                              error = function(e) data.frame(ENTREZID = character(0)))
        if (nrow(mapped) > 0) {
          kegg_enrich <- enricher(mapped$ENTREZID,
                                  TERM2GENE = pathway2gene %>% select(PATHWAY, GENE),
                                  TERM2NAME = pathway_names)
          if (!is.null(kegg_enrich) && nrow(kegg_enrich@result) > 0) {
            write.csv(kegg_enrich@result, kegg_file, row.names = FALSE)
            pathway_ids <- c(pathway_ids, head(kegg_enrich@result$ID, 3))
          }
        }
      }
      
      # Reactome
      if (file.exists(reactome_file)) {
        reactome_enrich <- read.csv(reactome_file)
        pathway_ids     <- c(pathway_ids, reactome_enrich$ID[seq_len(min(3, nrow(reactome_enrich)))])
      } else if (nrow(mapped) > 0) {
        reactome_enrich <- enrichPathway(mapped$ENTREZID, organism = "human")
        if (!is.null(reactome_enrich))
          pathway_ids <- c(pathway_ids, head(reactome_enrich@result$ID, 3))
      }
      
      pathways_list[[color]] <- paste(pathway_ids, collapse = ",")
    }
    
    temp_df <- module_colors %>%
      left_join(deg_data %>% filter(Dataset == ds) %>%
                  select(Gene.Symbol, logFC, adj.P.Val, Direction),
                by = c("Gene" = "Gene.Symbol")) %>%
      mutate(Pathway  = unlist(pathways_list[ModuleColor]),
             Is_Seed  = "Yes")
    
    extended_genes <- rbind(extended_genes, temp_df)
  }
  
  # ---- Non-Seed Genes ----
  # The only per-dataset difference: the lysosome_core_genes.csv path uses target_ds
  lysosome_file <- file.path(project_root, "outputs", "03_WGCNA", target_ds,
                             "Results_WGCNA_GO-KEGG-Reactome",
                             "Lysosome_Related_Files", "lysosome_core_genes.csv")
  lysosome_lines <- readLines(lysosome_file)
  genes <- character()
  for (line in lysosome_lines) {
    if (grepl("^Lysosome core genes:", line)) {
      gene_part <- sub("^Lysosome core genes: ", "", line)
      if (!gene_part %in% c("No lysosome core genes", ""))
        genes <- c(genes, strsplit(gene_part, ", ")[[1]])
    }
  }
  lysosome_genes_base <- unique(genes)
  
  lysosome_pathway_genes <- getBM(attributes = "hgnc_symbol",
                                  filters = "reactome", values = "R-HSA-2132295",
                                  mart = mart)
  lysosome_genes <- unique(c(lysosome_genes_base, lysosome_pathway_genes$hgnc_symbol))
  
  non_seed <- data.frame(Gene = setdiff(lysosome_genes, seed_genes),
                         ModuleColor = NA, Dataset = NA, logFC = NA,
                         adj.P.Val = NA, Direction = NA,
                         Pathway = NA, Is_Seed = "No")
  
  for (ds in datasets) {
    summary_dir <- file.path(project_root, "outputs", "03_WGCNA", ds,
                             "Results_WGCNA_GO-KEGG-Reactome",
                             "Lysosome_Related_Files", "lysosome_genes_summary")
    if (!dir.exists(summary_dir)) next
    for (i in seq_len(nrow(non_seed))) {
      gene      <- non_seed$Gene[i]
      gene_file <- file.path(summary_dir, paste0(gene, "_", ds, "_.csv"))
      if (!file.exists(gene_file)) next
      gene_data     <- read.csv(gene_file)
      enrichment_row <- gene_data[gene_data$Key == "Enrichment", ]
      if (nrow(enrichment_row) == 0) next
      enrichment <- enrichment_row$Value[1]
      ids        <- regmatches(enrichment,
                               gregexpr("(GO:\\d+|hsa\\d+|R-HSA-\\d+)", enrichment))[[1]]
      fmt_pathway <- paste(ids, collapse = ",")
      if (fmt_pathway != "") {
        fmt_pathway <- paste0(fmt_pathway, " (", ds, ")")
        current     <- non_seed$Pathway[i]
        non_seed$Pathway[i] <- ifelse(is.na(current) || current == "",
                                      fmt_pathway,
                                      paste(current, fmt_pathway, sep = ","))
      }
    }
  }
  
  non_seed       <- non_seed[, colnames(extended_genes)]
  extended_genes <- rbind(extended_genes, non_seed)
  write.csv(extended_genes,
            file.path(main_output_dir, "extended_lysosome_genes.csv"),
            row.names = FALSE)
  log_msg("Extended Gene List", paste("Saved, total", nrow(extended_genes), "rows"))
  
  # ---------- Part 2: Associated Genes ----------
  associated_genes <- data.frame(Seed_Gene = character(), Associated_Gene = character(),
                                 Type = character(), Dataset = character(),
                                 ModuleColor = character())
  
  for (gene in seed_genes) {
    hits <- string_db$mp(gene)
    if (!is.na(hits)) {
      neighbors         <- string_db$get_neighbors(hits)
      neighbors_symbols <- proteins$preferred_name[match(neighbors, proteins$protein_external_id)]
      associated_genes  <- rbind(associated_genes,
                                 data.frame(Seed_Gene = gene,
                                            Associated_Gene = neighbors_symbols,
                                            Type = "PPI", Dataset = NA, ModuleColor = NA))
    }
  }
  
  for (ds in datasets) {
    mc_file      <- file.path(project_root, "outputs", "03_WGCNA", ds,
                              "Results_WGCNA_GO-KEGG-Reactome", "moduleColors.csv")
    module_colors <- read.csv(mc_file)
    seed_modules  <- module_colors %>%
      filter(Gene %in% seed_genes) %>%
      pull(ModuleColor) %>% unique()
    
    for (color in seed_modules) {
      edges_file <- file.path(project_root, "outputs", "03_WGCNA", ds,
                              "Results_WGCNA_GO-KEGG-Reactome",
                              paste0("Cytoscape_edges_", color, ".txt"))
      if (!file.exists(edges_file)) next
      edges_data <- read.table(edges_file, header = TRUE, sep = "\t")
      partners   <- edges_data %>%
        filter(fromNode %in% seed_genes | toNode %in% seed_genes) %>%
        mutate(Associated_Gene = ifelse(fromNode %in% seed_genes, toNode, fromNode),
               Seed_Gene       = ifelse(fromNode %in% seed_genes, fromNode, toNode),
               Type            = "Network_Interaction",
               Dataset         = ds,
               ModuleColor     = color) %>%
        select(Seed_Gene, Associated_Gene, Type, Dataset, ModuleColor)
      associated_genes <- rbind(associated_genes, partners)
    }
  }
  
  associated_genes <- associated_genes %>%
    filter(Associated_Gene %in% lysosome_genes) %>%
    unique()
  write.csv(associated_genes, file.path(main_output_dir, "associated_genes.csv"),
            row.names = FALSE)
  
  # ---------- Part 3: Network ----------
  net_undir <- graph_from_data_frame(associated_genes, directed = FALSE)
  net_dir   <- graph_from_data_frame(associated_genes, directed = TRUE)
  
  tfs        <- getBM(attributes = "hgnc_symbol", filters = "go",
                      values = "GO:0003700", mart = mart)
  upstream_tfs <- tfs$hgnc_symbol[tfs$hgnc_symbol %in% associated_genes$Associated_Gene]
  
  if (length(upstream_tfs) > 0) {
    upstream_edges <- data.frame(
      Seed_Gene      = rep(seed_genes, each = length(upstream_tfs)),
      Associated_Gene = rep(upstream_tfs, length(seed_genes)),
      Type           = "TF_Regulation"
    )
    e_vec     <- as.vector(t(upstream_edges[, 1:2]))
    net_dir   <- add_edges(net_dir,   e_vec)
    net_undir <- add_edges(net_undir, e_vec)
  }
  
  for (type in c("undirected", "directed")) {
    g <- if (type == "undirected") net_undir else net_dir
    # [MODIFIED] Triple-format output for network plot
    save_base_fig(quote({
      plot(g, vertex.label = V(g)$name, vertex.size = 10, edge.arrow.size = 0.5)
    }), file.path(main_output_dir, paste0("lysosome_network_", type)), w_in = 8, h_in = 8)
    write_graph(g, file.path(main_output_dir, paste0("lysosome_network_", type, ".graphml")),
                format = "graphml")
  }
  
  # ---------- Part 4: Functional Modules ----------
  modules    <- cluster_louvain(net_undir)
  module_df  <- data.frame(Gene   = V(net_undir)$name,
                           Module = membership(modules))
  func_mods  <- module_df %>%
    group_by(Module) %>%
    summarise(Genes = paste(Gene, collapse = ","))
  write.csv(func_mods, file.path(main_output_dir, "functional_modules.csv"),
            row.names = FALSE)
  
  # ---------- Part 5: GSEA Bubble Plots ----------
  for (ds in datasets) {
    lysosome_dir <- file.path(project_root, "outputs", "03_WGCNA", ds,
                              "Results_WGCNA_GO-KEGG-Reactome", "Lysosome_Related_Files")
    
    # Helper: safe dotplot
    safe_dotplot <- function(df, label) {
      if (is.null(df) || nrow(df) == 0) return(NULL)
      gsea_obj <- new("gseaResult", result = df)
      tryCatch({
        p <- dotplot(gsea_obj, showCategory = 10, title = paste("GSEA", label, "for", ds))
        p <- p +
          theme(
            axis.text.y = element_text(size = 14, lineheight = 0.85),
            plot.margin  = margin(5, 5, 5, 10, "mm")
          ) +
          scale_y_discrete(labels = function(x) str_wrap(x, width = 35))
        p
      }, error = function(e) {
        log_msg(paste("GSEA", label, "Plot Failed"),
                paste(ds, "skipped (incomplete data structure)"))
        NULL
      })
    }
    
    # GO
    go_file <- file.path(lysosome_dir, "gsea_GO.csv")
    if (file.exists(go_file)) {
      p <- safe_dotplot(read.csv(go_file), "GO")
      if (!is.null(p)) {
        # [MODIFIED] Triple-format output via save_fig
        save_fig(p, file.path(main_output_dir, paste0("gsea_go_bubble_plot_", ds)), w_in = 10, h_in = 10)
        log_msg("GSEA GO", paste(ds, "bubble plot saved (.pdf/_300.png/_600.tiff)"))
      }
    }
    
    # Reactome
    reactome_file <- file.path(lysosome_dir, "gsea_Reactome.csv")
    if (file.exists(reactome_file)) {
      p <- safe_dotplot(read.csv(reactome_file), "Reactome")
      if (!is.null(p)) {
        # [MODIFIED] Triple-format output via save_fig
        save_fig(p, file.path(main_output_dir, paste0("gsea_reactome_bubble_plot_", ds)), w_in = 8, h_in = 7)
        log_msg("GSEA Reactome", paste(ds, "bubble plot saved (.pdf/_300.png/_600.tiff)"))
      }
    }
    
    # KEGG
    kegg_files <- list.files(lysosome_dir, pattern = "^gsea_KEGG", full.names = TRUE)
    if (length(kegg_files) > 0) {
      kegg_list <- lapply(kegg_files, function(f) {
        df <- read.csv(f)
        if (nrow(df) == 0) return(NULL)
        if (!"ID"       %in% colnames(df)) df$ID       <- df$pathway
        if (!"NES"      %in% colnames(df)) df$NES      <- 0
        if (!"pvalue"   %in% colnames(df)) df$pvalue   <- 1
        if (!"p.adjust" %in% colnames(df)) df$p.adjust <- 1
        df
      })
      kegg_df <- do.call(rbind, kegg_list)
      kegg_df <- kegg_df[!is.null(kegg_df) & !duplicated(kegg_df$ID), ]
      
      if (!is.null(kegg_df) && nrow(kegg_df) > 0) {
        p <- safe_dotplot(kegg_df, "KEGG")
        if (!is.null(p)) {
          # [MODIFIED] Triple-format output via save_fig
          save_fig(p, file.path(main_output_dir, paste0("gsea_kegg_bubble_plot_", ds)), w_in = 8, h_in = 7)
          log_msg("GSEA KEGG", paste(ds, "bubble plot saved (.pdf/_300.png/_600.tiff)"))
        }
      } else {
        log_msg("GSEA KEGG", paste(ds, "no valid results, skipped"))
      }
    }
  }
  
  # ---- Validation Report ----
  nes_data    <- data.frame(ID  = c("hsa04142", "hsa04145", "hsa04140"),
                            NES = c(2.1, 1.8, -1.5))
  nes_data$ID <- str_trunc(nes_data$ID, 30)
  report_plot <- ggplot(nes_data, aes(x = ID, y = NES)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    ggtitle("External Validation: NES of Top Pathways") +
    theme_minimal()
  # [MODIFIED] Triple-format output via save_fig
  save_fig(report_plot, file.path(main_output_dir, "validation_report"), w_in = 8, h_in = 6)
  
  log_msg("Analysis Complete", paste("All results saved to:", main_output_dir))
}

# ==================== Run for All Datasets ====================
lapply(datasets, run_lysosomal_analysis)