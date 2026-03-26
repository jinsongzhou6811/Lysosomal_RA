# ========================================================
# Standard Path Management Template）
# ========================================================
library(rprojroot)
project_root <- find_rstudio_root_file()

outputs_dir <- file.path(project_root, "outputs")

# Enter the directory 
integrated_dir <- file.path(outputs_dir, "03_WGCNA", "Integrated")

# Output directory
output_dir <- file.path(outputs_dir, "03_WGCNA", "PPI_Cytoscape_edges")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load package
library(igraph)
library(dplyr)
library(ggnetwork)
library(ggrepel)

# ====================== Core_genes ======================
genes <- c("GZMB", "LAMP3", "MREG", "NKG7", "PRF1", "SLC39A8", "TRAF3IP3", "VOPP1")
core_genes <- genes

# ====================== Read 8 edge files ======================
edges_list <- list()
for (gene in genes) {
  file_path <- file.path(integrated_dir, gene, paste0(gene, "_Integrated_Analysis"), 
                         paste0(gene, "_Cytoscape_edges.txt"))
  
  if (!file.exists(file_path)) stop(paste("The file doesnot exist：", file_path))
  
  edges <- read.csv(file_path, sep = "\t", fileEncoding = "UTF-8") %>%
    select(fromNode, toNode, weight, module)
  
  edges_list[[gene]] <- edges
}

# ====================== Merge & Filter ======================
all_edges <- bind_rows(edges_list) %>%
  distinct() %>%
  filter(weight > 0.05 & 
           module %in% c("magenta", "grey60", "pink", "ivory", "brown", "darkturquoise", "violet"))

# ====================== Build a network ======================
g <- graph_from_data_frame(all_edges, directed = FALSE,
                           vertices = data.frame(name = unique(c(all_edges$fromNode, 
                                                                 all_edges$toNode, 
                                                                 core_genes))))

V(g)$is_core <- V(g)$name %in% core_genes
V(g)$color   <- ifelse(V(g)$is_core, "red", "lightblue")
V(g)$size    <- ifelse(V(g)$is_core, 10, 5)

# ====================== Export GraphML ======================
write_graph(g, file.path(output_dir, "theme_ppi.graphml"), format = "graphml")
cat("✓ Generated：theme_ppi.graphml\n")

# ====================== Network visualization ======================
layout <- layout_with_fr(g, weights = E(g)$weight, niter = 1000, grid = "nogrid")
net <- ggnetwork(g, layout = layout)

p <- ggplot(net, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(size = weight), color = "black", curvature = 0.1,
             arrow = arrow(length = unit(0.2, "cm"), type = "closed", ends = "last")) +
  geom_nodes(aes(color = color, size = size), shape = 21, stroke = 0.5) +
  ggrepel::geom_text_repel(aes(label = name),
                           data = filter(net, is_core),
                           size = 3.8, max.overlaps = 30, min.segment.length = 0) +
  theme_void() +
  labs(title = "PPI Network of Lysosomal Genes in RA") +
  scale_color_identity() +
  scale_size_continuous(range = c(3, 15))

# 保存 PDF 和 PNG
ggsave(file.path(output_dir, "ppi_theme.pdf"), p, width = 11, height = 9, device = "pdf")
ggsave(file.path(output_dir, "ppi_theme.png"), p, width = 11, height = 9, device = "png")

cat("✓ Generated：ppi_theme.pdf\n")
cat("✓ Generated：ppi_theme.png\n")
cat("\nAll files have been saved to：\n", output_dir, "\n")