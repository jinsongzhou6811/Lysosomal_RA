# Lysosomal_RA

A **complete, reproducible multi-omics bioinformatics pipeline** for the study of lysosomal gene co-expression modules and their synergistic role with the cytoskeleton in rheumatoid arthritis (RA) pathogenesis.

**Paper**: Lysosomal-Cytoskeletal Co-expression Modules Drive Rheumatoid Arthritis Pathogenesis: Causal Evidence from Multi-Omics and Mendelian Randomization — Jinsong Zhou et al.

---

## What This Pipeline Does

- Bulk RNA-seq differential expression analysis (DEG) across 3 GEO datasets
- Lysosomal DEG intersection against a curated 719-gene UniProt reference
- Weighted gene co-expression network analysis (WGCNA) + Louvain module detection
- Protein–protein interaction (PPI) network construction and expansion
- Single-cell RNA-seq analysis (3 cohorts, Seurat + Monocle3)
- Mendelian randomization: single-gene (UVMR) and weighted-combination analyses
- TMT proteomics validation

---

## Repository Structure

```
Lysosomal_RA/
├── data/                          # Metadata, processed files, and gene lists
│   ├── DOWNLOAD_INSTRUCTIONS.md   # ⚠️ Read this first — ~54 GB raw data
│   ├── check_data_ready.R         # Verify all required files are present
│   ├── GWAS/                      # Per-gene eQTL / sQTL / GWAS association files
│   └── ...                        # Gene lists, module files, UVMR config files
├── scripts/                       # 35 analysis scripts across 7 modules (run in order)
└── outputs/                       # Auto-generated figures and tables (not committed)
```

---

## Data Preparation — Read First!

All large raw data files (~54 GB total) **must be downloaded manually** before running anything.

1. Open [`data/DOWNLOAD_INSTRUCTIONS.md`](data/DOWNLOAD_INSTRUCTIONS.md)
2. Download all 29 files:
   - 6 STRINGdb protein network files
   - 4 GEO microarray raw `.tar` files + 1 scRNA-seq `.rds`
   - 4 cisTarget `.feather` databases
   - 14 GWAS summary statistics `.vcf.gz` files
3. Extract the 4 `GSE*.tar` files into their corresponding subfolders
4. Verify everything is in place:
   ```r
   source("data/check_data_ready.R")
   ```

---

## Running the Pipeline

Open `Lysosomal_RA.Rproj` in RStudio, then run scripts **in numerical order**:

| Module | Scripts | Description |
|--------|---------|-------------|
| 01 | `101`, `102` | GEO preprocessing (ComBat batch correction) + DEG analysis |
| 02 | `201` | Lysosomal DEG enrichment intersection |
| 03 | `301`–`308` | WGCNA + enrichment + PPI + module unification |
| 04 | `401` | PPI network expansion (`401_PPI_Expansion 3 in 1.R` recommended) |
| 05 | `501`–`505` | scRNA-seq analysis — see Seurat version note below |
| 06 | `601`–`609` | SNP extraction + UVMR (single-gene and weighted-combination) |
| 07 | `701`–`703` | TMT proteomics validation + publication figures |

### ⚠️ Seurat Version Note

Scripts `501`, `502`, `503` have two variants each:
- **Standard** (`501_scRNA_Analysis_GSE159117.R`) — tested with **Seurat v4.4.0**
- **Font-swap variant** (`501_scRNA_Analysis_GSE159117 - 字体更换.R`) — alternate font settings for publication figures

Use the standard scripts unless you need to reproduce the exact publication fonts. Seurat v4/v5 API differences (especially `Assay5` / `Layers`) are the most common source of errors — **do not mix versions** within a session.

### Other Script Notes

- **`401`**: Two versions exist (`3 datasets` vs `3 in 1`). The `3 in 1` version is the consolidated pipeline script and is recommended.
- **`602`**: The file on disk is named `602_Combination_Weight_Tables_Generator.R.R` (double extension). Rename to `.R` before running if your system does not handle this automatically.

---

## System Requirements

- **R ≥ 4.4** (R 4.4.x recommended)
- **RStudio** strongly recommended (use the `.Rproj` file — never open scripts directly)
- **RAM**: 32 GB minimum; 64 GB+ recommended for scRNA-seq and MR steps
- **Disk**: ~54 GB for raw data + additional space for outputs
- All paths use `rprojroot` — portable across Windows / macOS / Linux

---

## Key R Packages

`GEOquery`, `limma`, `sva`, `WGCNA` (v1.72-1), `Seurat` (v4.4.0), `Monocle3`,  
`TwoSampleMR`, `MendelianRandomization`, `STRINGdb`, `clusterProfiler`,  
`ggplot2`, `ComplexHeatmap`, `patchwork`

---

## Citation

> **Zhou J, et al.** Lysosomal-Cytoskeletal Co-expression Modules Drive Rheumatoid Arthritis Pathogenesis: Causal Evidence from Multi-Omics and Mendelian Randomization. *(Journal & DOI to be updated upon publication)*

---

## Ethics & Data Reuse

All datasets used are publicly available (NCBI GEO, PRIDE, IEU OpenGWAS, BioBank Japan, UK Biobank) with prior ethical approvals by the original depositing institutions. No new human subjects data were collected.

## License

Released for academic and non-commercial research use only. Contact the corresponding author for other uses.

---

**Author**: Jinsong Zhou | **Version**: 1.0 | **Last updated**: March 2026
