# MANUAL DATA DOWNLOAD & EXTRACTION INSTRUCTIONS

All large raw data files (total ~54 GB) must be downloaded and extracted **manually** before running the pipeline.  
This avoids any automatic download issues.

### Step-by-step instructions
1. Download all files using your browser or **Motrix** (recommended for speed).
2. Save each file with the **exact filename** shown.
3. For the 4 GSE*.tar files: manually extract them into the correct subfolders.
4. After everything is ready, run `data/check_data_ready.R` to verify.
5. Then you can run the scripts in `scripts/` from 101 onwards.

**Total files: 29** (6 STRINGdb + 4 GEO tar + 1 RDS + 4 feather + 14 GWAS)

### 1. STRINGdb files (6 files) – save to `data/`
| Filename                              | Direct Download Link |
|---------------------------------------|----------------------|
| 9606.protein.aliases.v11.5.txt.gz     | https://stringdb-static.org/download/protein.aliases.v11.5/9606.protein.aliases.v11.5.txt.gz |
| 9606.protein.aliases.v12.0.txt.gz     | https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz |
| 9606.protein.info.v11.5.txt.gz        | https://stringdb-static.org/download/protein.info.v11.5/9606.protein.info.v11.5.txt.gz |
| 9606.protein.info.v12.0.txt.gz        | https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz |
| 9606.protein.links.v11.5.txt.gz       | https://stringdb-static.org/download/protein.links.v11.5/9606.protein.links.v11.5.txt.gz |
| 9606.protein.links.v12.0.txt.gz       | https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz |

### 2. GEO & cisTarget files (9 files) – save to `data/`
| Filename (save as)                                      | Direct Download Link |
|---------------------------------------------------------|----------------------|
| GSE159117_RAW.tar                                       | https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE159117&format=file |
| GSE55235_RAW.tar                                        | https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE55235&format=file |
| GSE55457_RAW.tar                                        | https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE55457&format=file |
| GSE55584_RAW.tar                                        | https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE55584&format=file |
| GSE296117_RA_geo.rds                                    | https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE296117&format=file&file=GSE296117%5FRA%5Fgeo%2Erds |
| hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather | https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather |
| hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather | https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather |
| hg38_screen_v10_clust.regions_vs_motifs.rankings.feather | https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather |
| hg38_screen_v10_clust.regions_vs_motifs.scores.feather  | https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.scores.feather |

### 3. GWAS Summary Statistics (14 files) – NEW – save to `data/GWAS/`
**Important**: First create the folder `data/GWAS/`.

| Filename (save as)          | Direct Download Link |
|-----------------------------|----------------------|
| bbj-a-151.vcf.gz            | https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/0EFJ1CMPyZH2uBqL2i1QT7X87hYbk9PBBY_m7bXy0PEAUBY1yUeH_YKjyjEYTVA1/n/ieup4/b/igd/o/bbj-a-151/bbj-a-151.vcf.gz |
| bbj-a-72.vcf.gz             | https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/XsldeQ6un94AnOSVSJ6TzPLTaqyG-csxmiPFI8yHGncbPv1m5Y7lCXpgss0LPkI1/n/ieup4/b/igd/o/bbj-a-72/bbj-a-72.vcf.gz |
| bbj-a-73.vcf.gz             | https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/ArGq2UvJnf2mzz8JlZVUbfsv6eTQiL9AK4J1l_I15Tmc-yQgqhXkJXu2l3aWTRq1/n/ieup4/b/igd/o/bbj-a-73/bbj-a-73.vcf.gz |
| bbj-a-74.vcf.gz             | https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/kQkdHF3NuNwWo1KkHPpVck6hWQA4E5rbHqRX6KC1N9LFiEiBJl5cbwlxwhoLlAGp/n/ieup4/b/igd/o/bbj-a-74/bbj-a-74.vcf.gz |
| ieu-a-831.vcf.gz            | https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/_zzfSjaW9UuAvw4SQZZJyD_cshO4Ubch12tdgVKDTnPbPOwdkgBD-sdAajS-kvHP/n/ieup4/b/igd/o/ieu-a-831/ieu-a-831.vcf.gz |
| ieu-a-832.vcf.gz            | https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/FGlx_42Fke4FfxmVbi4AsXCOD-10Ipg4kogZTgzLK20Ztbyc2-IS5wDCyWZHXikA/n/ieup4/b/igd/o/ieu-a-832/ieu-a-832.vcf.gz |
| ieu-a-833.vcf.gz            | https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/fTLiQyVfkaGwJZeM39XBagV2SkGB-LoTdHsa9D1_YJZAhz-2DpZaJh0eab2CZEy7/n/ieup4/b/igd/o/ieu-a-833/ieu-a-833.vcf.gz |
| ieu-a-834.vcf.gz            | https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/AL-ZbfB8F6y3jFVGxW2FFGR_EIP2TK4O82U1ObhA_ikKQCOg0S7WD4b4QHkIp95J/n/ieup4/b/igd/o/ieu-a-834/ieu-a-834.vcf.gz |
| ukb-a-105.vcf.gz            | https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/hLQdd9PJmFOYyc1C_swxK-z2-ByBrrcjBUVyqM0nagh-ttPka5StED6Lzmhq24M6/n/ieup4/b/igd/o/ukb-a-105/ukb-a-105.vcf.gz |
| ukb-b-11874.vcf.gz          | https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/y8yKaj9rZZV_jZNgCbYvADl5XjVSy13HUJUaUPMCAwvTkygfxWdbdypRHStVUpE9/n/ieup4/b/igd/o/ukb-b-11874/ukb-b-11874.vcf.gz |
| ukb-b-9125.vcf.gz           | https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/YeoUc0yjpJ7yS-b1RgncK9fCMMb9BPDT-pnUO8I2at_s35UD2fQIFWEfTDOCMOzc/n/ieup4/b/igd/o/ukb-b-9125/ukb-b-9125.vcf.gz |
| ukb-d-M06.vcf.gz            | https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/Uz5KJWcGORVz78rEB18VB8NvAVmQHYeXTVrdWR3hq_xATRfT8LgkNB1nMY2juGhw/n/ieup4/b/igd/o/ukb-d-M06/ukb-d-M06.vcf.gz |
| ukb-d-M13_RHEUMA.vcf.gz     | https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/Q6amxwXG18U4bB9sRkA5H3BdPSxPZBJVQ2WhsEW2Hr7G7tZ4iD3-U2XdE5mRsdq-/n/ieup4/b/igd/o/ukb-d-M13_RHEUMA/ukb-d-M13_RHEUMA.vcf.gz |
| ukb-d-RHEUMA_NOS.vcf.gz     | https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/7MaNoSLBPtoh7Rgg4DavfLrrTyONRIfE-jyc9qSdKWAx1k8r7QbJpLL0FGgRVnYm/n/ieup4/b/igd/o/ukb-d-RHEUMA_NOS/ukb-d-RHEUMA_NOS.vcf.gz |

### Motrix Batch Download (all 25 files)
Copy & paste the entire block below into Motrix → New Task → Batch:

https://stringdb-static.org/download/protein.aliases.v11.5/9606.protein.aliases.v11.5.txt.gz 9606.protein.aliases.v11.5.txt.gz
https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz 9606.protein.aliases.v12.0.txt.gz
https://stringdb-static.org/download/protein.info.v11.5/9606.protein.info.v11.5.txt.gz 9606.protein.info.v11.5.txt.gz
https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz 9606.protein.info.v12.0.txt.gz
https://stringdb-static.org/download/protein.links.v11.5/9606.protein.links.v11.5.txt.gz 9606.protein.links.v11.5.txt.gz
https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz 9606.protein.links.v12.0.txt.gz
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE159117&format=file GSE159117_RAW.tar
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE55235&format=file GSE55235_RAW.tar
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE55457&format=file GSE55457_RAW.tar
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE55584&format=file GSE55584_RAW.tar
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE296117&format=file&file=GSE296117%5FRA%5Fgeo%2Erds GSE296117_RA_geo.rds
https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather
https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather hg38_screen_v10_clust.regions_vs_motifs.rankings.feather
https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.scores.feather hg38_screen_v10_clust.regions_vs_motifs.scores.feather
https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/0EFJ1CMPyZH2uBqL2i1QT7X87hYbk9PBBY_m7bXy0PEAUBY1yUeH_YKjyjEYTVA1/n/ieup4/b/igd/o/bbj-a-151/bbj-a-151.vcf.gz bbj-a-151.vcf.gz
https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/XsldeQ6un94AnOSVSJ6TzPLTaqyG-csxmiPFI8yHGncbPv1m5Y7lCXpgss0LPkI1/n/ieup4/b/igd/o/bbj-a-72/bbj-a-72.vcf.gz bbj-a-72.vcf.gz
https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/ArGq2UvJnf2mzz8JlZVUbfsv6eTQiL9AK4J1l_I15Tmc-yQgqhXkJXu2l3aWTRq1/n/ieup4/b/igd/o/bbj-a-73/bbj-a-73.vcf.gz bbj-a-73.vcf.gz
https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/kQkdHF3NuNwWo1KkHPpVck6hWQA4E5rbHqRX6KC1N9LFiEiBJl5cbwlxwhoLlAGp/n/ieup4/b/igd/o/bbj-a-74/bbj-a-74.vcf.gz bbj-a-74.vcf.gz
https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/_zzfSjaW9UuAvw4SQZZJyD_cshO4Ubch12tdgVKDTnPbPOwdkgBD-sdAajS-kvHP/n/ieup4/b/igd/o/ieu-a-831/ieu-a-831.vcf.gz ieu-a-831.vcf.gz
https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/FGlx_42Fke4FfxmVbi4AsXCOD-10Ipg4kogZTgzLK20Ztbyc2-IS5wDCyWZHXikA/n/ieup4/b/igd/o/ieu-a-832/ieu-a-832.vcf.gz ieu-a-832.vcf.gz
https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/fTLiQyVfkaGwJZeM39XBagV2SkGB-LoTdHsa9D1_YJZAhz-2DpZaJh0eab2CZEy7/n/ieup4/b/igd/o/ieu-a-833/ieu-a-833.vcf.gz ieu-a-833.vcf.gz
https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/AL-ZbfB8F6y3jFVGxW2FFGR_EIP2TK4O82U1ObhA_ikKQCOg0S7WD4b4QHkIp95J/n/ieup4/b/igd/o/ieu-a-834/ieu-a-834.vcf.gz ieu-a-834.vcf.gz
https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/hLQdd9PJmFOYyc1C_swxK-z2-ByBrrcjBUVyqM0nagh-ttPka5StED6Lzmhq24M6/n/ieup4/b/igd/o/ukb-a-105/ukb-a-105.vcf.gz ukb-a-105.vcf.gz
https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/y8yKaj9rZZV_jZNgCbYvADl5XjVSy13HUJUaUPMCAwvTkygfxWdbdypRHStVUpE9/n/ieup4/b/igd/o/ukb-b-11874/ukb-b-11874.vcf.gz ukb-b-11874.vcf.gz
https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/YeoUc0yjpJ7yS-b1RgncK9fCMMb9BPDT-pnUO8I2at_s35UD2fQIFWEfTDOCMOzc/n/ieup4/b/igd/o/ukb-b-9125/ukb-b-9125.vcf.gz ukb-b-9125.vcf.gz
https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/Uz5KJWcGORVz78rEB18VB8NvAVmQHYeXTVrdWR3hq_xATRfT8LgkNB1nMY2juGhw/n/ieup4/b/igd/o/ukb-d-M06/ukb-d-M06.vcf.gz ukb-d-M06.vcf.gz
https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/Q6amxwXG18U4bB9sRkA5H3BdPSxPZBJVQ2WhsEW2Hr7G7tZ4iD3-U2XdE5mRsdq-/n/ieup4/b/igd/o/ukb-d-M13_RHEUMA/ukb-d-M13_RHEUMA.vcf.gz ukb-d-M13_RHEUMA.vcf.gz
https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/7MaNoSLBPtoh7Rgg4DavfLrrTyONRIfE-jyc9qSdKWAx1k8r7QbJpLL0FGgRVnYm/n/ieup4/b/igd/o/ukb-d-RHEUMA_NOS/ukb-d-RHEUMA_NOS.vcf.gz ukb-d-RHEUMA_NOS.vcf.gz


**Extraction tip**: Use 7-Zip (Windows) or built-in macOS Archive Utility. Extract each .tar file **directly** into its corresponding folder (e.g. extract `GSE159117_RAW.tar` → `data/GSE159117_RAW/`).

**After finishing**, run `data/check_data_ready.R`.

