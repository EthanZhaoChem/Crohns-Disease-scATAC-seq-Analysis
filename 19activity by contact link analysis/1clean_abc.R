dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(plyr)
library(dplyr)
library(ArchR)
source('~/yuzhao1/work/atac_gca2024/22abc/helper_abc.R')
out.dir <- '~/yuzhao1/work/atac_gca2024/22abc/'

### read abc-biosample enrichment
df_biosample <- read.table('~/yuzhao1/work/atac_gca2024/22abc/gwas_biosample_enrichment.txt', header = T)
biosamples_ibd <- unique(df_biosample[df_biosample$trait=='IBD' & df_biosample$FM.Enriched==T, 'Biosample'])

### read abc data (already mapped into hg38)
abc.raw <- read.table('~/yuzhao1/resource/abc/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3_hg38.txt', header = T)
cts.extract <- c("small_intestine_fetal-Roadmap", biosamples_ibd)
columns.extract <- c("chr", "start", "end", "name", "class", "activity_base", "TargetGene", "TargetGeneTSS", "ABC.Score", "CellType" )
abc.clean <- abc.raw[abc.raw$CellType %in% cts.extract, columns.extract]

### use this format
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1/")
links <- getPeak2GeneLinks(proj, corCutOff = 0.5)

### use archr TSS location
genes_all <- getGeneAnnotation(proj)
genes_all <- genes_all$genes
df_genes <- genes_all[!is.na(genes_all$symbol)]

abc.clean$TargetGeneTSS <- mapvalues(abc.clean$TargetGene, df_genes$symbol, start(df_genes), warn_missing = F)
abc.clean$TargetGeneTSS <- as.numeric(abc.clean$TargetGeneTSS)
abc.clean$TargetGeneTSS[is.na(abc.clean$TargetGeneTSS)] <- 0

### peak mid-point to TSS (min to max)
gr_abc <- GRanges(seqnames = abc.clean$chr, 
                  ranges = IRanges(start = pmin((abc.clean$start + abc.clean$end)/2, abc.clean$TargetGeneTSS),
                                   end = pmax((abc.clean$start + abc.clean$end)/2, abc.clean$TargetGeneTSS)),
                  value=abc.clean$ABC.Score,
                  FDR = 0)
links$Peak2GeneLinks <- gr_abc
saveRDS(links, '~/yuzhao1/work/atac_gca2024/22abc/abc_links_ibd.rds')
saveRDS(abc.clean, '~/yuzhao1/work/atac_gca2024/22abc/abc_df_ibd.rds')
