library(mapgen) # mapgen_0.5.9 susieR_0.12.35
library(tidyverse)
library(ggplot2)
library(data.table)

## 1. read necessary info 
results_folder <- '~/yuzhao1/work/final_GCAatac/8susie/results/'
n <- 40266
gwas_folder <- '~/yuzhao1/work/final_GCAatac/0gwas/b151_hg19/'
gwas_raw <- as.data.frame(
  fread(
    paste0(gwas_folder, "cd_build37_40266_20161107_hg19_b151_gwas_annotated_notFilterFor1KG.txt")
  )
)
LD_Blocks <- readRDS(system.file('extdata', 'LD.blocks.EUR.hg19.rds', package='mapgen'))
region_info <- get_UKBB_region_info(LD_Blocks,
                                    LDREF.dir = "/project2/mstephens/wcrouse/UKB_LDR_0.1_b37", 
                                    prefix = "ukb_b37_0.1")
LD_snp_info <- read_LD_SNP_info(region_info)

## 2. process gwas
gwas.sumstats <- process_gwas_sumstats(gwas_raw, 
                                       chr='chr', pos='pos', beta='beta', se='se', 
                                       a0='a0', a1='a1', snp='snp', pval='pval',
                                       LD_snp_info=LD_snp_info, 
                                       strand_flip=TRUE, 
                                       remove_strand_ambig=TRUE)
sig.loci <- gwas.sumstats %>% dplyr::filter(pval < 5e-8) %>% dplyr::pull(locus) %>% unique()
cat(length(sig.loci), "significant loci. \n")

## 3. run susie
selected.sumstats <- gwas.sumstats[gwas.sumstats$locus %in% sig.loci, ]
susie.results <- run_finemapping(selected.sumstats, 
                                 region_info = region_info, 
                                 priortype = 'uniform', 
                                 n = n, 
                                 L = 10,
                                 save = F,
                                 outname = "ukbb_L10")

saveRDS(selected.sumstats, paste0(results_folder, 'gwas.sumstats.sigloci.rds'))
saveRDS(susie.results, paste0(results_folder, 'cd_finemapping_unifprior_L10.rds'))

## 4, labeling HLA region
LD_Blocks <- read.csv('~/yuzhao1/work/final_GCAatac/8susie/LD_blocks_mapgen/LD_Blocks.csv')
gr_all <- GRanges(seqnames = LD_Blocks$chr, ranges = IRanges(start = LD_Blocks$start, end = LD_Blocks$end), locus_id = LD_Blocks$locus_id)
gr_hla <- GRanges(seqnames = 'chr6', ranges = IRanges(start = 25e6, end = 35e6))
overlaps <- GenomicRanges::findOverlaps(gr_all, gr_hla, ignore.strand = T)
overlapped_loci <- gr_all[overlaps@from, ]
write.csv(overlapped_loci, '~/yuzhao1/work/final_GCAatac/8susie/LD_blocks_mapgen/LD_Blocks_HLA.csv')






