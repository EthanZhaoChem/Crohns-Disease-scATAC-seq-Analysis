# https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html
# should test the suspicious locus first, use information in script 1 as well
library(bigsnpr)
library(mapgen)
library(susieR)
library(data.table)
library(dplyr)
library(stringr)
out.dir <- '~/yuzhao1/work/final_GCAatac/8susie/results/'
plots.dir <- '~/yuzhao1/work/final_GCAatac/8susie/diagnostic_plots/'

gwas_folder <- '~/yuzhao1/work/final_GCAatac/0gwas/b151_hg19/'
results_folder <- '~/yuzhao1/work/final_GCAatac/8susie/results/'
bigSNP <- snp_attach(rdsfile = '~/yuzhao1/resource/magma/g1000_eur.rds')
gwas <- readRDS('~/yuzhao1/work/final_GCAatac/8susie/results/gwas_processed.rds')
n = 40266
L <- 10

# 1. prepare
sig.loci <- gwas %>%
  group_by(locus) %>%
  summarise(max_mlogP = max(-log10(pval))) %>%
  filter(max_mlogP > -log10(5e-8)) %>% pull(locus)

gwas.sumstats.sigloci <- gwas[gwas$locus %in% sig.loci, ]

# 2. loop for each locus 
mismatch_rsID_list <- list()
for (locus in sig.loci) {
  cat('Calculating: \n')
  cat(locus, '\n')
  sumstats <- gwas.sumstats.sigloci
  sumstats <- sumstats[sumstats$locus == locus, ]
  X <- bigSNP$genotypes[, sumstats$bigSNP_index]
  X <- scale(X, center = T, scale = T)
  R <- cor(X)

  zflip <- sumstats$zscore
  ld <- R
  condz <- kriging_rss(zflip, ld, n=n)
  condz$conditional_dist$rsID <- sumstats$snp
  xx <- condz$conditional_dist
  
  locus_ID <- as.character(locus)
  mismatch_rsID_list[[locus_ID]] <- xx[which(abs(xx$z) > 2 & xx$logLR > 2), 'rsID']
  
  pdf(paste0(plots.dir, '', locus_ID, '.pdf'), width = 5, height = 5)
  print(condz$plot)
  dev.off()
}

saveRDS(mismatch_rsID_list, paste0(out.dir, 'diagnostic_mismatch_rsID_list.rds'))




