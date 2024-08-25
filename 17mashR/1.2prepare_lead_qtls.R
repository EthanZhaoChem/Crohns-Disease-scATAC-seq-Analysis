# get the effect size and sd for query QTLs
# missing values/ NA / infinite values will be set to a fixed values specified below
beta_default <- 0
sd_default <- 1e-6
q_default <- 1

dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(dplyr)
library(tidyverse)
library(data.table)
library(gtools)
library(limma)
library(ArchR)
library(ashr)
library(mashr)
"%&%" <- function(a, b) paste0(a, b)

##################  ################ ##################  ################
## presets
celltypes <- readLines('~/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt')
dir.raw <- '~/yuzhao1/work/atac_gca2024/19rasqual/5collect_results/5results/allPeaks_noDiseaseCo/'
out.dir <- '~/yuzhao1/work/atac_gca2024/26mash3/1select_qtl/'
qtlPair_all <- c()

for (ct in celltypes) {
  cat(ct, '\n')
  ct.raw <- readRDS(paste0(dir.raw, ct, '_lead.rds'))
  ct.raw$peak_snp_pair <- paste0(ct.raw$Feature_ID, '@', ct.raw$rsID)
  ct.raw <- ct.raw[!duplicated(ct.raw$peak_snp_pair),]
  qtlPair_all <- unique(c(ct.raw$peak_snp_pair, qtlPair_all))
}

qtlPair_all_flag <- setNames(rep(0, length(qtlPair_all)), qtlPair_all)

## initialize the effect and sd matrix, 
df_effect <- data.frame(matrix(beta_default, nrow = length(qtlPair_all), ncol = length(celltypes)))
df_sd <- data.frame(matrix(sd_default, nrow = length(qtlPair_all), ncol = length(celltypes)))
df_q <- data.frame(matrix(q_default, nrow = length(qtlPair_all), ncol = length(celltypes)))
rownames(df_effect) <- qtlPair_all
rownames(df_sd) <- qtlPair_all
rownames(df_q) <- qtlPair_all
colnames(df_effect) <- celltypes
colnames(df_sd) <- celltypes
colnames(df_q) <- celltypes


for (ct in celltypes) {
  cat(ct, '\n')
  ct.raw <- readRDS(paste0(dir.raw, ct, '.rds'))
  ct.raw$peak_snp_pair <- paste0(ct.raw$Feature_ID, '@', ct.raw$rsID)
  ct.raw <- ct.raw[!duplicated(ct.raw$peak_snp_pair),]
  ct.raw <- ct.raw[ct.raw$peak_snp_pair %in% qtlPair_all, ]
  
  ct.raw$beta <- as.numeric(ct.raw$Effect_size) - 0.5
  df_effect[ct.raw$peak_snp_pair, ct] <- ct.raw$beta
  ct.raw$Chisquare <- as.numeric(ct.raw$Chisquare)
  df_sd[ct.raw$peak_snp_pair, ct] <- abs(ct.raw$beta/sqrt(ct.raw$Chisquare))
  df_q[ct.raw$peak_snp_pair, ct] <- as.numeric(ct.raw$q)
  
  df_effect[is.na(df_effect[[ct]]), ct] <- beta_default
  df_effect[is.infinite(df_effect[[ct]]), ct] <- beta_default
  
  df_sd[is.na(df_sd[[ct]]), ct] <- sd_default
  df_sd[is.infinite(df_sd[[ct]]), ct] <- sd_default
}


saveRDS(df_effect, paste0(out.dir, 'df_effect_lead.rds'))
saveRDS(df_sd, paste0(out.dir, 'df_sd_lead.rds'))
saveRDS(df_q, paste0(out.dir, 'df_q_lead.rds'))


