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

##################  ################ ##################  ################
## initialize the data
## set a flag vector, to filter for qtl pairs that are properly tested in each cell type (start from a smallest cell type)
## some peaks are not existing because they are either SKIPPED in rasqual or did't converge during calculation
xx <- readRDS('~/yuzhao1/work/atac_gca2024/19rasqual/5collect_results/5results/allPeaks_noDiseaseCo/Neutrophil.rds')
qtlPair_all <- paste0(xx$Feature_ID, '@', xx$rsID) %>% unique()
qtlPair_all_flag <- setNames(rep(0, length(qtlPair_all)), qtlPair_all)
for (ct in celltypes) {
  cat(ct, '\n')
  ct.raw <- readRDS(paste0(dir.raw, ct, '.rds'))
  qtlPair_tmp <- paste0(ct.raw$Feature_ID, '@', ct.raw$rsID) %>% unique()
  qtlPair_all_flag[qtlPair_tmp] <- qtlPair_all_flag[qtlPair_tmp] + 1
}
qtlPair_all_celltype_tested <- unique(names(qtlPair_all_flag)[qtlPair_all_flag == 19])
qtlPair_all_celltype_tested <- qtlPair_all_celltype_tested[!is.na(qtlPair_all_celltype_tested)]


#################  ################ ##################  ################
## initialize the effect and sd matrix
qtlPair_all <- qtlPair_all_celltype_tested
df_effect <- data.frame(matrix(0, nrow = length(qtlPair_all), ncol = length(celltypes)))
df_sd <- data.frame(matrix(0, nrow = length(qtlPair_all), ncol = length(celltypes)))
rownames(df_effect) <- qtlPair_all
rownames(df_sd) <- qtlPair_all
colnames(df_effect) <- celltypes
colnames(df_sd) <- celltypes

for (ct in celltypes) {
  cat(ct, '\n')
  ct.raw <- readRDS(paste0(dir.raw, ct, '.rds'))
  ct.raw$peak_snp_pair <- paste0(ct.raw$Feature_ID, '@', ct.raw$rsID)
  ct.raw <- ct.raw[!duplicated(ct.raw$peak_snp_pair),]

  ct.raw$beta <- as.numeric(ct.raw$Effect_size) - 0.5
  df_effect[ct.raw$peak_snp_pair, ct] <- ct.raw$beta
  ct.raw$Chisquare <- as.numeric(ct.raw$Chisquare)
  df_sd[ct.raw$peak_snp_pair, ct] <- abs(ct.raw$beta/sqrt(ct.raw$Chisquare))
}

# remove NA
df_effect <- df_effect[complete.cases(df_effect), ]
df_sd <- df_sd[complete.cases(df_sd), ]

# remove INFINITE
df_effect <- df_effect[!is.infinite(rowSums(df_effect)),]
df_sd <- df_sd[!is.infinite(rowSums(df_sd)),]

# intersect to get those have both sd and beta values
qtlPair_all_celltype_tested_noNA <- intersect(rownames(df_effect), rownames(df_sd))
df_effect <- df_effect[qtlPair_all_celltype_tested_noNA, ]
df_sd <- df_sd[qtlPair_all_celltype_tested_noNA, ]

saveRDS(df_effect, paste0(out.dir, 'df_effect_qtlPair_all_celltype_tested.rds'))
saveRDS(df_sd, paste0(out.dir, 'df_sd_qtlPair_all_celltype_tested.rds'))
saveRDS(qtlPair_all_celltype_tested_noNA, paste0(out.dir, 'qtlPair_all_celltype_tested.rds'))


