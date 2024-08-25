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
source('~/yuzhao1/scripts/plot.R')
"%&%" <- function(a, b) paste0(a, b)

##################  ################ ##################  ################
## presets
celltypes <- readLines('~/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt')
dir.raw <- '~/yuzhao1/work/atac_gca2024/19rasqual/5collect_results/5results/allPeaks_noDiseaseCo/'
prep.dir <- '~/yuzhao1/work/atac_gca2024/26mash3/1select_qtl/'
out.dir <- '~/yuzhao1/work/atac_gca2024/26mash3/2mash/'

## read significant QTLs
FDR_results_list <- readRDS('~/yuzhao1/work/atac_gca2024/24rasqual2/1rasqual_results/FDR_results_list.rds')
qtls_rasqual <- bind_rows(FDR_results_list, .id = 'celltype') %>% as.data.frame()
qtls_rasqual$peak_snp_pair <- paste0(qtls_rasqual$Feature_ID, '@', qtls_rasqual$rsID)

df_effect <- readRDS(paste0(prep.dir, 'df_effect_qtlPair_all_celltype_tested.rds'))
df_sd <- readRDS(paste0(prep.dir, 'df_sd_qtlPair_all_celltype_tested.rds'))
qtlPair_all_celltype_tested_noNA <- readRDS(paste0(prep.dir, 'qtlPair_all_celltype_tested.rds'))


##################  ################ ##################  ################
# standard error minumum set to 1e-6 to avoid numerical problems, only fewer than 100/18million qtls have this problem
zero_threshold <- 1e-6 # for minimum standard error
for (ct in celltypes) {
  cat(ct, '\n')
  xx <- which(df_sd[, ct] < zero_threshold)
  df_sd[xx, ct] <- zero_threshold
}

##################  ################ ##################  ################
## select 2million random QTLs, and strong QTLs that are passed FDR
qtls_pool <- qtlPair_all_celltype_tested_noNA

random.subset <- sample(qtls_pool, 2e6)
strong.subset <- intersect(qtls_rasqual$peak_snp_pair, qtls_pool) %>% unique()

simdata = list(Bhat=as.matrix(df_effect), Shat=as.matrix(df_sd))

data.temp = mash_set_data(simdata$Bhat[random.subset,],simdata$Shat[random.subset,])
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)

data.random = mash_set_data(simdata$Bhat[random.subset,],simdata$Shat[random.subset,],V=Vhat)
data.strong = mash_set_data(simdata$Bhat[strong.subset,],simdata$Shat[strong.subset,], V=Vhat)

U.pca = cov_pca(data.strong, 5)
U.ed = cov_ed(data.strong, U.pca)

U.c = cov_canonical(data.random)
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)

saveRDS(m, paste0(out.dir, 'mash_m_model_random_2million.rds'))
















