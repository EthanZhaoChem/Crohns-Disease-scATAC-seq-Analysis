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

## read query data
df_effect <- readRDS(paste0(prep.dir, 'df_effect_query.rds'))
df_sd <- readRDS(paste0(prep.dir, 'df_sd_query.rds'))

## read model
m <- readRDS(paste0(out.dir, 'mash_m_model_random_2million.rds'))

##################  ################ ##################  ################
# standard error minumum set to 1e06 to avoid numerical problems, only fewer than 100/18million qtls have this problem
zero_threshold <- 1e-6 # for minimum standard error
for (ct in celltypes) {
  cat(ct, '\n')
  xx <- which(df_sd[, ct] < zero_threshold)
  df_sd[xx, ct] <- zero_threshold
}

##################  ################ ##################  ################
## get all the information for QTLs to be tested
simdata = list(Bhat=as.matrix(df_effect), Shat=as.matrix(df_sd))
data.query = mash_set_data(simdata$Bhat,simdata$Shat)


##################  ################ ##################  ################
m2 = mash(data.query, g=get_fitted_g(m), fixg=TRUE)
df_lfsr <- get_lfsr(m2)

saveRDS(m2, paste0(out.dir, 'mash_m2_results_random_2million.rds'))
saveRDS(df_lfsr, paste0(out.dir, 'df_lfsr.rds'))
















