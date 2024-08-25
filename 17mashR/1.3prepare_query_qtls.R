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
df_q <- readRDS(paste0(out.dir, 'df_q_lead.rds'))
df_effect <- readRDS(paste0(out.dir, 'df_effect_lead.rds'))
df_sd <- readRDS(paste0(out.dir, 'df_sd_lead.rds'))

##################  ################ ##################  ################
### prioritize one snp for one peak
df_q$peak <- rownames(df_q) %>% strsplit(split = '@', fixed=T) %>% sapply(.,`[[`,1)
df_q$snp <- rownames(df_q) %>% strsplit(split = '@', fixed=T) %>% sapply(.,`[[`,2)
df_q$min_q <- apply(df_q[, celltypes], 1, min)
df_q <- df_q[order(df_q$min_q, decreasing=F),]
df_q <- df_q[!duplicated(df_q$peak),]
qtlPair_all <- rownames(df_q)

df_effect <- df_effect[qtlPair_all,]
df_sd <- df_sd[qtlPair_all,]
saveRDS(df_effect, paste0(out.dir, 'df_effect_query.rds'))
saveRDS(df_sd, paste0(out.dir, 'df_sd_query.rds'))





















