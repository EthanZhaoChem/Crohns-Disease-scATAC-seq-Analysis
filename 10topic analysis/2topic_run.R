args <- commandArgs(trailingOnly = TRUE) 
lineage <- args[1]
nTopic <- as.numeric(args[2])
nIteration1 <- as.numeric(args[3])
nIteration2 <- as.numeric(args[4])
flag_continue <- as.numeric(args[5])
iterations_before <- as.numeric(args[6])
iterations_add <- as.numeric(args[7])

if(flag_continue == 0){
  cat('nIteration1:', nIteration1, '\n')
  cat('nIteration2:', nIteration2, '\n')
  nIteration_aiming <- nIteration1 + nIteration2
}

if(flag_continue == 1){
  cat('nIterations before:', iterations_before, '\n')
  cat('nIterations add:', iterations_add, '\n')
  nIteration_aiming <- iterations_before + iterations_add
}

library(fastTopics)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(plyr)
library(dplyr)
library(stringr)
library(ArchR)

source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/scripts/helper_archr.R')


plot.dir <- '~/yuzhao1/work/atac_gca2024/13fasttopic/plots/'
data.dir <- '~/yuzhao1/work/atac_gca2024/13fasttopic/rds/'
load(paste0('~/yuzhao1/work/atac_gca2024/13fasttopic/rds/', lineage, '_raw_sub100.RData'))


################################## run ##############################################
cat(paste0('fitting ', lineage, ' lineage with N topics: ', nTopic, ' for ', nIteration_aiming, ' iterations\n'))
set.seed(6)
t0 <- proc.time()

if(flag_continue == 0){
  fit0 <- fit_poisson_nmf(counts,k = nTopic,numiter = nIteration1, method = "em",
                        init.method = "random",verbose = "detailed",
                        control = list(numiter = 4,extrapolate = FALSE,
                                       nc = 24))
  
  fit <- fit_poisson_nmf(counts,fit0 = fit0,numiter = nIteration2, method = "scd",
                         control = list(numiter = 4,extrapolate = TRUE,nc = 24),
                         verbose = "detailed")
}

if(flag_continue == 1){
  file_fasttopic0 <- paste0('~/yuzhao1/work/atac_gca2024/13fasttopic/rds/fit_', lineage, '_sub100_k', nTopic, '_', iterations_before, 'iterations.rds')
  fit0 <- readRDS(file_fasttopic0)
  fit <- fit_poisson_nmf(counts,fit0 = fit0,numiter = iterations_add, method = "scd",
                         control = list(numiter = 4,extrapolate = TRUE,nc = 24),
                         verbose = "detailed")
}

t1 <- proc.time()
print(t1 - t0)

cat(paste0('fitting completed with N topics: ', nTopic))
saveRDS(fit, paste0('~/yuzhao1/work/atac_gca2024/13fasttopic/rds/fit_', lineage, '_sub100_k', nTopic, '_', nIteration_aiming, 'iterations.rds'))








