library(fastTopics)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(plyr)
library(dplyr)
library(stringr)
library(ArchR)


plot.dir <- '~/yuzhao1/work/atac_gca2024/13fasttopic/plots/'
data.dir <- '~/yuzhao1/work/atac_gca2024/13fasttopic/rds/'
file_fasttopic <- '~/yuzhao1/work/atac_gca2024/13fasttopic/rds/fit_union_sub100_k45_converged.rds'
file_de <- '~/yuzhao1/work/atac_gca2024/13fasttopic/rds/fit_union_sub100_k45_converged_de_vsnull.rds'

################################## run ##############################################
load("~/yuzhao1/work/atac_gca2024/13fasttopic/rds/union_raw_sub100.RData")
fit <- readRDS(file_fasttopic)
fit <- poisson2multinom(fit)


# de_analysis
set.seed(6)
# de <- de_analysis(fit,counts,pseudocount = 0.01, shrink.method = "none",
#                   control = list(ns = 1e4,nc = 12))
de <- de_analysis(fit,counts,pseudocount = 0.01, shrink.method = "none", lfc.stat = 'vsnull',
                  control = list(ns = 1e4,nc = 12))

saveRDS(de, file_de)



