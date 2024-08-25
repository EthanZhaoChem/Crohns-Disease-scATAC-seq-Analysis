library(fastTopics)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(plyr)
library(dplyr)
library(stringr)

source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/scripts/helper_archr.R')

file_fasttopic <- '~/yuzhao1/work/atac_gca2024/13fasttopic/rds/fit_union_sub100_k50_750iterations.rds'


################################## run ##############################################
fit <- readRDS(file_fasttopic)
fit <- poisson2multinom(fit)

plot_progress(fit,x = "iter",add.point.every = 5,colors = "black") +
  theme_cowplot(font_size = 10)

