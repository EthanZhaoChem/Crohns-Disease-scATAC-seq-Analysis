library(ggplot2)
library(dplyr)
library(plyr)
library(stringr)

FDR_results_list <- readRDS('~/yuzhao1/work/atac_gca2024/24rasqual2/1rasqual_results/FDR_results_list.rds')
rsID_all <- lapply(FDR_results_list, function(x) x$rsID) %>% unlist() %>% unique()
writeLines(rsID_all, '~/yuzhao1/work/atac_gca2024/24rasqual2/2annotateRSid/snps_raw.txt')




