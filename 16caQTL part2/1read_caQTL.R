dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')

library(stringr)
library(dplyr)
library(tidyverse)
library(data.table)
library(limma)
"%&%" <- function(a, b) paste0(a, b)

##################  ################ ##################  ################
# summarize the called snps, this is expaneded in foler7 in mashR
## presets
celltypes <- readLines('~/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt')
caQTL_dir <- '~/yuzhao1/work/atac_gca2024/19rasqual/5collect_results/5results/allPeaks_noDiseaseCo/'
filtered_peaks_dir <- '~/yuzhao1/work/atac_gca2024/19rasqual/4calculation_filtered_peaks/4ctSpecific_rasqualInput/'
out.dir <- '~/yuzhao1/work/atac_gca2024/24rasqual2/1rasqual_results/'
FDR_results_list <- list()

for(ct in celltypes){

  filename_rds <- paste0(caQTL_dir, ct, '_Perm_three_times_FDR_0.1.rds')
  FDR_results <- readRDS(filename_rds)
  
  # filter by genotype corr R2
  FDR_results <- FDR_results[FDR_results$Sq_corr_rSNP>=0.8, ]
  length(unique(FDR_results$Feature_ID))
  
  # filter peaks by insertion counts
  filtered_peaks <- read.columns(paste0(filtered_peaks_dir, ct, '.accessibility.txt'), required.col = 1)[[1]]
  FDR_results <- FDR_results[FDR_results$Feature_ID %in% filtered_peaks, ]
  length(unique(FDR_results$Feature_ID))
  
  # lead variants
  FDR_results <- FDR_results[order(FDR_results$P, decreasing = F),]
  FDR_results <- FDR_results[!duplicated(FDR_results$Feature_ID),]
  
  # save to list
  FDR_results_list[[ct]] <- FDR_results
}

saveRDS(FDR_results_list, out.dir %&% 'FDR_results_list.rds')

##################  ################ ##################  ################
# a few quick codes to read and format the data to use

# lead snp chr_pos_pos
FDR_results_list <- readRDS('~/yuzhao1/work/atac_gca2024/24rasqual2/1rasqual_results/FDR_results_list.rds')
qtls_lead_snps_list <- lapply(FDR_results_list, function(x) paste0(x$Chromosome, '_', x$SNP_position, '_', x$SNP_position))

# lead snp chr:pos:REF:ALT
rsID_all <- lapply(FDR_results_list, function(x) x$rsID) %>% unlist() %>% unique()

# peak - snp pair
xx <- lapply(FDR_results_list, function(x) paste0(x$FeatureID, '@', x$rsID)) %>% unlist() %>% unique()

# table
qtls_rasqual <- bind_rows(FDR_results_list, .id = 'celltype') %>% as.data.frame()



































