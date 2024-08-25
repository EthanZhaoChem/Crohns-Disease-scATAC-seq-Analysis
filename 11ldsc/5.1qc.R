library(ggplot2)
library(cowplot)
library(tidyverse)
library(plyr)
library(dplyr)
library(stringr)
library(rmeta)
library(data.table)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/work/final_GCArna/scripts/gca_colors.R')

########################### part 1 ##############################
# check parameter M: all snps used for heritability analysis
Mref <- 0
for (chr in 1:22){
  snps_tmp <- readLines(paste0('~/spott/yuzhao1/ldsc/minimum_required/1000G_Phase3_ldscores/LDscore.', chr, '.l2.M_5_50'))
  snps_tmp <- as.numeric(snps_tmp)
  Mref <- Mref + snps_tmp
}
# after calculation, it is 5961159


########################### part 2 ##############################
annot_cell <- '~/yuzhao1/work/atac_gca2024/14ldsc/results/peaks_anno1_5kmin_6TSS_DoubletRatio2_filtered1/annots/'
annot_names <- list.dirs(annot_cell, recursive = F, full.names = F)
Nsnps_annotated <- matrix(0, nrow = length(annot_names), ncol = 2)
rownames(Nsnps_annotated) <- annot_names
colnames(Nsnps_annotated) <- c('Nsnps', 'ratio')

for(aa in 1:length(annot_names)){
  cell_path = paste0(annot_cell, annot_names[aa], "/")
  ll <- list.files(cell_path, pattern = ".annot.gz")
  num <- 0
  
  cat('working on ', annot_names[aa], '\n')
  for(m in 1:length(ll)){
    dat <- data.frame(fread(cmd=paste0("zcat ", cell_path, ll[m])))
    num = num  + sum(dat$ANNOT)
    rm(dat)
  }
  Nsnps_annotated[aa, 1]<- num
  Nsnps_annotated[aa, 2]<- num/Mref
}


saveRDS(Nsnps_annotated, '~/yuzhao1/work/atac_gca2024/14ldsc/QC/Nsnps_annotated_anno1Group.rds')


xx <- readRDS('~/yuzhao1/work/atac_gca2024/14ldsc/QC/Nsnps_annotated_anno1Group.rds')




