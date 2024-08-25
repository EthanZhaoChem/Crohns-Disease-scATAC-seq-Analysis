library(ggplot2)
library(cowplot)
library(tidyverse)
library(plyr)
library(dplyr)
library(stringr)
library(GenomicRanges)


# need to modify
bedfilesNewDir <- '/home/yuzhao1/yuzhao1/work/atac_gca2024/14ldsc/results/union_sub100_k45_daPeaks_positive_flexibleLpval30k_vsnull/bed_lifted_extended/'
bedfiles_path <- list.files(bedfilesNewDir)

for (i in 1:length(bedfiles_path)){
  df <- read.table(paste0(bedfilesNewDir, bedfiles_path[[i]]))
  
  ## sort the peaks
  df <- df[df[,1] %in% c(paste0('chr', 1:22), 'chrX'), ]
  df$chr_sorted <- df[, 1]
  df[df[, 1]=='chrX', 'chr_sorted'] <- 'chr998' # to help sort chrX
  df$chr_sorted <- as.numeric(gsub("chr", "", df[, 'chr_sorted']))
  df <- df[order(df$chr_sorted, df[, 2]), ]
  df <- df[, 1:3] # Remove the auxiliary column
  colnames(df) <- c('chr', 'start', 'end')
  
  # # merge overlapped peaks
  # # Create a GRanges object
  # gr <- GRanges(seqnames=df$chr, ranges=IRanges(start=df$start, end=df$end))
  # 
  # # Merge overlapping intervals
  # merged <- reduce(gr)
  # 
  # # Convert back to a dataframe
  # merged_df <- data.frame(chr = seqnames(merged), start = start(merged), end = end(merged))
  # 
  
  write.table(df, paste0(bedfilesNewDir, bedfiles_path[[i]]), 
              row.names = F,
              col.names = F, 
              sep="\t", 
              quote=FALSE)
}
























                            
                                