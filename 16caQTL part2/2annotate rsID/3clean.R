library(ggplot2)
library(dplyr)
library(plyr)
library(stringr)

snps_raw <- readLines('~/yuzhao1/work/atac_gca2024/24rasqual2/2annotateRSid/snps_raw.txt')
snp_annotation <- read.table('~/yuzhao1/work/atac_gca2024/24rasqual2/2annotateRSid/snps_annotated.txt')
snp_annotation <- snp_annotation[, c(1,2,3)]
colnames(snp_annotation) <- c('chr', 'pos', 'rsID')

snp_annotation[!grepl('rs', snp_annotation$rsID), 'rsID'] <- paste0(snp_annotation[!grepl('rs', snp_annotation$rsID), c('chr')],
                                                                    ':',
                                                                    snp_annotation[!grepl('rs', snp_annotation$rsID), c('pos')])

write.csv(snp_annotation, '~/yuzhao1/work/atac_gca2024/24rasqual2/2annotateRSid/snps_rsID_clean.txt', row.names = F, quote = F)
