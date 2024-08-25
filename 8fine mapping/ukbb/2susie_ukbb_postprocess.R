library(bigsnpr)
library(mapgen)
library(susieR)
library(data.table)
library(dplyr)
library(stringr)
library(Repitools)
results_folder <- '~/yuzhao1/work/final_GCAatac/8susie/results/'


## function to process susie raw results
merge_susie_sumstats_ethan <- function (susie_results, sumstats) 
{
  sumstats$susie_pip <- 0
  sumstats$CS <- '0'
  loci <- names(susie_results)
  for (l in loci) {
    n.snps <- length(susie_results[[l]]$pip)
    sumstats[sumstats$locus == as.numeric(l), "susie_pip"] <- susie_results[[l]]$pip
    snps.in.cs <- rep('0', n.snps)
    if (!is.null(susie_results[[l]]$sets$cs)) {
      csNames <- names(susie_results[[l]]$sets$cs)
      for (csName in csNames) {
        snps.in.cs[unlist(susie_results[[l]]$sets$cs[[csName]])] <- csName
      }
    }
    sumstats[sumstats$locus == as.numeric(l), "CS"] <- snps.in.cs
  }
  return(sumstats)
}

########################### 1. post process susie results ####################
gwas.sumstats.sigloci <- readRDS('~/yuzhao1/work/final_GCAatac/8susie/results/gwas.sumstats.sigloci.rds')
susie_raw_L10 <- readRDS(paste0(results_folder, 'cd_finemapping_unifprior_L10.rds'))
susie_merged_L10 <- merge_susie_sumstats_ethan(susie_results = susie_raw_L10, sumstats = gwas.sumstats.sigloci)
susie_clean_L10 <- susie_merged_L10[susie_merged_L10$CS!='0',]
susie_clean_L10 <- as.data.frame(susie_clean_L10)
susie_clean_L10 <- susie_clean_L10[susie_clean_L10$snp!='.', ]
susie_clean_L10$signal_ID <- paste0(susie_clean_L10$locus, '_', susie_clean_L10$CS)
rownames(susie_clean_L10) <- susie_clean_L10$snp

# statistics of susie results
df <- susie_clean_L10
df$nSNPs_perCS <- 0
for (i in 1:nrow(df)) {
  df$nSNPs_perCS[[i]] <- sum(df$signal_ID==df[i, 'signal_ID'])
}

# save post processed result
write.csv(df, paste0(results_folder, 'cd_finemapping_unifprior_L10.csv'))

table(table(df$signal_ID)) # nSNPs per signal
table(table(strsplit(unique(df$signal_ID), split = '_') %>% sapply(.,`[[`,1))) # nSignals per locus

# write bed file for finemapped SNPs
finemapped_snps_L10_bed <- data.frame(chr = paste0('chr', susie_clean_L10$chr),
                                      start = susie_clean_L10$pos,
                                      end = susie_clean_L10$pos)
write.table(finemapped_snps_L10_bed, paste0(results_folder, 'finemapped_snps_L10.bed'), 
            row.names = F, col.names = F, sep="\t", quote=FALSE)



