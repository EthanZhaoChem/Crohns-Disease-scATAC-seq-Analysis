library(bigsnpr)
library(mapgen)
library(susieR)
library(data.table)
library(dplyr)
library(stringr)
library(Repitools)


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
gwas.sumstats.sigloci <- readRDS('~/yuzhao1/work/final_GCAatac/8susie/results/gwas.sumstats.sigloci_95loci.rds')
susie_raw_L10 <- readRDS('~/yuzhao1/work/final_GCAatac/8susie/results/cd_finemapping_unifprior_95loci_L10.rds')
susie_merged_L10 <- merge_susie_sumstats_ethan(susie_results = susie_raw_L10, sumstats = gwas.sumstats.sigloci)
susie_clean_L10 <- susie_merged_L10[susie_merged_L10$CS!='0',]
susie_clean_L10 <- as.data.frame(susie_clean_L10)
susie_clean_L10 <- susie_clean_L10[susie_clean_L10$snp!='.', ]
susie_clean_L10$signal_ID <- paste0(susie_clean_L10$locus, '_', susie_clean_L10$CS)
rownames(susie_clean_L10) <- susie_clean_L10$snp
leadingSNPSs_susieL10 <- susie_clean_L10[susie_clean_L10 $susie_pip > 0.5, 'snp']

# statistics of susie results
df <- susie_clean_L10
df$nSNPs_perCS <- 0
for (i in 1:nrow(df)) {
  df$nSNPs_perCS[[i]] <- sum(df$signal_ID==df[i, 'signal_ID'])
}
df$nSNPs_perCS <- sum(df$signal_ID==df)

# save post processed result
write.csv(df, '~/yuzhao1/work/final_GCAatac/8susie/results/cd_finemapping_unifprior_95loci_L10.csv')

table(table(df$signal_ID)) # nSNPs per signal
table(table(strsplit(unique(df$signal_ID), split = '_') %>% sapply(.,`[[`,1))) # nSignals per locus

# write bed file for finemapped SNPs
finemapped_snps_L10_bed <- data.frame(chr = paste0('chr', susie_clean_L10$chr),
                                      start = susie_clean_L10$pos,
                                      end = susie_clean_L10$pos)
write.table(finemapped_snps_L10_bed, '~/yuzhao1/work/final_GCAatac/8susie/results/finemapped_snps_L10.bed', 
            row.names = F, col.names = F, sep="\t", quote=FALSE)

########################## 2. read published results (nature 2017) ####################
credibleSets <- read.table('~/yuzhao1/work/final_GCAatac/12snps/1nature2017cs/nature2017cs.csv', sep = ',', header = T)
snps_pip <- read.table('~/yuzhao1/work/final_GCAatac/12snps/1nature2017cs/snps.csv', sep = ',', header = T)
snps_pip <- snps_pip[!is.na(snps_pip$chr),] # remove NA rows based on chr

# there is an error in the snps_pip (rs3838334 should be in TNFSF8, remove the wrong one)
snps_pip <- snps_pip[!(snps_pip$variant == 'rs3838334' & snps_pip$Gene.refGene != 'TNFSF8'),]
rownames(snps_pip) <- snps_pip$variant

credibleSets <- strsplit(credibleSets$all.variant, split = ',')
# credibleSets <- credibleSets[sapply(credibleSets, length) <= 50] # choose whether to filter out large cs
nt2017_rsIDs <- unique(unlist(credibleSets))
setdiff(nt2017_rsIDs, snps_pip$variant)
setdiff(snps_pip$variant, nt2017_rsIDs)

# confirmed this is the final snp df
nt2017_snps <- snps_pip

# create a signal ID (154 signals in 94 loci)
nt2017_snps$signal_ID <- paste0(nt2017_snps$HD, '_', nt2017_snps$signal) 

# compare the leading snps with our results
leadingSNPSs_nt2017 <- nt2017_snps[nt2017_snps$P_mean_95 > 0.5, 'variant']
intersect(leadingSNPSs_nt2017, leadingSNPSs_susieL10)

# compare all snps with our results and show their pip distribution
shared_snps <- intersect(nt2017_snps$variant, susie_clean_L10$snp)
hist(susie_clean_L10[shared_snps, 'susie_pip'], breaks = 50)
hist(nt2017_snps[shared_snps, 'P_mean_95'], breaks = 50)


############################ 3. read published results (NG 2023) ####################

snp320_hg38 <- readRDS('~/yuzhao1/work/final_GCAatac/12snps/data/00summarized_hg38.rds')

intersect(susie_clean_L5$snp, snp320_hg38$rsID)
intersect(nt2017_rsIDs, snp320_hg38$rsID)
intersect(nt2017_rsIDs, snp320_hg38$rsID) %>% intersect(., susie_clean_L5$snp)


########################## 4. check the overlapping SNPs ####################
# for each signal in susie results, check whether there is a signal in published results that shares at least one SNP
grouped_list1 <- split(susie_clean_L10$snp, susie_clean_L10$signal_ID)
grouped_list2 <- split(nt2017_snps$variant, nt2017_snps$signal_ID)
matched <- c()
for (i in 1:length(grouped_list1)){
  # check whether there is a matched signal in list2 (overlapped at least 1 snp)
  if(length(grouped_list1[[i]]) > 1){
    next
  }
  print(grouped_list1[[i]])
  flag <- 0
  snps_in_signal <- grouped_list1[[i]]
  if(length(intersect(snps_in_signal, unlist(grouped_list2))) > 0){
    flag <- 1
  }
  matched <- c(matched, flag)
}
table(matched)


########################## 5, check with known in deLange ####################
## 5.1 check how many SNPs are within known loci (241 loci from de Lange)
signals_known <- read.table('~/yuzhao1/work/final_GCAatac/12snps/3natureGenetics2017/ng2017_knownSNPs.csv', header = T, skip = 8, sep = ',')
signals_known <- signals_known[!is.na(signals_known$Chr), ]
signals_known_bed <- data.frame(chr = paste0('chr', signals_known$Chr),
                                            start = signals_known$LD_left,
                                            end = signals_known$LD_right)
write.table(signals_known_bed, '~/yuzhao1/work/final_GCAatac/12snps/3natureGenetics2017/signals_known.bed', 
            row.names = F, col.names = F, sep="\t", quote=FALSE)


# bedtools intersect -wa -a ~/yuzhao1/work/final_GCAatac/8susie/results/finemapped_snps_L10.bed -b ~/yuzhao1/work/final_GCAatac/12snps/3natureGenetics2017/signals_known.bed > ~/yuzhao1/work/final_GCAatac/8susie/results/finemapped_snps_L10_checkKnownWithDelange.bed


## 5.2 check how many loci are overlapping with known loci (241 loci from de Lange)
loci_identified <- LD_Blocks[LD_Blocks$X4 %in% unique(susie_clean_L10$locus),]
loci_identified_bed <- data.frame(chr = paste0('chr', loci_identified$X1),
                                start = loci_identified$X2,
                                end = loci_identified$X3)
write.table(loci_identified_bed, '~/yuzhao1/work/final_GCAatac/8susie/results/loci_identified_L10.bed', 
            row.names = F, col.names = F, sep="\t", quote=FALSE)

# bedtools intersect -u -a ~/yuzhao1/work/final_GCAatac/8susie/results/loci_identified_L10.bed -b ~/yuzhao1/work/final_GCAatac/12snps/3natureGenetics2017/signals_known.bed > ~/yuzhao1/work/final_GCAatac/8susie/results/finemapped_loci_L10_checkKnownWithDelange.bed



######################## 6, check with known in Liu ####################
## 6.1 check how many SNPs are within known loci (320 loci)
signals_known <- read.table('~/yuzhao1/work/final_GCAatac/12snps/4natureGenetics2023/IBDloci_EASandEUR.csv', sep = ',')
colnames(signals_known) <- signals_known[2,]
signals_known <- signals_known[-c(1,2),]
signals_known <- signals_known[signals_known$Chr_index!='', ]
signals_known_bed <- data.frame(chr = paste0('chr', signals_known$Chr),
                                start = signals_known$BP_left_loci,
                                end = signals_known$BP_right_loci)
write.table(signals_known_bed, '~/yuzhao1/work/final_GCAatac/12snps/4natureGenetics2023/IBDloci_EASandEUR.bed', 
            row.names = F, col.names = F, sep="\t", quote=FALSE)

# bedtools intersect -wa -a ~/yuzhao1/work/final_GCAatac/8susie/results/finemapped_snps_L10.bed -b ~/yuzhao1/work/final_GCAatac/12snps/4natureGenetics2023/IBDloci_EASandEUR.bed > ~/yuzhao1/work/final_GCAatac/8susie/results/finemapped_snps_L10_checkKnownWithLiu.bed

# bedtools intersect -u -a ~/yuzhao1/work/final_GCAatac/8susie/results/loci_identified_L10.bed -b ~/yuzhao1/work/final_GCAatac/12snps/4natureGenetics2023/IBDloci_EASandEUR.bed > ~/yuzhao1/work/final_GCAatac/8susie/results/finemapped_loci_L10_checkKnownWithLiu.bed


 


