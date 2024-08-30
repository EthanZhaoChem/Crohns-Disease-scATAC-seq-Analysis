library(bigsnpr)
library(mapgen)
library(susieR)
library(data.table)
library(dplyr)
library(stringr)

## 1. read necessary info 
out.dir <- '~/yuzhao1/work/final_GCAatac/8susie/results/'
diagnostic_folder <- '~/yuzhao1/work/final_GCAatac/8susie/diagnostic_plots_ukbb_filtered/'
n <- 40266
gwas_folder <- '~/yuzhao1/work/final_GCAatac/0gwas/b151_hg19/'
gwas_raw <- as.data.frame(
  fread(
    paste0(gwas_folder, "cd_build37_40266_20161107_hg19_b151_gwas_annotated_notFilterFor1KG.txt")
  )
)
LD_Blocks <- readRDS(system.file('extdata', 'LD.blocks.EUR.hg19.rds', package='mapgen'))
region_info <- get_UKBB_region_info(LD_Blocks,
                                    LDREF.dir = "/project2/mstephens/wcrouse/UKB_LDR_0.1_b37", 
                                    prefix = "ukb_b37_0.1")
LD_snp_info <- read_LD_SNP_info(region_info)

## 2. process gwas
gwas.sumstats <- process_gwas_sumstats(gwas_raw, 
                                       chr='chr', pos='pos', beta='beta', se='se', 
                                       a0='a0', a1='a1', snp='snp', pval='pval',
                                       LD_snp_info=LD_snp_info, 
                                       strand_flip=TRUE, 
                                       remove_strand_ambig=TRUE)
sig.loci <- gwas.sumstats %>% dplyr::filter(pval < 5e-8) %>% dplyr::pull(locus) %>% unique()
cat(length(sig.loci), "significant loci. \n")

## 3. diagnostic 
mismatch_rsID_list <- list()
for (locus in sig.loci) {
  cat('Calculating: \n')
  cat(locus, '\n')
  sumstats.locus <- gwas.sumstats[gwas.sumstats$locus == locus, ]
  # load LD matrix and SNP info for this locus
  LD_ref <- load_UKBB_LDREF(LD_Blocks, 
                            locus, 
                            LDREF.dir = "/project2/mstephens/wcrouse/UKB_LDR_0.1_b37", 
                            prefix = "ukb_b37_0.1")
  
  # Match GWAS sumstats with LD reference, only keep variants included in LD reference.
  matched.sumstat.LD <- match_gwas_LDREF(sumstats.locus, LD_ref$R, LD_ref$snp_info)
  sumstats.locus <- matched.sumstat.LD$sumstats
  z.locus <- sumstats.locus$zscore
  R.locus <- matched.sumstat.LD$R
  
  # diagnostic
  condz <- kriging_rss(z.locus, R.locus, n=n)
  condz$conditional_dist$rsID <- sumstats.locus$snp
  xx <- condz$conditional_dist
  
  # save mismatch snps
  locus_ID <- as.character(locus)
  mismatch_rsID_list[[locus_ID]] <- xx[which(abs(xx$z) > 2 & xx$logLR > 2), 'rsID']
  
  # save diagnostic plots
  pdf(paste0(diagnostic_folder, '', locus, '.pdf'), width = 5, height = 5)
  print(condz$plot)
  dev.off()
}

saveRDS(mismatch_rsID_list, paste0(out.dir, 'diagnostic_mismatch_rsID_list.rds'))








