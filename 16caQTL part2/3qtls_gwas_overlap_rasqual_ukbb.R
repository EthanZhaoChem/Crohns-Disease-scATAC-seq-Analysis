library(stringr)
library(dplyr)
library(plyr)
library(tidyverse)
library(data.table)
library(limma)

## called snps
FDR_results_list <- readRDS('~/yuzhao1/work/atac_gca2024/24rasqual2/1rasqual_results/FDR_results_list.rds')
qtls_lead_snps_list <- lapply(FDR_results_list, function(x) paste0(x$Chromosome, '_', x$SNP_position, '_', x$SNP_position))
qtls_lead_snps <- unique(unlist(qtls_lead_snps_list))

# read gwas lead snps
gwas_lead_df <- read.csv('~/yuzhao1/work/atac_gca2024/19rasqual/8LDoverlap/gwas_cd/cd_finemapping_unifprior_93loci_L10_annotated_ukbb.csv')
gwas_lead_snps <- gwas_lead_df$name_chr_pos_hg38

# read ld corr 
ld_corr <- read.table('~/yuzhao1/work/atac_gca2024/19rasqual/8LDoverlap/ld_1kg_and_ourSample/ld_0.8/all_sample_1KGP_ld.ld', header = T)
ld_corr$name_chr_pos_hg38_A <- paste0('chr', ld_corr$CHR_A, '_', ld_corr$BP_A, '_', ld_corr$BP_A)
ld_corr$name_chr_pos_hg38_B <- paste0('chr', ld_corr$CHR_B, '_', ld_corr$BP_B, '_', ld_corr$BP_B)


###################################################################
# overlapped gwas with qtls: exactly same, or in LD with each other
snps_shared <- intersect(qtls_lead_snps, gwas_lead_snps)
ld_corr_sub1 <- ld_corr[ld_corr$name_chr_pos_hg38_A %in% gwas_lead_snps & ld_corr$name_chr_pos_hg38_B %in% qtls_lead_snps,]
ld_corr_sub2 <- ld_corr[ld_corr$name_chr_pos_hg38_A %in% qtls_lead_snps & ld_corr$name_chr_pos_hg38_B %in% gwas_lead_snps,]
overlapped_snp_pairs <- data.frame(gwas_lead = c(snps_shared, ld_corr_sub1$name_chr_pos_hg38_A, ld_corr_sub2$name_chr_pos_hg38_B),
                                   qtls_lead = c(snps_shared, ld_corr_sub1$name_chr_pos_hg38_B, ld_corr_sub2$name_chr_pos_hg38_A),
                                   R2 = c(rep(1, length(snps_shared)), ld_corr_sub1$R2, ld_corr_sub2$R2))

# add cell type specificity
overlapped_snp_pairs$celltypes <- ''
overlapped_snp_pairs$n_celltypes <- 0
for (ct in names(qtls_lead_snps_list)) {
  ct_qtls <- qtls_lead_snps_list[[ct]]
  for (i in 1:nrow(overlapped_snp_pairs)) {
    tmp.qtl <- overlapped_snp_pairs[i, 'qtls_lead']
    if(tmp.qtl %in% ct_qtls){
      overlapped_snp_pairs[i, 'n_celltypes'] <- overlapped_snp_pairs[i, 'n_celltypes'] + 1
      if(overlapped_snp_pairs[i, 'celltypes'] != ''){
        overlapped_snp_pairs[i, 'celltypes'] <- paste0(overlapped_snp_pairs[i, 'celltypes'], '-', ct)
      }else{
        overlapped_snp_pairs[i, 'celltypes'] <- ct
      }
    }
  }
}

# mash QTLs
write.csv(overlapped_snp_pairs, '/project/gca/yuzhao1/work/atac_gca2024/24rasqual2/3qtls_gwas_overlap/gwas_rasqual_caQTL_LD08_ukbbFinemapping.csv')


