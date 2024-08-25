########################### prepare finemapped snps #######################
library(plyr)
gwas_annotated_hg19 <- read.table('~/yuzhao1/work/final_GCAatac/0gwas/b151_hg19/cd_build37_40266_20161107_hg19_b151_gwas_annotated_notFilterFor1KG.txt', header = T)
gwas_annotated_hg19$name_chr_pos_hg19 <- paste0('chr', gwas_annotated_hg19$chr, '_', gwas_annotated_hg19$pos, '_', gwas_annotated_hg19$pos)

gwas_annotated_hg38 <- read.table('~/yuzhao1/work/final_GCAatac/0gwas/b151_hg38/cd_40266_20161107_hg38_b151_gwas_annotated_notFilterFor1KG.txt', header = F)
colnames(gwas_annotated_hg38)[1:3] <- c('snp', 'chr', 'pos')
gwas_annotated_hg38$name_chr_pos_hg38 <- paste0(gwas_annotated_hg38$chr, '_', gwas_annotated_hg38$pos, '_', gwas_annotated_hg38$pos)


## 1. UKBB finemapped results
gwas_lead <- read.csv('~/yuzhao1/work/final_GCAatac/8susie/results/cd_finemapping_unifprior_L10.csv', row.names = 1)
gwas_annotated_hg19 <- gwas_annotated_hg19[gwas_annotated_hg19$snp %in% gwas_lead$snp, ]
gwas_annotated_hg38 <- gwas_annotated_hg38[gwas_annotated_hg38$snp %in% gwas_lead$snp, ]
gwas_lead$name_chr_pos_hg19 <- mapvalues(gwas_lead$snp, from = gwas_annotated_hg19$snp, to = gwas_annotated_hg19$name_chr_pos_hg19)
gwas_lead$name_chr_pos_hg38 <- mapvalues(gwas_lead$snp, from = gwas_annotated_hg38$snp, to = gwas_annotated_hg38$name_chr_pos_hg38)
write.csv(gwas_lead, '~/yuzhao1/work/atac_gca2024/19rasqual/8LDoverlap/gwas_cd/cd_finemapping_unifprior_93loci_L10_annotated_ukbb.csv')

## 2. 1KG finemapped results
gwas_lead <- read.csv('~/yuzhao1/work/final_GCAatac/8susie/results/cd_finemapping_unifprior_95loci_L10.csv', row.names = 1)
gwas_annotated_hg19 <- gwas_annotated_hg19[gwas_annotated_hg19$snp %in% gwas_lead$snp, ]
gwas_annotated_hg38 <- gwas_annotated_hg38[gwas_annotated_hg38$snp %in% gwas_lead$snp, ]
gwas_lead$name_chr_pos_hg19 <- mapvalues(gwas_lead$snp, from = gwas_annotated_hg19$snp, to = gwas_annotated_hg19$name_chr_pos_hg19)
gwas_lead$name_chr_pos_hg38 <- mapvalues(gwas_lead$snp, from = gwas_annotated_hg38$snp, to = gwas_annotated_hg38$name_chr_pos_hg38)
write.csv(gwas_lead, '~/yuzhao1/work/atac_gca2024/19rasqual/8LDoverlap/gwas_cd/cd_finemapping_unifprior_95loci_L10_annotated.csv')


########################### save snps nominated by p, not important #######################
# use EUR LD ref because this GWAS was mainly studying EUR samples
library(gwasglue)
library(ieugwasr)

# read top snps (p < 5e-08)
dat_eur_r2_0.8 <- ieugwasr::tophits("ebi-a-GCST004132", pop='EUR', r2 = 0.8, force_server=T)
saveRDS(dat_eur_r2_0.8, '~/yuzhao1/work/atac_gca2024/19rasqual/8LDoverlap/gwas_cd/eur_r2_0.8.rds')

# save gwas regions (in each region, determine an independent snp within 1Mb region (see LD R2 threshold above))
pop <- 'EUR'
regionfile <- system.file("extdata", "ldetect", paste0(pop, ".bed"), package="gwasglue")
regions <- data.table::fread(regionfile, header=TRUE) %>%
  dplyr::mutate(
    chr=as.numeric(gsub("chr", "", chr)),
    start=as.numeric(start),
    stop=as.numeric(stop)
  ) %>% dplyr::as_tibble()
saveRDS(regions, '~/yuzhao1/work/atac_gca2024/19rasqual/8LDoverlap/gwas_cd/gwas_regions.rds')


# annotate
library(plyr)
gwas_annotated_hg19 <- read.table('~/yuzhao1/work/final_GCAatac/0gwas/b151_hg19/cd_build37_40266_20161107_hg19_b151_gwas_annotated_notFilterFor1KG.txt', header = T)
gwas_annotated_hg19$name_chr_pos_hg19 <- paste0('chr', gwas_annotated_hg19$chr, '_', gwas_annotated_hg19$pos, '_', gwas_annotated_hg19$pos)

gwas_annotated_hg38 <- read.table('~/yuzhao1/work/final_GCAatac/0gwas/b151_hg38/cd_40266_20161107_hg38_b151_gwas_annotated_notFilterFor1KG.txt', header = F)
colnames(gwas_annotated_hg38)[1:3] <- c('snp', 'chr', 'pos')
gwas_annotated_hg38$name_chr_pos_hg38 <- paste0(gwas_annotated_hg38$chr, '_', gwas_annotated_hg38$pos, '_', gwas_annotated_hg38$pos)

gwas_lead <- readRDS('~/yuzhao1/work/atac_gca2024/19rasqual/8LDoverlap/gwas_cd/eur_r2_0.8.rds')
gwas_annotated_hg19 <- gwas_annotated_hg19[gwas_annotated_hg19$snp %in% gwas_lead$rsid, ]
gwas_annotated_hg38 <- gwas_annotated_hg38[gwas_annotated_hg38$snp %in% gwas_lead$rsid, ]
gwas_lead$name_chr_pos_hg19 <- mapvalues(gwas_lead$rsid, from = gwas_annotated_hg19$snp, to = gwas_annotated_hg19$name_chr_pos_hg19)
gwas_lead$name_chr_pos_hg38 <- mapvalues(gwas_lead$rsid, from = gwas_annotated_hg38$snp, to = gwas_annotated_hg38$name_chr_pos_hg38)
write.csv(gwas_lead, '~/yuzhao1/work/atac_gca2024/19rasqual/8LDoverlap/gwas_cd/eur_r2_0.8_annotated.csv')












