# calculate LD information for the snps we are interested in:
library(dplyr)
library(plyr)
results <- readRDS('~/yuzhao1/work/atac_gca2024/0manu/plots/3snps_susie_L10_prioritization/results.rds')
snps <- names(results$snps_list_full)

gwas <- read.csv('~/yuzhao1/work/atac_gca2024/0manu/plots/3snps_susie_L10/gwas_finemapped_processed.csv', row.names = 1)
gwas <- gwas[snps,]


xx <- read.table('~/yuzhao1/work/atac_gca2024/19rasqual/8LDoverlap/ld_1kg_and_ourSample/all_sample_1KGP.bim')
yy <- paste0(strsplit(xx$V2, split = '_', fixed=T) %>% sapply(.,`[[`,1), '_',
             strsplit(xx$V2, split = '_', fixed=T) %>% sapply(.,`[[`,2), '_',
             strsplit(xx$V2, split = '_', fixed=T) %>% sapply(.,`[[`,2))
xx[which(yy %in% gwas$name_chr_pos_hg38), 'V2']


sum(rownames(gwas) %in% bigSNP$map$marker.ID)
