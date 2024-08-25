# read vcf info, AF in reference panel
df <- read.table('~/yuzhao1/work/final_GCAatac/18glimpse/7.1snps_filter/info_extracted.txt', header = F)
colnames(df) <- c('SNP', 'RAF', 'INFO')

# read MAF from sample cohort
maf_cohort <- read.table('~/yuzhao1/work/final_GCAatac/18glimpse/6.2maf/plink.frq', header = T)
if(sum(df$SNP == maf_cohort$SNP) == nrow(df)){
  df$MAF_cohort <- maf_cohort$MAF
}

# calculate MAF in reference panel
df$MAF_ref <- pmin(df$RAF, 1-df$RAF)

# filter by two MAF cutoffs and one INFO cutoff
snps_filtered <- df$SNP[df$INFO > 0.9 & df$MAF_ref > 0.01 & df$MAF_cohort > 0.05]
snps_filtered <- unique(snps_filtered)


sink("~/yuzhao1/work/final_GCAatac/18glimpse/7.1snps_filter/snps_filtered.txt", append = F)
cat(snps_filtered, sep = '\n')
sink()

# read filtered results
snps_test <- read.table('~/yuzhao1/work/final_GCAatac/18glimpse/7.1snps_filter/test.txt', sep = '\t')
colnames(snps_test) <- c('SNP', 'CHROM', 'REF', 'ALT')
xx <- table(snps_test$V1)


df_noname <- snps_test[snps_test$SNP=='.',]






