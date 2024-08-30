library(data.table)
library(dplyr)
library(stringr)
gwas_folder <- '~/yuzhao1/work/final_GCAatac/0gwas/b151_hg19/'

####################### A, cd #####################
gwas<- as.data.frame(
  fread(
    paste0(gwas_folder, "cd_build37_40266_20161107_annotated_hg19_b151_gwas_annotated.txt")
  )
)

# further remove indels, these are not SNV class in 1kg, but indels in gwas
gwas$Allele1 <- toupper(gwas$Allele1)
gwas$Allele2 <- toupper(gwas$Allele2)
nucs <- c("A", "C", "T", "G")
bola1 <- (gwas$Allele1 %in% nucs)
bola2 <- (gwas$Allele2 %in% nucs)
gwas <- gwas[(bola1 & bola2),]

# add chr and pos
gwas$chr <- gsub(pattern = "(.*):(.*)_(.*)_(.*)", gwas$MarkerName, replacement = "\\1") %>% as.numeric()
gwas$pos <- gsub(pattern = "(.*):(.*)_(.*)_(.*)", gwas$MarkerName, replacement = "\\2") %>% as.numeric()

# use common column names
gwas_cleaned <- gwas[, c('rsID', "chr", "pos", "Effect", "StdErr", "P.value", "Allele2", "Allele1")]
colnames(gwas_cleaned) <- c('snp', "chr", "pos", "beta", "se", "pval", "a0", "a1")

# zscore
gwas_cleaned$zscore <- gwas_cleaned$beta / gwas_cleaned$se
nrow(gwas_cleaned)

# Remove SNPs without z-score, file unchanged
gwas_cleaned <- gwas_cleaned[!is.na(gwas_cleaned$zscore),]
nrow(gwas_cleaned)  

# Remove multi-allelic SNPs, file unchanged
chrpos <- paste0(gwas_cleaned$chr, "_", gwas_cleaned$pos)
gwas_cleaned <- gwas_cleaned[!duplicated(chrpos),]
nrow(gwas_cleaned)  

# add ld block by mapgen package, optional step
library(mapgen)
data('Euro_LD_Chunks', package='mapgen')
gwas_cleaned <- assign_locus_snp(cleaned.sumstats = gwas_cleaned, ld = LD_Blocks)
nrow(gwas_cleaned)  

# add sample size, assuming it is from the file name
gwas_cleaned$N <- 40266

fwrite(gwas_cleaned, 
       paste0(gwas_folder, 'cd_build37_40266_20161107_hg19_b151_gwas_annotated_notFilterFor1KG.txt'),
       quote=FALSE, sep='\t', row.names=FALSE)




##################### B, uc #####################
gwas<- as.data.frame(
  fread(
    paste0(gwas_folder, "uc_build37_45975_20161107_annotated_hg19_b151_gwas_annotated.txt")
  )
)

# further remove indels, these are not SNV class in 1kg, but indels in gwas
gwas$Allele1 <- toupper(gwas$Allele1)
gwas$Allele2 <- toupper(gwas$Allele2)
nucs <- c("A", "C", "T", "G")
bola1 <- (gwas$Allele1 %in% nucs)
bola2 <- (gwas$Allele2 %in% nucs)
gwas <- gwas[(bola1 & bola2),]

# add chr and pos
gwas$chr <- gsub(pattern = "(.*):(.*)_(.*)_(.*)", gwas$MarkerName, replacement = "\\1") %>% as.numeric()
gwas$pos <- gsub(pattern = "(.*):(.*)_(.*)_(.*)", gwas$MarkerName, replacement = "\\2") %>% as.numeric()

# use common column names
gwas_cleaned <- gwas[, c('rsID', "chr", "pos", "Effect", "StdErr", "P.value", "Allele2", "Allele1")]
colnames(gwas_cleaned) <- c('snp', "chr", "pos", "beta", "se", "pval", "a0", "a1")

# zscore
gwas_cleaned$zscore <- gwas_cleaned$beta / gwas_cleaned$se
nrow(gwas_cleaned)

# Remove SNPs without z-score, file unchanged
gwas_cleaned <- gwas_cleaned[!is.na(gwas_cleaned$zscore),]
nrow(gwas_cleaned)  

# Remove multi-allelic SNPs, file unchanged
chrpos <- paste0(gwas_cleaned$chr, "_", gwas_cleaned$pos)
gwas_cleaned <- gwas_cleaned[!duplicated(chrpos),]
nrow(gwas_cleaned)  

# add ld block by mapgen package, optional step
library(mapgen)
data('Euro_LD_Chunks', package='mapgen')
gwas_cleaned <- assign_locus_snp(cleaned.sumstats = gwas_cleaned, ld = LD_Blocks)
nrow(gwas_cleaned)  

# add sample size, assuming it is from the file name
gwas_cleaned$N <- 45975

fwrite(gwas_cleaned, 
       paste0(gwas_folder, 'uc_build37_45975_20161107_hg19_b151_gwas_annotated_notFilterFor1KG.txt'),
       quote=FALSE, sep='\t', row.names=FALSE)




####################### C, ibd #####################
gwas<- as.data.frame(
  fread(
    paste0(gwas_folder, "ibd_build37_59957_20161107_annotated_hg19_b151_gwas_annotated.txt")
  )
)

# further remove indels, these are not SNV class in 1kg, but indels in gwas
gwas$Allele1 <- toupper(gwas$Allele1)
gwas$Allele2 <- toupper(gwas$Allele2)
nucs <- c("A", "C", "T", "G")
bola1 <- (gwas$Allele1 %in% nucs)
bola2 <- (gwas$Allele2 %in% nucs)
gwas <- gwas[(bola1 & bola2),]

# add chr and pos
gwas$chr <- gsub(pattern = "(.*):(.*)_(.*)_(.*)", gwas$MarkerName, replacement = "\\1") %>% as.numeric()
gwas$pos <- gsub(pattern = "(.*):(.*)_(.*)_(.*)", gwas$MarkerName, replacement = "\\2") %>% as.numeric()

# use common column names
gwas_cleaned <- gwas[, c('rsID', "chr", "pos", "Effect", "StdErr", "P.value", "Allele2", "Allele1")]
colnames(gwas_cleaned) <- c('snp', "chr", "pos", "beta", "se", "pval", "a0", "a1")

# zscore
gwas_cleaned$zscore <- gwas_cleaned$beta / gwas_cleaned$se
nrow(gwas_cleaned)

# Remove SNPs without z-score, file unchanged
gwas_cleaned <- gwas_cleaned[!is.na(gwas_cleaned$zscore),]
nrow(gwas_cleaned)  

# Remove multi-allelic SNPs, file unchanged
chrpos <- paste0(gwas_cleaned$chr, "_", gwas_cleaned$pos)
gwas_cleaned <- gwas_cleaned[!duplicated(chrpos),]
nrow(gwas_cleaned)  

# add ld block by mapgen package, optional step
library(mapgen)
data('Euro_LD_Chunks', package='mapgen')
gwas_cleaned <- assign_locus_snp(cleaned.sumstats = gwas_cleaned, ld = LD_Blocks)
nrow(gwas_cleaned)  

# add sample size, assuming it is from the file name
gwas_cleaned$N <- 59957

fwrite(gwas_cleaned, 
       paste0(gwas_folder, 'ibd_build37_59957_20161107_hg19_b151_gwas_annotated_notFilterFor1KG.txt'),
       quote=FALSE, sep='\t', row.names=FALSE)


