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

# only consider SNV class
gwas <- gwas[gwas$VC == 'SNV',]

# only consider those with 1000G MAF
gwas <- gwas[!is.na(gwas$MAF),]

# further remove indels, these are not SNV class in 1kg, but indels in gwas
gwas$Allele1 <- toupper(gwas$Allele1)
gwas$Allele2 <- toupper(gwas$Allele2)
nucs <- c("A", "C", "T", "G")
bola1 <- (gwas$Allele1 %in% nucs)
bola2 <- (gwas$Allele2 %in% nucs)
gwas <- gwas[(bola1 & bola2),]

# extract AF_ref AF_alt
gwas_info_elements <- strsplit(gwas$INFO, split = ';', fixed = T) %>% unlist()
caf_vector <- grep("^CAF=", gwas_info_elements, value = TRUE)
caf_vector_numbers <- str_extract_all(caf_vector, "\\b[0-9.]+\\b")
vector_lengths <- sapply(caf_vector_numbers, length)

# keep those 1000Genomes alternate allele is in the dbSNPs alternate allele set (2 numbers in CAF)
gwas <- gwas[which(vector_lengths == 2),]
caf_vector_numbers <- caf_vector_numbers[which(vector_lengths == 2)]
gwas$AF_ref <- caf_vector_numbers %>% sapply(.,`[[`,1) # ref allele
gwas$AF_alt <- caf_vector_numbers %>% sapply(.,`[[`,2) # alt allele
gwas$AF_ref <- as.numeric(gwas$AF_ref)
gwas$AF_alt <- as.numeric(gwas$AF_alt)

# add chr and pos
gwas$chr <- gsub(pattern = "(.*):(.*)_(.*)_(.*)", gwas$MarkerName, replacement = "\\1") %>% as.numeric()
table(gwas$chr, useNA = "ifany")
gwas$pos <- gsub(pattern = "(.*):(.*)_(.*)_(.*)", gwas$MarkerName, replacement = "\\2") %>% as.numeric()

# use common column names
gwas_cleaned <- gwas[, c('rsID', "chr", "pos", "Effect", "StdErr", "P.value", "Allele2", "Allele1", 'MAF', 'AF_alt', 'AF_ref')]
colnames(gwas_cleaned) <- c('snp', "chr", "pos", "beta", "se", "pval", "a0", "a1", 'MAF', 'AF_a1', 'AF_a0')

# final format correction and sanity check
gwas_cleaned$MAF <- as.numeric(gwas_cleaned$MAF)
sum(abs(pmin(gwas_cleaned$AF_a1, gwas_cleaned$AF_a0) - gwas_cleaned$MAF) == 0) == nrow(gwas_cleaned)

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
       paste0(gwas_folder, 'cd_build37_40266_20161107_hg19_b151_gwas_annotated_filtered_SNVs.txt'),
       quote=FALSE, sep='\t', row.names=FALSE)



##################### B, uc #####################
gwas<- as.data.frame(
  fread(
    paste0(gwas_folder, "uc_build37_45975_20161107_annotated_hg19_b151_gwas_annotated.txt")
  )
)

# only consider SNV class
gwas <- gwas[gwas$VC == 'SNV',]

# only consider those with 1000G MAF
gwas <- gwas[!is.na(gwas$MAF),]

# further remove indels, these are not SNV class in 1kg, but indels in gwas
gwas$Allele1 <- toupper(gwas$Allele1)
gwas$Allele2 <- toupper(gwas$Allele2)
nucs <- c("A", "C", "T", "G")
bola1 <- (gwas$Allele1 %in% nucs)
bola2 <- (gwas$Allele2 %in% nucs)
gwas <- gwas[(bola1 & bola2),]

# extract AF_ref AF_alt
gwas_info_elements <- strsplit(gwas$INFO, split = ';', fixed = T) %>% unlist()
caf_vector <- grep("^CAF=", gwas_info_elements, value = TRUE)
caf_vector_numbers <- str_extract_all(caf_vector, "\\b[0-9.]+\\b")
vector_lengths <- sapply(caf_vector_numbers, length)

# keep those 1000Genomes alternate allele is in the dbSNPs alternate allele set (2 numbers in CAF)
gwas <- gwas[which(vector_lengths == 2),]
caf_vector_numbers <- caf_vector_numbers[which(vector_lengths == 2)]
gwas$AF_ref <- caf_vector_numbers %>% sapply(.,`[[`,1) # ref allele
gwas$AF_alt <- caf_vector_numbers %>% sapply(.,`[[`,2) # alt allele
gwas$AF_ref <- as.numeric(gwas$AF_ref)
gwas$AF_alt <- as.numeric(gwas$AF_alt)

# add chr and pos
gwas$chr <- gsub(pattern = "(.*):(.*)_(.*)_(.*)", gwas$MarkerName, replacement = "\\1") %>% as.numeric()
table(gwas$chr, useNA = "ifany")
gwas$pos <- gsub(pattern = "(.*):(.*)_(.*)_(.*)", gwas$MarkerName, replacement = "\\2") %>% as.numeric()

# use common column names
gwas_cleaned <- gwas[, c('rsID', "chr", "pos", "Effect", "StdErr", "P.value", "Allele2", "Allele1", 'MAF', 'AF_alt', 'AF_ref')]
colnames(gwas_cleaned) <- c('snp', "chr", "pos", "beta", "se", "pval", "a0", "a1", 'MAF', 'AF_a1', 'AF_a0')

# final format correction and sanity check
gwas_cleaned$MAF <- as.numeric(gwas_cleaned$MAF)
sum(abs(pmin(gwas_cleaned$AF_a1, gwas_cleaned$AF_a0) - gwas_cleaned$MAF) == 0) == nrow(gwas_cleaned)

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
       paste0(gwas_folder, 'uc_build37_45975_20161107_hg19_b151_gwas_annotated_filtered_SNVs.txt'),
       quote=FALSE, sep='\t', row.names=FALSE)


####################### C, ibd #####################
gwas<- as.data.frame(
  fread(
    paste0(gwas_folder, "ibd_build37_59957_20161107_annotated_hg19_b151_gwas_annotated.txt")
  )
)

# only consider SNV class
gwas <- gwas[gwas$VC == 'SNV',]

# only consider those with 1000G MAF
gwas <- gwas[!is.na(gwas$MAF),]

# further remove indels, these are not SNV class in 1kg, but indels in gwas
gwas$Allele1 <- toupper(gwas$Allele1)
gwas$Allele2 <- toupper(gwas$Allele2)
nucs <- c("A", "C", "T", "G")
bola1 <- (gwas$Allele1 %in% nucs)
bola2 <- (gwas$Allele2 %in% nucs)
gwas <- gwas[(bola1 & bola2),]

# extract AF_ref AF_alt
gwas_info_elements <- strsplit(gwas$INFO, split = ';', fixed = T) %>% unlist()
caf_vector <- grep("^CAF=", gwas_info_elements, value = TRUE)
caf_vector_numbers <- str_extract_all(caf_vector, "\\b[0-9.]+\\b")
vector_lengths <- sapply(caf_vector_numbers, length)

# keep those 1000Genomes alternate allele is in the dbSNPs alternate allele set (2 numbers in CAF)
gwas <- gwas[which(vector_lengths == 2),]
caf_vector_numbers <- caf_vector_numbers[which(vector_lengths == 2)]
gwas$AF_ref <- caf_vector_numbers %>% sapply(.,`[[`,1) # ref allele
gwas$AF_alt <- caf_vector_numbers %>% sapply(.,`[[`,2) # alt allele
gwas$AF_ref <- as.numeric(gwas$AF_ref)
gwas$AF_alt <- as.numeric(gwas$AF_alt)

# add chr and pos
gwas$chr <- gsub(pattern = "(.*):(.*)_(.*)_(.*)", gwas$MarkerName, replacement = "\\1") %>% as.numeric()
table(gwas$chr, useNA = "ifany")
gwas$pos <- gsub(pattern = "(.*):(.*)_(.*)_(.*)", gwas$MarkerName, replacement = "\\2") %>% as.numeric()

# use common column names
gwas_cleaned <- gwas[, c('rsID', "chr", "pos", "Effect", "StdErr", "P.value", "Allele2", "Allele1", 'MAF', 'AF_alt', 'AF_ref')]
colnames(gwas_cleaned) <- c('snp', "chr", "pos", "beta", "se", "pval", "a0", "a1", 'MAF', 'AF_a1', 'AF_a0')

# final format correction and sanity check
gwas_cleaned$MAF <- as.numeric(gwas_cleaned$MAF)
sum(abs(pmin(gwas_cleaned$AF_a1, gwas_cleaned$AF_a0) - gwas_cleaned$MAF) == 0) == nrow(gwas_cleaned)

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
       paste0(gwas_folder, 'ibd_build37_59957_20161107_hg19_b151_gwas_annotated_filtered_SNVs.txt'),
       quote=FALSE, sep='\t', row.names=FALSE)

