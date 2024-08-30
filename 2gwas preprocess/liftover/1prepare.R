dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(dplyr)
library(tidyverse)
library(data.table)
library(magrittr)
"%&%" <- function(a, b) paste0(a, b)

# preset
dir1 <- '~/yuzhao1/work/final_GCAatac/0gwas/b151_hg19/'
filename1 <- 'cd_build37_40266_20161107_hg19_b151_gwas_annotated_notFilterFor1KG.txt'
dir2 <- '~/yuzhao1/work/final_GCAatac/0gwas/b151_hg38/'
filename2 <- 'cd_40266_20161107_hg38_b151_gwas_annotated_notFilterFor1KG.txt'
dir_tmp <- '~/yuzhao1/work/final_GCAatac/0gwas/tmp/'

# write bed file
gwas <- fread(dir1 %&% filename1)
gwas$chr <- as.integer(gwas$chr)
gwas$pos <- as.integer(gwas$pos)
df <- data.frame(chr = paste0('chr', gwas$chr),
                 start = gwas$pos - 1,
                 end = gwas$pos,
                 rsID = gwas$snp)
df <- format(df, scientific = FALSE)
write.table(df, dir_tmp%&%'hg19.bed', 
            row.names = F,
            col.names = F, 
            sep="\t", 
            quote=FALSE)

