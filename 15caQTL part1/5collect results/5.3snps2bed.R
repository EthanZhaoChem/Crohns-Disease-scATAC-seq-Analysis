dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(dplyr)
library(tidyverse)
library(data.table)
library(readxl)    

## presets
out.dir <- '~/yuzhao1/work/atac_gca2024/19rasqual/5results/snps/'

## read snps
calculated.rootDir <- '~/yuzhao1/work/atac_gca2024/19rasqual/5results/calculated/'
celltypes <- readLines('~/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt')

filename <- paste0(calculated.rootDir, 'Perm_five_times_FDR_0.05.xlsx')
sheets <- readxl::excel_sheets(filename)
caQTL.list <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
caQTL.list <- lapply(caQTL.list, as.data.frame)
names(caQTL.list) <- sheets

snps.list <- list()
for (ct in celltypes) {
  ct.caqtl <- caQTL.list[[ct]]
  snps.list[[ct]] <- ct.caqtl$rsID
}

snps.all <- unique(unlist(snps.list))

## prepare df snps
snps_hg38_bed_df <- data.frame(chr=snps.all %>% strsplit(split = ':', fixed=T) %>% sapply(.,`[[`,1),
                           start=snps.all %>% strsplit(split = ':', fixed=T) %>% sapply(.,`[[`,2) %>% as.numeric() - 1,
                           end=snps.all %>% strsplit(split = ':', fixed=T) %>% sapply(.,`[[`,2) %>% as.numeric())


write.table(snps_hg38_bed_df, file=paste0(out.dir, 'hg38/all_30385.bed'), quote=F, sep="\t", row.names=F, col.names=F)


## prepare nature 2017 finemapping results
snps_pip_nt2017 <- read.table('~/yuzhao1/work/final_GCAatac/12snps/1nature2017cs/snps.csv', sep = ',', header = T)
snps_pip_nt2017 <- snps_pip_nt2017[!is.na(snps_pip_nt2017$chr),] # remove NA rows based on chr

xx <- data.frame(chr=paste0('chr', snps_pip_nt2017$chr),
                 start=snps_pip_nt2017$position-1,
                 end=snps_pip_nt2017$position)


write.table(xx, file=paste0(out.dir, 'hg19/bed_lifted/snps_nt2017_4312.bed'), quote=F, sep="\t", row.names=F, col.names=F)










