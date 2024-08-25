dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(rasqualTools)
library(ArchR)
library(DESeq2)
celltypes <- readLines('~/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt')
raw.dir <- '~/yuzhao1/work/atac_gca2024/19rasqual/9boxplot/9.2genotype_raw/'
out.dir <- '~/yuzhao1/work/atac_gca2024/19rasqual/9boxplot/9.3genotype_clean/'

for(ct in celltypes){
  samples <- readLines(paste0(raw.dir, ct, '.samples'))
  snps = read.table(paste0(raw.dir, ct, '.genotype'), header=F, stringsAsFactors = F)
  snps <- replace(snps, snps == "/", "|")
  snps <- replace(snps, snps == "1|1", 2)
  snps <- replace(snps, snps == "0|1", 1)
  snps <- replace(snps, snps == "1|0", 1)
  snps <- replace(snps, snps == "0|0", 0)
  colnames(snps) = c("CHR", "POS", "ID", 'REF', 'ALT', samples) 
  genotypes <- as.matrix(snps[,samples])
  rownames(genotypes) <- snps$ID
  write.table(genotypes, paste0(out.dir, ct, '_genotypes.tsv'), sep="\t", quote=F)
}











