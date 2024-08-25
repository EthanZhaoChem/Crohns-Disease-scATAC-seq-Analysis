dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(rasqualTools)
library(ArchR)
library(DESeq2)
celltypes <- readLines('~/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt')
raw.dir <- '~/yuzhao1/work/atac_gca2024/19rasqual/2.1mtx_peak_patient_raw/'
out.dir <- '~/yuzhao1/work/atac_gca2024/19rasqual/9boxplot/9.1mtx/'

for(ct in celltypes){
  mtx.raw <- read.table(paste0(raw.dir, ct, '.txt')) # we calculated all peaks
  mtx <- varianceStabilizingTransformation(as.matrix(mtx.raw), blind = F)
  write.table(mtx, paste0(out.dir, ct, ".vst.count_matrix"), sep="\t", quote=F, row.names=T)
}











