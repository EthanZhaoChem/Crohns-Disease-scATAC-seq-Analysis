dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)

source('~/yuzhao1/scripts/plot.R')
addArchRThreads(24)

# # read files
proj_union <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1/")


temp_lineage <- 'union'
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)
proj <- addPeakMatrix(proj, force = T, ceiling=10^9)
saveArchRProject(proj)


peak.mtx <- getMatrixFromProject(proj, useMatrix = "PeakMatrix", binarize = F)
saveRDS(peak.mtx, paste0(proj@projectMetadata$outputDirectory, '/peakMat_unbinarized.rds'))









