dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)

source('~/yuzhao1/scripts/plot.R')
addArchRThreads(1)

# # read files
proj_union <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_subSampling_max100/")
proj_epithelial <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_subSampling_max100_topic_epithelial/")
proj_immune <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_subSampling_max100_topic_immune/")
proj_stromal <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_subSampling_max100_topic_stromal/")


temp_lineage <- 'union'
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)
peak.mtx <- getMatrixFromProject(proj, useMatrix = "PeakMatrix", binarize = T)
saveRDS(peak.mtx, paste0(proj@projectMetadata$outputDirectory, '/peakMat_binarized.rds'))

temp_lineage <- 'epithelial'
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)
peak.mtx <- getMatrixFromProject(proj, useMatrix = "PeakMatrix", binarize = T)
saveRDS(peak.mtx, paste0(proj@projectMetadata$outputDirectory, '/peakMat_binarized.rds'))

temp_lineage <- 'immune'
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)
peak.mtx <- getMatrixFromProject(proj, useMatrix = "PeakMatrix", binarize = T)
saveRDS(peak.mtx, paste0(proj@projectMetadata$outputDirectory, '/peakMat_binarized.rds'))

temp_lineage <- 'stromal'
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)
peak.mtx <- getMatrixFromProject(proj, useMatrix = "PeakMatrix", binarize = T)
saveRDS(peak.mtx, paste0(proj@projectMetadata$outputDirectory, '/peakMat_binarized.rds'))




















