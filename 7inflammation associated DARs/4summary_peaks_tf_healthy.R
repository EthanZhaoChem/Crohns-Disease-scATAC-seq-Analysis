dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(ggplot2)
library(ComplexHeatmap)
library(ggpubr)
library(openxlsx)
source('~/yuzhao1/scripts/helper_archr.R')
addArchRThreads(1)
plot.dir <- '~/yuzhao1/work/atac_gca2024/7dreamlet/plots/'

############################################################
# rowID reference table to peaks
lineage <- 'union'
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1/")
peaks_all <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/', 'peakMtx_unbinarized_rowranges.rds'))

raw_statistics <- readRDS('~/yuzhao1/work/atac_gca2024/7dreamlet/differential_test/statistics_inf_vs_control.rds')
celltype_names <- names(raw_statistics)
outFile.dir <- '~/yuzhao1/work/atac_gca2024/7dreamlet/differential_test/'
plot.dir <- '~/yuzhao1/work/atac_gca2024/7dreamlet/plots/'

####################### prepare peaks in each contrast in each celltype ####################################
peaks_allCelltypes <- list()
direction <- 0

for (celltype in celltype_names) {
  t_statistics <- raw_statistics[[celltype]]
  
  if(direction == 1){
    peaks_test <- peaks_all[as.numeric(rownames(t_statistics)[t_statistics$logFC > 0.5 & t_statistics$adj.P.Val < 0.1]),] 
  }
  
  else if(direction == 0){
    peaks_test <- peaks_all[as.numeric(rownames(t_statistics)[t_statistics$logFC < -0.5 & t_statistics$adj.P.Val < 0.1]),] 
  }
  
  peaks_test_vectors <- c()
  if(length(peaks_test) > 0){
    peaks_test_vectors <- paste0(seqnames(peaks_test), '_', start(peaks_test), '_', end(peaks_test))
  }
  peaks_allCelltypes[[celltype]] <- peaks_test_vectors
}

write.xlsx(peaks_allCelltypes, file = paste0(outFile.dir, 'peaks_inf_vs_control_allCelltypes_FC05_adjP01_healthy.xlsx'))



####################### tfs ####################################
tfs_allCelltypes <- list()
for (celltype in celltype_names) {
  peaks_test_vectors <- peaks_allCelltypes[[celltype]]
  if(length(peaks_test_vectors) > 20){
    peaks_test_vectors <- peaks_allCelltypes[[celltype]]
    tfs <- archr_customized_motif_enrichment(proj,
                                             peakAnnotation = "Motif",
                                             candidate_peaks_vector = peaks_test_vectors)
    tfs_allCelltypes[[celltype]] <- tfs
  }
}
write.xlsx(tfs_allCelltypes, file = paste0(outFile.dir, 'tfs_inf_vs_control_allCelltypes_FC05_adjP01_healthy.xlsx'))





















