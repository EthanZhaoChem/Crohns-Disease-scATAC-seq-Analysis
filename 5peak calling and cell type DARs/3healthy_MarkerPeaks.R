dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
# https://github.com/satijalab/seurat/issues/6746
Csparse_validate = "CsparseMatrix_validate"

library(stringr)
library(ArchR)
library(Seurat)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)
source('~/yuzhao1/scripts/plot.R')
addArchRThreads(1)
# addArchRLocking(locking = F)
# # read files
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_healthy/")


############# anno1 #############
set.seed(6)
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  normBy = 'ReadsInTSS',
  maxCells = 500,
  groupBy = "anno1",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

saveRDS(markersPeaks, paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno1_1vsAll.rds'))



markersPeaks <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_anno1_1vsAll.rds'))
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

saveRDS(enrichMotifs, paste0(proj@projectMetadata$outputDirectory, '/', 'enrichMotifs_anno1_1vsAll.rds'))











