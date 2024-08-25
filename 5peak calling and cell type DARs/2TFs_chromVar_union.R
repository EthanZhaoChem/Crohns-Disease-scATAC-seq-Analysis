dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)

source('~/yuzhao1/scripts/plot.R')
addArchRThreads(36)
addArchRLocking(locking = F)

# # read files
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1/")

out.dir <- '~/yuzhao1/work/atac_gca2024/4peaks/plots/'
dir.create(out.dir, showWarnings = F, recursive = T)


########################## 1. calculation ###############################
proj <- addMotifAnnotations(ArchRProj = proj,  motifSet = "cisbp", force = TRUE)
proj <- addBgdPeaks(proj, force = TRUE)
## save
saveArchRProject(ArchRProj = proj, load = T)


proj <- addDeviationsMatrix(
  ArchRProj = proj,
  matrixName = 'cisbp',
  out = c("z", "deviations"),
  binarize = FALSE,
  verbose = TRUE,
  force = T
)

saveArchRProject(ArchRProj = proj, load = T)





