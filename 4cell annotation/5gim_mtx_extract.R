# filter again based on unconstrained predicted labels: at least they don't work based on gene score
dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(ggplot2)
library(Matrix)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/scripts/gca_markers.R')

addArchRThreads(1)
pathToMacs2 <- '/home/yuzhao1/.local/bin/macs2'

# read files
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1/")

cat(paste0('Geting gim: ', '\n'))
gim <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "GeneIntegrationMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

# cat(paste0('Imputing GIM: ', '\n'))
# gim_imputed <- imputeMatrix(
#   mat = gim@assays@data$GeneIntegrationMatrix,
#   imputeWeights = getImputeWeights(proj),
#   threads = getArchRThreads(),
#   verbose = T,
#   logFile = createLogFile("imputeMatrix")
# )
# 
# rownames(gim_imputed) <- gim@elementMetadata$name

cat(paste0('Saving gim\n'))
saveRDS(gim, paste0(proj@projectMetadata$outputDirectory, '/GeneIntegrationMatrix.rds'))




