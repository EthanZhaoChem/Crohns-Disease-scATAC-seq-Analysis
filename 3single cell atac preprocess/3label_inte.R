dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/scripts/gca_markers.R')

addArchRThreads(4)
# read files
setwd('~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2')
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2/projdir/")
out.dir <- '/project/gca/yuzhao1/work/atac_gca2024/2preprocess/3label_inte/'
dir.create(out.dir, showWarnings = F)

# ################# Section 1: unconstrained Integration ###########
# 
# # read files
# rna_ref <- readRDS('~/yuzhao1/work/final_GCArna/annotation/rds/gca_combined_final.rds')
# 
# # ## remove old labels
# # meta <- proj@cellColData
# # idx <- grep('predicted', colnames(meta))
# # proj@cellColData[, idx] <- NULL
# 
# df_rna <- rna_ref@meta.data
# df_rna$embedding1 <- data.frame(rna_ref@reductions$pca_umap@cell.embeddings)$UMAP_1
# df_rna$embedding2 <- data.frame(rna_ref@reductions$pca_umap@cell.embeddings)$UMAP_2
# df_rna <- df_rna[, c('embedding1', 'embedding2')]
# 
# df_atac <- data.frame((proj@cellColData))
# df_atac$embedding1 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_1`
# df_atac$embedding2 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_2`
# df_atac <- df_atac[, c('embedding1', 'embedding2')]
# 
# 
# proj<- addGeneIntegrationMatrix(
#   ArchRProj = proj,
#   useMatrix = "GeneScoreMatrix",
#   matrixName = "GeneIntegrationMatrix",
#   reducedDims = "IterativeLSI",
#   seRNA = rna_ref,
#   addToArrow = F,
#   groupList = NULL,
#   groupRNA = "anno1",
#   nameCell = "predictedCell_unconstrained",
#   nameGroup = "predictedGroup_unconstrained",
#   nameScore = "predictedScore_unconstrained",
#   sampleCellsATAC = 20000,
#   sampleCellsRNA = 20000,
#   dimsToUse = 1:50,
#   corCutOff = 0.75,
#   embeddingATAC = df_atac,
#   embeddingRNA = df_rna,
#   force = T
# )
# 
# saveArchRProject(ArchRProj = proj, load = T)



################# Section 2: plot ###########


png(paste0(out.dir, 'umap_predictedGroup_unconstrained.png'),width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_2`
df$cluster_name <- proj$predictedGroup_unconstrained
plot_df_umap_custom(df, show.label = 'name')
dev.off()



















