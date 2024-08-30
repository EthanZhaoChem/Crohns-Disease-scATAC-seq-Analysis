dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/work/final_GCArna/scripts/gca_colors.R')
source('~/yuzhao1/work/final_GCArna/scripts/gca_markers.R')

addArchRThreads(1)
# read files
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_stromal2/")
out.dir <- '~/yuzhao1/work/atac_gca2024/3annotation/plots/stromal2/'
dir.create(out.dir, showWarnings = F, recursive = T)
addArchRLocking(locking = F)


# # ############################# Section 0: annotation ################################

## filter from first iteration
# df_annotation_res0.5 <- list(
#   'C1' ='Glial',
#   'C2' ='Endothelium',
#   'C3' ='Endothelium',
#   'C4' ='Endothelium',
#   'C5' ='Endothelium',
#   'C6' ='Endothelium',
#   'C7' ='Endothelium',
#   'C8' ='Fibroblast',
#   'C9' ='Pericyte',
#   'C10' ='Fibroblast',
#   'C11' ='Fibroblast',
#   'C12' ='Fibroblast',
#   'C13' ='Fibroblast',
#   'C14' ='Fibroblast',
#   'C15' ='LQ',
#   'C16' ='Pericyte',
#   'C17' ='Pericyte')
# 
# proj$anno1 <- unlist(mapvalues(proj$LSI_Clusters_res0.5, names(df_annotation_res0.5), df_annotation_res0.5))
# 
# # remove cells identified as not belonging to the correct lineage during unconstrained matrix integration
# proj$anno1[!(proj$predictedGroup_unconstrained %in% unique(rna_ref$anno1))] <- 'LQ'
# saveArchRProject(ArchRProj = proj, load = T)
# 
# # subset
# idxPass <- which(!(proj$anno1 %in% c("LQ")))
# cellsPass.stromal <- proj$cellNames[idxPass]
# proj.stromal <- subsetArchRProject(
#   ArchRProj = proj,
#   cells = cellsPass.stromal,
#   outputDirectory = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_stromal2",
#   dropCells = TRUE,
#   logFile = NULL,
#   threads = getArchRThreads(),
#   force = T
# )


# ############################# Section1: Dim-Red ################################
# Dim Reduction
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 6,
  dimsToUse = 2:100,
  varFeatures = 30000,
  clusterParams = list(
    sampleCells = 5000,
    resolution = c(0.5, 1, 1.5, 2, 2),
    n.start = 10
  ),
  sampleCellsPre = 5000,
  force = TRUE
)

proj <- addImputeWeights(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  dimsToUse = 2:30,
  corCutOff = 0.75,
  sampleCells = 5000)


saveArchRProject(ArchRProj = proj, load = T)

# 
# LSI umap and clusters
# proj <- addUMAP(
#   ArchRProj = proj,
#   reducedDims = "IterativeLSI",
#   name = "LSI_UMAP",
#   nNeighbors = 50,
#   minDist = 0.5,
#   force = T,
#   dimsToUse = 2:30,
#   metric = "cosine"
# )
# 
# proj <- addClusters(
#   input = proj,
#   reducedDims = "IterativeLSI",
#   method = "Seurat",
#   name = "LSI_Clusters_res0.5",
#   force = T,
#   dimsToUse = 2:30,
#   resolution = 0.5,
#   maxClusters = 100,
# )
# 
# saveArchRProject(ArchRProj = proj, load = T)



################# Section 2: unconstrained Integration ###########
# read files
rna_ref <- readRDS('~/yuzhao1/work/final_GCArna/annotation/rds/stromal2.rds')

# ## remove old labels
# meta <- proj@cellColData
# idx <- grep('predicted', colnames(meta))
# proj@cellColData[, idx] <- NULL

# lineage constrained
df_rna <- rna_ref@meta.data
df_rna$embedding1 <- data.frame(rna_ref@reductions$pca_umap@cell.embeddings)$UMAP_1
df_rna$embedding2 <- data.frame(rna_ref@reductions$pca_umap@cell.embeddings)$UMAP_2
df_rna <- df_rna[, c('embedding1', 'embedding2')]

df_atac <- data.frame((proj@cellColData))
df_atac$embedding1 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_1`
df_atac$embedding2 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_2`
df_atac <- df_atac[, c('embedding1', 'embedding2')]


proj<- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = rna_ref,
  addToArrow = F,
  groupList = NULL,
  groupRNA = "anno1",
  nameCell = "predictedCell_lineage_constrained",
  nameGroup = "predictedGroup_lineage_constrained",
  nameScore = "predictedScore_lineage_constrained",
  sampleCellsATAC = 10000,
  sampleCellsRNA = 10000,
  dimsToUse = 1:30,
  embeddingATAC = df_atac,
  embeddingRNA = df_rna,
  force = T
)


saveArchRProject(ArchRProj = proj, load = T)







# # ############################# Section3: LSI plots ################################
png(paste0(out.dir, 'umap_anno1.png'),width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_2`
df$cluster_name <- proj$anno1
plot_df_umap_custom(df, show.label = 'name', custom_colors = gca_colors_atac_stromal_anno1)
dev.off()

png(paste0(out.dir, 'umap_transferred_anno1.png'),width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_2`
df$cluster_name <- proj$transferred_anno1
plot_df_umap_custom(df, show.label = 'name')
dev.off()


png(paste0(out.dir, 'umap_transferred_anno1_group.png'),width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_2`
df$cluster_name <- proj$transferred_anno1_group
plot_df_umap_custom(df, show.label = 'name')
dev.off()

png(paste0(out.dir, 'umap_predictedGroup_unconstrained.png'),width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_2`
df$cluster_name <- proj$predictedGroup_unconstrained
plot_df_umap_custom(df, show.label = 'name')
dev.off()

png(paste0(out.dir, 'umap_predictedGroup_lineage_constrained.png'),width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_2`
df$cluster_name <- proj$predictedGroup_lineage_constrained
plot_df_umap_custom(df, show.label = 'name')
dev.off()

# umap location
png(paste0(out.dir, 'umap_location.png'),width = 5600, height = 3600,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_2`
df$cluster_name <- proj$anno1
p<-plot_df_umap_custom(df, show.label = 'name', custom_colors = gca_colors_atac_stromal_anno1)+
  facet_wrap(~ biopsy_location) +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(size = 12)
  )
print(p)
dev.off()

# umap disease
png(paste0(out.dir, 'umap_disease.png'),width = 5600, height = 3600,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_2`
df$cluster_name <- proj$anno1
p<-plot_df_umap_custom(df, show.label = 'name', custom_colors = gca_colors_atac_stromal_anno1)+
  facet_wrap(~ disease_status) +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(size = 12)
  )
print(p)
dev.off()

# umap inf
png(paste0(out.dir, 'umap_inf.png'),width = 5600, height = 6000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_2`
df$cluster_name <- proj$anno1
p<-plot_df_umap_custom(df, show.label = 'name', custom_colors = gca_colors_atac_stromal_anno1)+
  facet_wrap(~ inflammation_status) +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(size = 12)
  )
print(p)
dev.off()


# umap sample
png(paste0(out.dir, 'umap_sample.png'),res = 300, height = 7200, width = 4500)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_2`
df$cluster_name <- proj$anno1
p<-plot_df_umap_custom(df, show.label = 'na', custom_colors = gca_colors_atac_stromal_anno1)+
  facet_wrap(~ Sample, ncol = 6) +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(size = 12)
  )
print(p)
dev.off()



















