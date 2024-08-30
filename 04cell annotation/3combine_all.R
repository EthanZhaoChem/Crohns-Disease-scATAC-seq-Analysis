dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/work/final_GCArna/scripts/gca_colors.R')
source('~/yuzhao1/work/final_GCArna/scripts/gca_markers.R')

addArchRThreads(1)
# read files
# proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2/projdir")
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1/")
out.dir <- '~/yuzhao1/work/atac_gca2024/3annotation/plots/filtered1/'
dir.create(out.dir, showWarnings = F, recursive = T)
addArchRLocking(locking = F)

# # ############################# Section0: subset ################################
# proj.epithelial <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_epithelial2/")
# proj.immune <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_immune2/")
# proj.stromal <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_stromal2/")
# 
# 
# df_annotation <- data.frame(cellnames = c(proj.epithelial$cellNames, proj.immune$cellNames, proj.stromal$cellNames),
#                             anno1s = c(proj.epithelial$anno1, proj.immune$anno1, proj.stromal$anno1))
# 
# # proj.filtered <- subsetArchRProject(
# #   ArchRProj = proj,
# #   cells = df_annotation$cellnames,
# #   outputDirectory = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1",
# #   dropCells = TRUE,
# #   logFile = NULL,
# #   threads = getArchRThreads(),
# #   force = T
# # )
# 
# idx <- match(df_annotation$cellnames, proj$cellNames)
# proj$anno1 <- 'NA'
# proj$anno1[idx] <- df_annotation$anno1s
# 
# saveArchRProject(ArchRProj = proj, load = T)
# 
# ############################ Section1: Dim-Red ################################
# # Dim Reduction
# proj <- addIterativeLSI(
#   ArchRProj = proj,
#   useMatrix = "TileMatrix",
#   name = "IterativeLSI",
#   iterations = 6,
#   dimsToUse = 2:100,
#   varFeatures = 50000,
#   clusterParams = list(
#     sampleCells = 50000,
#     resolution = c(0.5, 1, 1.5, 2, 2),
#     n.start = 10
#   ),
#   sampleCellsPre = 50000,
#   force = TRUE
# )
# 
# saveArchRProject(ArchRProj = proj, load = T)
# 
# 
# # LSI umap and clusters
# proj <- addUMAP(
#   ArchRProj = proj,
#   reducedDims = "IterativeLSI",
#   name = "LSI_UMAP",
#   nNeighbors = 50,
#   minDist = 0.5,
#   force = T,
#   dimsToUse = 2:100,
#   metric = "cosine"
# )
# 
# 
# proj <- addImputeWeights(
#   ArchRProj = proj,
#   reducedDims = "IterativeLSI",
#   dimsToUse = 2:50,
#   sampleCells = 50000)
# 
# saveArchRProject(ArchRProj = proj, load = T)


################# Section 2: lineage constrained Integration ###########
# read files
rna_ref <- readRDS('~/yuzhao1/work/final_GCArna/annotation/rds/gca_combined_final.rds')

# constrain by lineage
rna_epithelial <- colnames(rna_ref)[grep('epithelial', rna_ref$category1)]
rna_immune <- colnames(rna_ref)[grep('immune', rna_ref$category1)]
rna_stromal <- colnames(rna_ref)[grep('stromal', rna_ref$category1)]

atac_epithelial <- proj$cellNames[proj$category1 == 'epithelial']
atac_immune <- proj$cellNames[proj$category1 == 'immune']
atac_stromal <- proj$cellNames[proj$category1 == 'stromal']

groupList <- SimpleList(
  epithelial = SimpleList(
    ATAC = atac_epithelial,
    RNA = rna_epithelial
  ),
  immune = SimpleList(
    ATAC = atac_immune,
    RNA = rna_immune
  ),
  stromal = SimpleList(
    ATAC = atac_stromal,
    RNA = rna_stromal
  )
)

# embedding
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
  addToArrow = T,
  groupList = groupList,
  groupRNA = "anno1",
  nameCell = "predictedCell_lineage_constrained",
  nameGroup = "predictedGroup_lineage_constrained",
  nameScore = "predictedScore_lineage_constrained",
  sampleCellsATAC = 20000,
  sampleCellsRNA = 20000,
  dimsToUse = 2:50,
  corCutOff = 0.75,
  embeddingATAC = df_atac,
  embeddingRNA = df_rna,
  force = T
)

saveArchRProject(ArchRProj = proj, load = T)

# # ############################# Section1: LSI plots ################################
png(paste0(out.dir, 'umap_anno1_number.png'),width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_2`
df$cluster_name <- proj$anno1
plot_df_umap_custom(df, show.label = 'number', custom_colors = gca_colors_atac_union_anno1)
dev.off()


png(paste0(out.dir, 'umap_predicted_lineage_constrained.png'),width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_2`
df$cluster_name <- proj$predictedGroup_lineage_constrained
plot_df_umap_custom(df, show.label = 'name')
dev.off()


# umap loc
png(paste0(out.dir, 'umap_loc_number.png'),width = 5600, height = 3600,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_2`
df$cluster_name <- proj$anno1
p<-plot_df_umap_custom(df, show.label = 'number', custom_colors = gca_colors_atac_union_anno1)+
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
p<-plot_df_umap_custom(df, show.label = 'name', custom_colors = gca_colors_atac_union_anno1)+
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
p<-plot_df_umap_custom(df, show.label = 'name', custom_colors = gca_colors_atac_union_anno1)+
  facet_wrap(~ inflammation_status) +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(size = 12)
  )
print(p)
dev.off()


# umap sample
png(paste0(out.dir, 'umap_sample.png'),res = 300, height = 13000, width = 5600)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_2`
df$cluster_name <- proj$anno1
p<-plot_df_umap_custom(df, show.label = 'name', custom_colors = gca_colors_atac_union_anno1)+
  facet_wrap(~ Sample, ncol = 6) +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(size = 12)
  )
print(p)
dev.off()




















