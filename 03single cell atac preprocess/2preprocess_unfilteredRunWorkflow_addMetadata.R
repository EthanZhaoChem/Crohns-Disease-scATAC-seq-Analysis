dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/scripts/gca_markers.R')

addArchRThreads(1)
# read files
setwd('~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2')
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2/projdir/")
out.dir <- '/project/gca/yuzhao1/work/atac_gca2024/2preprocess/2preprocess_LSIplots/'
dir.create(out.dir, showWarnings = F)

# ############################# Section0: add metadata ################################
# metadata <- read.table('~/yuzhao1/work/atac_gca2024/0metadata/meta_Ethan_curated_20240211.csv',
#                        header = T, sep = ',')
# # writeLines(metadata$sample, '~/yuzhao1/work/atac_gca2024/0metadata/meta_Ethan_curated_20240211_sampleIDs.txt' )
# 
# proj$Patient_ID <- mapvalues(proj$Sample, from = metadata$sample, metadata$patient)
# proj$Patient_ID_masked <- mapvalues(proj$Sample, from = metadata$sample, metadata$patient_masked)
# proj$biopsy_location <- mapvalues(proj$Sample, from = metadata$sample, metadata$biopsy_location)
# proj$disease_status <- mapvalues(proj$Sample, from = metadata$sample, metadata$disease_status)
# proj$inflammation_status <- mapvalues(proj$Sample, from = metadata$sample, metadata$inflammation_status)
# proj$inflammation_status[proj$inflammation_status==''] <- 'Control'
# 
# ############################# Section1: Dim-Red ################################
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
# proj <- addClusters(
#   input = proj,
#   reducedDims = "IterativeLSI",
#   method = "Seurat",
#   name = "LSI_Clusters_res0.5",
#   force = T,
#   dimsToUse = 2:100,
#   resolution = 0.5,
#   maxClusters = 100,
# )
# 
# proj <- addImputeWeights(proj)
# saveArchRProject(ArchRProj = proj, load = T)
# 
# 
# ############################# Section1: LSI plots ################################


png(paste0(out.dir, 'umap_res0.5.png'),width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_2`
df$cluster_name <- proj$LSI_Clusters_res0.5
plot_df_umap_custom(df, show.label = 'name')
dev.off()


# umap loc
png(paste0(out.dir, 'umap_loc.png'),width = 5600, height = 3600,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_2`
df$cluster_name <- proj$LSI_Clusters_res0.5
p<-plot_df_umap_custom(df, show.label = 'name')+
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
df$cluster_name <- proj$LSI_Clusters_res0.5
p<-plot_df_umap_custom(df, show.label = 'name')+
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
df$cluster_name <- proj$LSI_Clusters_res0.5
p<-plot_df_umap_custom(df, show.label = 'name')+
  facet_wrap(~ inflammation_status) +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(size = 12)
  )
print(p)
dev.off()


# umap sample
png(paste0(out.dir, 'umap_sample.png'),res = 300, height = 15000, width = 5600)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_2`
df$cluster_name <- proj$LSI_Clusters_res0.5
p<-plot_df_umap_custom(df, show.label = 'name')+
  facet_wrap(~ Sample, ncol = 6) +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(size = 12)
  )
print(p)
dev.off()



# feature distribution umap: CA2 APOA4
pdf(paste0(out.dir, 'umap_markerGeneUmaps.pdf'), width = 5, height = 6, pointsize = 1)
p <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "GeneScoreMatrix",
  pal = paletteContinuous("solarExtra"),
  name = c('CA2', 'APOA4', 'CEACAM5'),
  embedding = "LSI_UMAP",
  imputeWeights = getImputeWeights(proj))
print(p)
dev.off()


# sample check proportion
metadata <- read.table('~/yuzhao1/work/atac_gca2024/0metadata/meta_Ethan_curated_20240211.csv',
                       header = T, sep = ',')

png(paste0(out.dir, 'check_sample_proportion.png'),width = 3600, height = 2500,res = 300)
x <- table(as.data.frame(proj@cellColData)[,c('Sample', 'LSI_Clusters_res0.5')])
xx <- as.data.frame(x[, c('C14', 'C18')])
xx$biopsy_location <- mapvalues(xx$Sample, from = metadata$sample, to = metadata$biopsy_location)
xx$patient <- mapvalues(xx$Sample, from = metadata$sample, to = metadata$patient)
xx <- as.data.frame(arrange(xx, biopsy_location, patient))
xx$Sample <- factor(xx$Sample, levels = unique(xx$Sample))

ggplot(xx, aes(x = Sample, y = Freq, fill = LSI_Clusters_res0.5)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c('#1f78b4', '#fccde5'))+
  theme(
    panel.grid.major = element_line(linetype = 'dashed'),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 10, face = 'bold', colour = 'black', angle = 90, vjust = 0.5, hjust=1),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    plot.margin = margin(2, 2, 2, 2, "cm"),
    plot.title = element_blank()
  )
dev.off()






















