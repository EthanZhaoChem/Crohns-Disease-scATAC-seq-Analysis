dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(ggplot2)
library(dplyr)
library(plyr)
library(purrr)
library(stringr)
library(harmony)
library(Seurat)
source('~/yuzhao1/scripts/plot.R')

# 1. plot data
seurat <- readRDS('~/gca/yuzhao1/work/final_GCArna/preprocessing/GCArna_55PairedSamples_removedAmbientRNA_calculatedDoubletScores_seurat_filtered_processed.rds')
out.dir <- '/project/gca/yuzhao1/work/final_GCArna/preprocessing/6plot_corrected/'

xx <- FindMarkers(seurat, ident.1 = '8', ident.2 = '0', group.by = 'seurat_clusters_res0.5')
xx <- FindMarkers(seurat, ident.1 = '24', ident.2 = '0',group.by = 'seurat_clusters_res0.5')

# umap
png(paste0(out.dir, 'umap_res0.5.png'),res = 300, height = 3600, width = 3000)
df <- data.frame(seurat@meta.data)
df$embedding1 <- data.frame(seurat@reductions$pca_umap@cell.embeddings)$UMAP_1
df$embedding2 <- data.frame(seurat@reductions$pca_umap@cell.embeddings)$UMAP_2
df$cluster_name <- seurat$seurat_clusters_res0.5
p<-plot_df_umap_custom(df, show.label = 'name')
print(p)
dev.off()

# umap loc
png(paste0(out.dir, 'umap_loc_res0.5.png'),res = 300, height = 3600, width = 5600)
df <- data.frame(seurat@meta.data)
df$embedding1 <- data.frame(seurat@reductions$pca_umap@cell.embeddings)$UMAP_1
df$embedding2 <- data.frame(seurat@reductions$pca_umap@cell.embeddings)$UMAP_2
df$cluster_name <- seurat$seurat_clusters_res0.5
p<-plot_df_umap_custom(df, show.label = 'name')+
  facet_wrap(~ biopsy_location) +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(size = 12)
  )
print(p)
dev.off()

# umap disease
png(paste0(out.dir, 'umap_disease_res0.5.png'),res = 300, height = 3600, width = 5600)
df <- data.frame(seurat@meta.data)
df$embedding1 <- data.frame(seurat@reductions$pca_umap@cell.embeddings)$UMAP_1
df$embedding2 <- data.frame(seurat@reductions$pca_umap@cell.embeddings)$UMAP_2
df$cluster_name <- seurat$seurat_clusters_res0.5
p<-plot_df_umap_custom(df, show.label = 'name')+
  facet_wrap(~ disease_status) +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(size = 12)
  )
print(p)
dev.off()

# umap inf
png(paste0(out.dir, 'umap_inf_res0.5.png'),res = 300, height = 6000, width = 5600)
df <- data.frame(seurat@meta.data)
df$embedding1 <- data.frame(seurat@reductions$pca_umap@cell.embeddings)$UMAP_1
df$embedding2 <- data.frame(seurat@reductions$pca_umap@cell.embeddings)$UMAP_2
df$cluster_name <- seurat$seurat_clusters_res0.5
p<-plot_df_umap_custom(df, show.label = 'name')+
  facet_wrap(~ inflammation_status) +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(size = 12)
  )
print(p)
dev.off()

# pca umap sample
png(paste0(out.dir, 'pca_umap_sample_res0.5.png'),res = 300, height = 13000, width = 5600)
df <- data.frame(seurat@meta.data)
df$embedding1 <- data.frame(seurat@reductions$pca_umap@cell.embeddings)$UMAP_1
df$embedding2 <- data.frame(seurat@reductions$pca_umap@cell.embeddings)$UMAP_2
df$cluster_name <- seurat$seurat_clusters_res0.5
p<-plot_df_umap_custom(df, show.label = 'name')+
  facet_wrap(~ Sample_ID, ncol = 6) +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(size = 12)
  )
print(p)
dev.off()

# umap gene feature dist
gene <- 'DEFA5'
png(paste0(out.dir, 'gene_umap_', gene,'.png'), res = 300, height = 2000, width = 2000)
feature <- gene
df <- data.frame(seurat@meta.data)
df$embedding1 <- data.frame(seurat@reductions$pca_umap@cell.embeddings)$UMAP_1
df$embedding2 <- data.frame(seurat@reductions$pca_umap@cell.embeddings)$UMAP_2
df$feature_to_plot <- FetchData(seurat, feature) %>% unlist()
p <- plot_df_umap_custom(df, plot_feature = T)
print(p)
dev.off()

# umap CA2
png(paste0(out.dir, 'umap_CA2.png'), res = 300, height = 2000, width = 2000)
feature <- 'CA2'
df <- data.frame(seurat@meta.data)
df$embedding1 <- data.frame(seurat@reductions$pca_umap@cell.embeddings)$UMAP_1
df$embedding2 <- data.frame(seurat@reductions$pca_umap@cell.embeddings)$UMAP_2
df$feature_to_plot <- FetchData(seurat, feature) %>% unlist()
p <- plot_df_umap_custom(df, plot_feature = T)
print(p)
dev.off()

png(paste0(out.dir, 'umap_CA2_sample.png'), res = 300, height = 13000, width = 5600)
feature <- 'CA2'
df <- data.frame(seurat@meta.data)
df$embedding1 <- data.frame(seurat@reductions$pca_umap@cell.embeddings)$UMAP_1
df$embedding2 <- data.frame(seurat@reductions$pca_umap@cell.embeddings)$UMAP_2
df$feature_to_plot <- FetchData(seurat, feature) %>% unlist()
p <- plot_df_umap_custom(df, plot_feature = T)+
  facet_wrap(~ Sample_ID, ncol = 6) +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(size = 12)
  )
print(p)
dev.off()

# umap APOA4
png(paste0(out.dir, 'umap_APOA4.png'),res = 300, height = 2000, width = 2000)
feature <- 'APOA4'
df <- data.frame(seurat@meta.data)
df$embedding1 <- data.frame(seurat@reductions$pca_umap@cell.embeddings)$UMAP_1
df$embedding2 <- data.frame(seurat@reductions$pca_umap@cell.embeddings)$UMAP_2
df$feature_to_plot <- FetchData(seurat, feature) %>% unlist()
p <- plot_df_umap_custom(df, plot_feature = T)
print(p)
dev.off()


# sample check proportion
metadata <- read.table('~/yuzhao1/work/final_GCAatac/0metadata/meta_Ethan_curated_20230316.csv', 
                       header = T, sep = ',')

png(paste0(out.dir, 'check_sample_proportion.png'),width = 3600, height = 2500,res = 300)
x <- table(seurat@meta.data[,c('Sample_ID', 'seurat_clusters_res0.5')])
xx <- as.data.frame(x[, c('0', '13')])
xx$biopsy_location <- mapvalues(xx$Sample_ID, from = metadata$Sample_ID, to = metadata$biopsy_location)
xx$patient <- mapvalues(xx$Sample_ID, from = metadata$Sample_ID, to = metadata$Patient_ID)
xx <- as.data.frame(arrange(xx, biopsy_location, patient))
xx$Sample_ID <- factor(xx$Sample_ID, levels = unique(xx$Sample_ID))

ggplot(xx, aes(x = Sample_ID, y = Freq, fill = seurat_clusters_res0.5)) +
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




