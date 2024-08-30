dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(ggplot2)
library(dplyr)
library(plyr)
library(purrr)
library(stringr)
library(harmony)
library(Seurat)
source('~/yuzhao1/scripts/plot.R')
seurat <- readRDS('~/gca/yuzhao1/work/final_GCArna/preprocessing/GCArna_all80samples_removedAmbientRNA_calculatedDoubletScores_seurat_filtered.rds')

# 1. subset (removing confusing samples, we are not sure whether it's cause by low quality of samples or metaplasia)
samples.confusing <- c('HA14-AC', 'HA16-AC', 'HA17-AC', 'HA30-AC', 'HA14-TI',
                       'HA19-TI-non', 'HA70TI-inf', 'HA19-TI-adj', 'HA19-TI-inf',
                       'HA56TI-inf')
metadata <- read.table('~/yuzhao1/work/final_GCAatac/0metadata/meta_Ethan_curated_20230316.csv', 
                       header = T, sep = ',')

# get a seurat subset of validated samples
validate_samples <- setdiff(unique(metadata$Sample_ID), samples.confusing)
paired_cells <- Cells(seurat)[which(seurat$Sample_ID%in%validate_samples)]
seurat_paired <- subset(seurat, cells = paired_cells)


# 2. correct swapped labels in the subset
df_correct_swappedLabels <- data.frame(
  old = c('HA02-AC', 'HA07-AC', 'HA08-AC', 'HA11-AC', 'HA02-TI', 'HA07-TI', 'HA08-TI', 'HA11-TI'),
  new = c('HA02-TI', 'HA07-TI', 'HA08-TI', 'HA11-TI', 'HA02-AC', 'HA07-AC', 'HA08-AC', 'HA11-AC')
)
seurat_paired$Sample_ID <- mapvalues(seurat_paired$Sample_ID, from=df_correct_swappedLabels$old,
          to=df_correct_swappedLabels$new)
seurat_paired$submission_date <- NULL
seurat_paired$Patient_ID <- mapvalues(seurat_paired$Sample_ID, from = metadata$Sample_ID, to = metadata$Patient_ID)
seurat_paired$disease_status <- mapvalues(seurat_paired$Sample_ID, from = metadata$Sample_ID, to = metadata$disease_status)
seurat_paired$inflammation_status <- mapvalues(seurat_paired$Sample_ID, from = metadata$Sample_ID, to = metadata$inflammation_status)
seurat_paired$biopsy_location <- mapvalues(seurat_paired$Sample_ID, from = metadata$Sample_ID, to = metadata$biopsy_location)


# 3. run analysis workflow
seurat <- seurat_paired
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat <- CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
seurat$CC.Difference <- seurat$S.Score - seurat$G2M.Score

seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
seurat <- ScaleData(seurat, vars.to.regress = c('percent.mt', 'nCount_RNA', 'CC.Difference'))
seurat <- RunPCA(seurat, npcs = 50)
ElbowPlot(seurat,ndims = 50)


# running harmony will make non-functioning cells into colonocyte
# seurat <- RunHarmony(seurat, group.by.vars = 'Patient_ID', max.iter.harmony = 20) 
seurat <- FindNeighbors(seurat, reduction = 'pca', dims = 1:30)

# seurat <- RunUMAP(seurat,  dims = 1:50, reduction = 'harmony', reduction.name = 'harmony_umap',
#                   min.dist = 0.5, n.neighbors = 50, seed.use = 5,
#                   reduction.key = 'UMAP_')

seurat <- RunUMAP(seurat,  dims = 1:30, reduction = 'pca', reduction.name = 'pca_umap',
                  min.dist = 0.5, n.neighbors = 30, seed.use = 5,
                  reduction.key = 'UMAP_')

seurat <- FindClusters(seurat, resolution = 0.5)
seurat$seurat_clusters_res0.5 <- Idents(seurat)

seurat <- FindClusters(seurat, resolution = 0.8)
seurat$seurat_clusters_res0.8 <- Idents(seurat)

seurat <- FindClusters(seurat, resolution = 1)
seurat$seurat_clusters_res1 <- Idents(seurat)

seurat <- FindClusters(seurat, resolution = 1.5)
seurat$seurat_clusters_res1.5 <- Idents(seurat)


saveRDS(seurat,
        '~/gca/yuzhao1/work/final_GCArna/preprocessing/GCArna_55PairedSamples_removedAmbientRNA_calculatedDoubletScores_seurat_filtered_processed.rds')























