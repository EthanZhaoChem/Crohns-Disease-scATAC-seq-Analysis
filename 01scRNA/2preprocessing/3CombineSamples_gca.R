dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(Seurat)
library(Matrix)
library(magrittr)
library(dplyr)
library(plyr)
SampleIDs <- read.csv("~/yuzhao1/work/final_GCArna/metadata/SampleIDs.csv")
SampleIDs <- SampleIDs$SampleID
meta <- read.table(file = '~/yuzhao1/work/final_GCArna/metadata/meta_Ethan_curated_20221108.csv', header = T, sep = ',')
meta$Sample_ID[which(meta$Sample_ID%in%SampleIDs == 0)]
SampleIDs[which(SampleIDs%in%meta$Sample_ID == 0)]

mtx_outdir <- '~/yuzhao1/work/final_GCArna/preprocessing/matrix_removedAmbientRNA/'
gene_outdir <- '~/yuzhao1/work/final_GCArna/preprocessing/gene_removedAmbientRNA/'
barcodes_outdir <- '~/yuzhao1/work/final_GCArna/preprocessing/barcodes_outdir/'
doublet_outdir <- '~/yuzhao1/work/final_GCArna/preprocessing/doublet_scores/'
out.srt.obj.path <- '~/yuzhao1/work/final_GCArna/preprocessing/GCArna_all80samples_removedAmbientRNA_calculatedDoubletScores_seurat.rds'

# read each sample matrix and combine as a seurat object
srat = list()
sample.col <- c()
doubletscore.all <- c()
for (i in 1:80) {
  SampleID <- SampleIDs[[i]]
  mtx <- readMM(paste0(mtx_outdir, SampleID,'.mtx'))
  genes <- read.table(file = paste0(gene_outdir, SampleID,'.tsv'),
                      sep = '\t', header = F) %>% unlist()
    
  
  cells <- read.table(file = paste0(barcodes_outdir, SampleID,'.tsv'),
                      sep = '\t', header = F) %>% unlist()
  doublet.scores <- read.table(file = paste0(doublet_outdir, SampleID,'.csv'),
                      sep = ',', header = F) %>% unlist()
  names(doublet.scores) <- cells
  colnames(mtx) <- cells
  rownames(mtx) <- genes
  
  doubletscore.all <- c(doubletscore.all, doublet.scores)
  sample.col <- c(sample.col, rep(SampleID, length(cells)))
  srat[[SampleID ]] <- mtx
}


# Combine all count matricies into one matrix
srat.mtx = do.call(cbind,srat[1:80])
srat.obj = CreateSeuratObject(srat.mtx)

# add metadata from processing to seurat object
srat.obj$Sample_ID <- sample.col
srat.obj$Doublet_score <- doubletscore.all

# add metadata from meta file to seurat object
metaColsTodo <- c("Patient_ID", "biopsy_location", "disease_status", "inflammation_status", "submission_date")
for (column.ToDo in metaColsTodo) {
  temp.col <- mapvalues(srat.obj$Sample_ID, from = meta$Sample_ID, to = meta[,column.ToDo])
  srat.obj@meta.data[[column.ToDo]] <- temp.col
}

saveRDS(srat.obj, out.srt.obj.path)














