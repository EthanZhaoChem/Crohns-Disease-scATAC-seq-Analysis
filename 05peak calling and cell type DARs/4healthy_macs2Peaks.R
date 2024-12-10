dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)
source('~/yuzhao1/scripts/plot.R')
addArchRThreads(1)
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_healthy/")
celltypes <- unique(proj$anno1)
grpath <- '~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_healthy/PeakCalls/anno1/'
out.dir <- '~/yuzhao1/work/atac_gca2024/4peaks/macs2_peaks/healthy/'

# write a bed file for each celltype
for (t.celltype in celltypes){
  peaks <- readRDS(paste0(grpath, t.celltype, '-reproduciblePeaks.gr.rds'))

  peaks_chr <- seqnames(peaks)
  peaks_start <- start(peaks)
  peaks_end <- end(peaks)
  
  df <- data.frame(chr = peaks_chr,
                   start = peaks_start,
                   end = peaks_end)
  
  filename <- paste0(out.dir, t.celltype,'.bed')
  write.table(df, filename, 
              row.names = F, col.names = F, sep="\t", quote=FALSE)
}









