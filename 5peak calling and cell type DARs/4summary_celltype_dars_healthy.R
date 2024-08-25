dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')

library(stringr)
library(ArchR)
library(Seurat)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)

source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/scripts/helper_archr.R')
addArchRThreads(12)
addArchRLocking(locking = F)

# # read files
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_healthy/")
celltypes <- unique(proj$anno1)
markersPeaks <- readRDS(paste0(proj@projectMetadata$outputDirectory,
                               '/', 'MarkersPeaks_anno1_1vsAll.rds'))
                        
# customize
cutoff_FDR <- 0.1
cutoff_Log2FC <- 0.5
out.dir <- '~/yuzhao1/work/atac_gca2024/4peaks/DARs/healthy_union_anno1_peaks/FDR0_1FC0_5/'

##########################
peaks_filtered <- list()

df_stat_ByCluster <- archr_helper_markerPeaks_converter_multiple(markersPeaks)

for (cluster.name in celltypes){ 
  df_stat <- df_stat_ByCluster[[cluster.name]]
  df_stat_significant <- df_stat[which(df_stat$FDR < cutoff_FDR & df_stat$Log2FC > cutoff_Log2FC), ]
  peaks_filtered[[cluster.name]] <- rownames(df_stat_significant)
}


# write a bed file for each celltype
for (i in 1:length(peaks_filtered)){
  t.celltype <- names(peaks_filtered)[[i]]
  nPeaks <- length(peaks_filtered[[i]])
  peaks <- peaks_filtered[[i]]
  
  if(nPeaks < 1){
    next
  }
  
  peaks_chr <- peaks %>% strsplit(split = ':', fixed=T) %>% sapply(.,`[[`,1)
  
  start_end <- peaks %>% strsplit(split = ':', fixed=T) %>% sapply(.,`[[`,2) 
  peaks_start <- start_end%>% strsplit(split = '-', fixed=T) %>% sapply(.,`[[`,1)
  peaks_end <- start_end%>% strsplit(split = '-', fixed=T) %>% sapply(.,`[[`,2)
  
  df <- data.frame(chr = peaks_chr,
                   start = peaks_start,
                   end = peaks_end)
  
  filename <- paste0(out.dir, t.celltype,'.bed')
  write.table(df, filename, 
              row.names = F, col.names = F, sep="\t", quote=FALSE)
}


# save summarized peaks
peaks_filtered_all <- unique(unlist(peaks_filtered))
peaks <- peaks_filtered_all
peaks_chr <- peaks %>% strsplit(split = ':', fixed=T) %>% sapply(.,`[[`,1)
start_end <- peaks %>% strsplit(split = ':', fixed=T) %>% sapply(.,`[[`,2) 
peaks_start <- start_end%>% strsplit(split = '-', fixed=T) %>% sapply(.,`[[`,1)
peaks_end <- start_end%>% strsplit(split = '-', fixed=T) %>% sapply(.,`[[`,2)

peaks_filtered_all_bed <- data.frame(chr = peaks_chr,
                                     start = peaks_start,
                                     end = peaks_end)
write.table(peaks_filtered_all_bed,
            paste0(out.dir,'all_celltype_Specific.bed'), 
            row.names = F, col.names = F, sep="\t", quote=FALSE)


##########################
# build peak-celltype table
peaks_anno1_list <- list()
dir_peak_dar_anno1 <- '~/yuzhao1/work/atac_gca2024/4peaks/DARs/healthy_union_anno1_peaks/FDR0_1FC0_5/'
celltypes <- names(gca_colors_atac_union_anno1)
for(ct in celltypes){
  xx <- read.table(paste0(dir_peak_dar_anno1, ct, '.bed'))
  peaks_anno1_list[[ct]] <- paste0(xx$V1, '_', xx$V2, '_', xx$V3)
}

peaks_anno1_df <- data.frame(peak = unique(unlist(peaks_anno1_list)), cellType = 'TBD')
rownames(peaks_anno1_df) <- peaks_anno1_df$peak
for (ct in celltypes) {
  peaks_ct_all <- peaks_anno1_list[[ct]]
  peaks_anno1_df[peaks_ct_all, 'cellType'] <- paste0(peaks_anno1_df[peaks_ct_all, 'cellType'], '&', ct)
}
peaks_anno1_df$cellType <- gsub('TBD&', '', peaks_anno1_df$cellType)
rownames(peaks_anno1_df) <- NULL

saveRDS(peaks_anno1_df, '~/yuzhao1/work/atac_gca2024/4peaks/DARs/healthy_union_anno1_peaks/FDR0_1FC0_5_peak_cellType_table.rds')





