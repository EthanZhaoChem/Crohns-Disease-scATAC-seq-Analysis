dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(plyr)
library(dplyr)
library(ArchR)
source('~/yuzhao1/work/atac_gca2024/22abc/helper_abc.R')
out.dir <- '~/yuzhao1/work/atac_gca2024/22abc/'

df_abc <- readRDS('~/yuzhao1/work/atac_gca2024/22abc/abc_df_ibd.rds')
peakset <- readRDS('~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1/peakMtx_unbinarized_rowranges.rds')
peaknames <- paste0(seqnames(peakset), '_', start(peakset), '_', end(peakset))

caPeak_abcGene <- c()
count <- 0
df_sub_list <- list()
for (tmppeak in peaknames) {
  df_abc_sub <- helper_abcMaxLine(name_chr_pos_hg38 = tmppeak, df_abc = df_abc)
  if (length(df_abc_sub) > 0) {
    count <- count+1
    if(count%%100 == 0){
      print(paste0(count, ' finished.. \n'))
    }
    
    df_abc_sub$archr_peak <- tmppeak
    df_sub_list[[count]] <- df_abc_sub
  }else{
    next
  }
}

df <- bind_rows(df_sub_list)
rownames(df) <- df$archr_peak
saveRDS(df, '~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1/peakMtx_unbinarized_rowranges_abc_annotated.rds')