dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP150.GRCh38) 
library(BSgenome.Hsapiens.UCSC.hg38)     
library(BSgenome)
library(plyr)
library(dplyr)
library(stringr)
out.dir <- '~/yuzhao1/work/atac_gca2024/26motif_disruption3/3results_summary/'
reslts.chunk.dir <- '~/yuzhao1/work/atac_gca2024/26motif_disruption3/2results_allHumanTF/'
results <- list()

# 1. read each chunk and clean for dataframe summary
for (bedChunkID in 1:200) {
  if(bedChunkID%%10 ==0){
    cat(bedChunkID, '\n')
  }
  
  df <- readRDS(paste0(reslts.chunk.dir, bedChunkID, '.rds')) %>% as.data.frame(row.names = NULL)
  df <- df[, c("seqnames", "start", "end", "width", "strand", "SNP_id", "REF", "ALT",
               "varType", "motifPos", "geneSymbol", "scoreRef", "scoreAlt", "alleleDiff", "effect" )]
  df$snp.score <- abs(df$scoreAlt - df$scoreRef)
  results[[bedChunkID]] <- df
}

motif.breaks.all <- bind_rows(results)
motif.breaks.all <- motif.breaks.all[order(motif.breaks.all$snp.score, decreasing = T), ]
motif.breaks.all_unique <- motif.breaks.all[!duplicated(paste0(motif.breaks.all$SNP_id, motif.breaks.all$geneSymbol)), ]
saveRDS(motif.breaks.all_unique, paste0(out.dir, 'motif.breaks.all_unique.rds'))

# # 2. read and save raw results
# results_raw <- NULL
# for (bedChunkID in 1:200) {
#   if(bedChunkID%%10 ==0){
#     cat(bedChunkID, '\n')
#   }
#   xx <- readRDS(paste0(reslts.chunk.dir, bedChunkID, '.rds'))
#   
#   # to initialize the correct format
#   if(bedChunkID == 1){
#     results_raw <- xx
#     next
#   }
#   
#   results_raw <- c(results_raw, xx)
# }


# save
# saveRDS(results_raw, paste0(out.dir, 'results_raw.rds'))







