# supposed to calculate all snps
dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP150.GRCh38) 
library(BSgenome.Hsapiens.UCSC.hg38)     
library(BSgenome)
library(plyr)
library(dplyr)
library(stringr)
bedDir <- '~/yuzhao1/work/atac_gca2024/26motif_disruption3/1rsID_beds/'

# 1. significant QTLs from mash or rasqual
## 1.1 mash
celltypes <- readLines('~/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt')
df_lfsr <- readRDS('~/yuzhao1/work/atac_gca2024/26mash3/2mash/df_lfsr.rds')
qtls_ctList <- list()
for (ct in celltypes) {
  qtls_ctList[[ct]] <- rownames(df_lfsr)[df_lfsr[, ct] <= 0.05]
}
qtls_mash <- unique(unlist(qtls_ctList))
rsID_mash <- qtls_mash %>% strsplit(split = '@', fixed=T) %>% sapply(.,`[[`,2) %>% unique()

## 1.2 rasqual
filename_rds <- paste0('~/yuzhao1/work/atac_gca2024/24rasqual2/1rasqual_results/FDR_results_list.rds')
FDR_results_list <- readRDS(filename_rds)
qtls_rasqual_rsIDs <- lapply(FDR_results_list, function(x) paste0(x$rsID)) %>% unlist() %>% unique()

## 1.3 rasqual q<=0.1
qtls_q_0_1_list <- readRDS('/home/yuzhao1/yuzhao1/work/atac_gca2024/25rasqual3/1prepareQTLs/qtls_q_0_1_list.rds')
qtls_rasqual_q01_rsIDs <- lapply(qtls_q_0_1_list, function(x) x$rsID) %>% unlist() %>% unique()


## collect both
rsID_all <- unique(c(rsID_mash, qtls_rasqual_rsIDs, qtls_rasqual_q01_rsIDs))
df_bed <- data.frame(chr = strsplit(rsID_all, ':', fixed=T) %>% sapply(.,`[[`,1),
                     start = format(strsplit(rsID_all, ':', fixed=T) %>% sapply(.,`[[`,2) %>% as.numeric() - 1, scientific = FALSE),
                     end = format(strsplit(rsID_all, ':', fixed=T) %>% sapply(.,`[[`,2) %>% as.numeric(), scientific = FALSE),
                     ID = rsID_all)

# 2. split into 200 jobs
n_jobs <- 200 # number of chunks/jobs
ids_to_split <- 1:length(rsID_all)
ids_chunk <- split(ids_to_split, cut(seq_along(ids_to_split), n_jobs, labels = FALSE))

# 3. write snps to test into files
df_bed$chunkID <- 0
for (i in 1:n_jobs) {
  bedPath <- paste0(bedDir, i, '.bed')
  tmpIDs <- ids_chunk[[i]]
  df_bed[tmpIDs, 'chunkID'] <- i
  write.table(df_bed[tmpIDs,],
              bedPath,
              row.names = F,
              col.names = F,
              sep="\t",
              quote=FALSE)
}
saveRDS(df_bed, paste0(bedDir, '00rsID_chunkID_reference.rds'))









