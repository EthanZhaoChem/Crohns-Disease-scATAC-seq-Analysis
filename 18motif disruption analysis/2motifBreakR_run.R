args <- commandArgs(trailingOnly = T)
bedChunkID <- args[1]

# supposed to calculate all snps
dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP150.GRCh38) 
library(BSgenome.Hsapiens.UCSC.hg38)     
library(BSgenome)
library(plyr)
library(stringr)

bedDir <- '~/yuzhao1/work/atac_gca2024/26motif_disruption3/1rsID_beds/'
out.dir <- '~/yuzhao1/work/atac_gca2024/26motif_disruption3/2results_allHumanTF/'

##################  ################ ##################  ################
bedPath <- paste0(bedDir, bedChunkID, '.bed')
snps.mb.frombed <- snps.from.file(file = bedPath,
                                  search.genome = BSgenome.Hsapiens.UCSC.hg38,
                                  format = "bed")

##################  ################ ##################  ################
# 2. prepare TF pwm subset to test
motifs_pool <- MotifDb::query(MotifDb, andStrings = c("hsapiens"))
motifs_pool <- motifs_pool[!is.na(mcols(motifs_pool)$geneSymbol)]

results <- motifbreakR(snpList = snps.mb.frombed,
                       verbose = T,
                       pwmList = motifs_pool, 
                       filterp = TRUE,
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25), 
                       BPPARAM = BiocParallel::MulticoreParam(workers = 30))

saveRDS(results, paste0(out.dir, bedChunkID, '.rds'))












