dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)

source('~/yuzhao1/scripts/plot.R')
addArchRThreads(4)
addArchRLocking(locking = F)

out.dir <- '~/yuzhao1/work/atac_gca2024/16cCRE/rds/'

# read files
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1/")

# ############# compute ############
# proj <- addCoAccessibility(
#   ArchRProj = proj,
#   reducedDims = "IterativeLSI"
# )
# 
# proj <- addPeak2GeneLinks(
#   ArchRProj = proj,
#   reducedDims = "IterativeLSI"
# )
# saveArchRProject(proj, load = T)

############# 1, get link annotation ############
# enhancer or promoter region
archr_conns <- getPeak2GeneLinks(
  ArchRProj = proj,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = FALSE
)

peakDic <- paste(archr_conns@metadata$peakSet@seqnames,
                 start(archr_conns@metadata$peakSet@ranges),
                 end(archr_conns@metadata$peakSet@ranges),
                 sep="_")

geneDic <- archr_conns@metadata$geneSet$name

archr_conns_updated <- data.frame(matrix(nrow = nrow(archr_conns),
                                         ncol = 3))
archr_conns_updated$X1 <- unlist(mapvalues(archr_conns$idxATAC, 1:length(peakDic), to = peakDic, warn_missing = F))
archr_conns_updated$X2 <- unlist(mapvalues(archr_conns$idxRNA, 1:length(geneDic), to = geneDic, warn_missing = F))
archr_conns_updated$X3 <- archr_conns$Correlation
colnames(archr_conns_updated) <- c('peak', 'gene', 'corr')

cCREs <- unique(archr_conns_updated$peak)
length(cCREs)


############# 2, get peak annotation ############
peakset_unionDic <- data.frame(chr = proj@peakSet@seqnames,
                               ranges = proj@peakSet@ranges,
                               peakType = proj@peakSet$peakType,
                               nearestGene = proj@peakSet$nearestGene)
rownames(peakset_unionDic) <- paste0(peakset_unionDic$chr, '_', peakset_unionDic$ranges.start, '_', peakset_unionDic$ranges.end)

#save for later use
saveRDS(peakset_unionDic, paste0(out.dir, 'allPeaks.rds'))



# filter for cCRE

archr_conns_final <- cbind(peakset_unionDic[archr_conns_updated$peak, ], archr_conns_updated)
colnames(archr_conns_final) <- gsub('ranges.', '', colnames(archr_conns_final))
colnames(archr_conns_final)[colnames(archr_conns_final)=='gene'] <- 'linkedGene'
# 34677/89525 peaks' linked gene and nearest gene are



############# 3, get tf annotation ################
matches <- readRDS(proj@peakAnnotation$Motif$Matches)
matches.mtx <- matches@assays@data$matches
rownames(matches.mtx) <- paste(matches@rowRanges@seqnames,
                               start(matches@rowRanges@ranges),
                               end(matches@rowRanges@ranges),
                               sep="_")

peak_tf_mtx <- matches.mtx[cCREs,]
colnames(peak_tf_mtx) %<>%
  strsplit(., split = '_', fixed=T) %>%
  sapply(.,`[[`,1)

peak_tf_df <- as.data.frame(peak_tf_mtx)
peak_tf_df <- peak_tf_df %>%
  mutate(across(everything(), ~ ifelse(.x == TRUE, 1.0, 0.0)))

############# 4, combine everything ################

cCREs_df <- cbind(archr_conns_final, peak_tf_df[archr_conns_final$peak, ])
saveRDS(cCREs_df, paste0(out.dir, 'cCREs_enhancer&promoter.rds'))












