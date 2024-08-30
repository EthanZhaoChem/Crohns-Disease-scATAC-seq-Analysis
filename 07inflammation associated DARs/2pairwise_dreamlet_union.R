# notes: in pair 4(adj vs non), we did not include age (explained as random effect in patient factor) 
# (need to modify two starred lines of code for pair4)


# dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')

args <- commandArgs(trailingOnly = TRUE)
testpairID <- as.integer(args[1]) # ranging from 1 to 6

library(stringr)
library(ArchR)
library(SingleCellExperiment)
library(dreamlet)
library(muscat)
library(ExperimentHub)
library(zenith)
library(scater)
library(Matrix)
source('~/yuzhao1/scripts/helper_archr.R')
addArchRThreads(1)
out.dir <- '~/yuzhao1/work/atac_gca2024/7dreamlet/differential_test/'

# read files
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1/")
lineage <- 'union'

############################################################
# add patient metadata ###########
df1 <- read.csv('~/yuzhao1/work/atac_gca2024/0metadata/meta_Ethan_curated_20240311.csv')
proj$age <- mapvalues(proj$Sample, from = df1$sample, df1$Age) %>% unlist() %>% as.numeric()
proj$sex <- mapvalues(proj$Sample, from = df1$sample, df1$Sex) %>% unlist()


############################################################
# create original single cell experimen object
peakMtx <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/', 'peakMat_unbinarized.rds'))
peakMtx$age <- mapvalues(peakMtx$Sample %>% as.character(), from = df1$sample, df1$Age) %>% unlist() %>% as.numeric()
peakMtx$sex <- mapvalues(peakMtx$Sample %>% as.character(), from = df1$sample, df1$Sex) %>% unlist()


############################################################
# pairwise combinations in TI and AC separately
pairs <- data.frame(contrast1 = c("Control", "Control", "Control", "nonInf", "nonInf", "adjInf"),
                    contrast2 = c("nonInf", "adjInf", "inf", "adjInf", "inf", "inf"))


############################################################
conditions_test <- pairs[testpairID, c('contrast1', 'contrast2')]
cat('testing: ', ' ', conditions_test[[1]], ' ', conditions_test[[2]], '\n')

# subset cells
subset_cells <- rownames(colData(peakMtx))[colData(peakMtx)$inflammation_status %in% conditions_test]
peakMtx_sub <- peakMtx[, colnames(peakMtx) %in% subset_cells]

## clean env
rm(peakMtx)
rm(proj)
##

se <- SummarizedExperiment(assays=list(counts=peakMtx_sub@assays@data$PeakMatrix),
                           rowRanges=peakMtx_sub@rowRanges,
                           colData=peakMtx_sub@colData)
sce <- as(se, "SingleCellExperiment")
sce$anno1 <- factor(sce$anno1, levels = unique(sce$anno1))
sce$Sample <- factor(sce$Sample, levels = unique(sce$Sample))
rm(se)

# create pseudo-bulk
pb <- aggregateToPseudoBulk(sce,
                            assay = "counts",
                            cluster_id = "anno1",
                            sample_id = "Sample",
                            verbose = T)

# Normalize and apply voom/voomWithDreamWeights
# the resulting object of class dreamletProcessedData stores
# normalized data and other information


# starred:  
# res.proc = processAssays(pb, ~ inflammation_status + (1|Patient_ID_masked), # this is for pair 4 only, adj vs non 
res.proc = processAssays(pb, ~ inflammation_status + sex + age + (1|Patient_ID_masked),
                         min.cells=3, # remove samples a fewer than N cells
                         min.count=1, # minimum number of reads for a gene to be considered expressed in a sample
                         min.samples = 6, # remove cell types shown in less than 6 samples
                         min.prop = 0, # minimum proportion of retained samples with non-zero counts for a peak to be retained
                         isCounts = T,
                         normalize.method = "TMM",
                         useCountsWeights = TRUE,
                         quiet = FALSE)
cat('\n------------------------ end of part 1 ------------------------- \n')

# evaluated on the voom normalized data
# starred:  
# res.dl = dreamlet(res.proc, ~ inflammation_status + (1|Patient_ID_masked)) # this is for pair 4 only, adj vs non 
res.dl = dreamlet(res.proc, ~ inflammation_status + sex + age + (1|Patient_ID_masked))

####################################################################
saveRDS(res.dl, paste0(out.dir, 'dl_', lineage, '_', conditions_test[[1]], '_', conditions_test[[2]]))
cat('------------------------ end of script ------------------------- \n')


