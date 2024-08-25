dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)

source('~/yuzhao1/scripts/plot.R')
addArchRLocking(locking = F)

out.dir <- '~/yuzhao1/work/atac_gca2024/4peaks/plots/'
union_peakset_filePath <- '~/yuzhao1/work/atac_gca2024/4peaks/DARs/union_peakset.bed'
dir.create(out.dir, showWarnings = F, recursive = T)

# # read files
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1/")

############################ 4. union ###############################
# addArchRThreads(48)
# # add pseudo group replicates
# proj <- addGroupCoverages(ArchRProj = proj,
#                           groupBy = "anno1",
#                           useLabels = T,
#                           minReplicates = 6,
#                           maxReplicates = 100,
#                           minCells = 50,
#                           maxCells = 10000,
#                           force = T)
# ## save
# saveArchRProject(ArchRProj = proj, load = T)

addArchRThreads(2)
pathToMacs2 <- '/home/yuzhao1/.local/bin/macs2'
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "anno1",
  pathToMacs2 = pathToMacs2,
  reproducibility = "3",
  peaksPerCell = 500,
  maxPeaks = 150000,
  minCells = 50,
  excludeChr = c("chrM", "chrY"),
  extsize = 150,
  cutOff = 0.05,
  extendSummits = 250,
  force = T
)
## save
saveArchRProject(ArchRProj = proj, load = T)

## add peak matrix to object
addArchRThreads(1)
proj <- addPeakMatrix(proj, force = T, ceiling=10^9)

## save
saveArchRProject(ArchRProj = proj, load = T)

# save union peaks as bed
peaks_chr <- proj@peakSet@seqnames
peaks_start <- start(proj@peakSet@ranges)
peaks_end <- end(proj@peakSet@ranges)
peaks_filtered_all_bed <- data.frame(chr = peaks_chr,
                                     start = peaks_start,
                                     end = peaks_end)
write.table(peaks_filtered_all_bed,
            union_peakset_filePath, 
            row.names = F, col.names = F, sep="\t", quote=FALSE)


