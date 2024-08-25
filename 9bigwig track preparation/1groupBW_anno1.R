library(stringr)
library(ArchR)
addArchRThreads(1)

############################ 1. union ###############################
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1")
getGroupBW(
  ArchRProj = proj,
  groupBy = "anno1",
  normMethod = "ReadsInTSS",
  tileSize = 25,
  maxCells = 2000,
  ceiling = 4,
  verbose = TRUE, 
  threads = 1
)

