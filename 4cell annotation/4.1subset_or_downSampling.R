dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/work/final_GCArna/scripts/gca_colors.R')
source('~/yuzhao1/work/final_GCArna/scripts/gca_markers.R')

addArchRThreads(24)


######### healthy #########
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1/")
cells_healthy <- proj$cellNames[proj$disease_status == 'Control']
proj.healthy <- subsetArchRProject(
  ArchRProj = proj,
  cells = cells_healthy,
  outputDirectory = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_healthy",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = T
)

######### sample 100 #########
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1/")
set.seed(6)
cells_metadata <- as.data.frame(proj@cellColData)
df <- cells_metadata[, c('Sample', 'anno1')]


## max 100
max_cells_perSample_perCluster <- 100

cells_sampled_union <- c()
for (temp.sample in unique(df$Sample)){
  df1 <- df[df$Sample==temp.sample,]
  for(temp.anno1 in unique(df1$anno1)){
    df2 <- df1[df1$anno1==temp.anno1,]
    nCells_temp <- min(max_cells_perSample_perCluster, nrow(df2))
    cells_sampled_temp <- sample(rownames(df2), nCells_temp)
    cells_sampled_union <- c(cells_sampled_union, cells_sampled_temp)
  }
}
length(cells_sampled_union)

proj.union.subSample100 <- subsetArchRProject(
  ArchRProj = proj,
  cells = cells_sampled_union,
  outputDirectory = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_subSampling_max100",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = T
)

######### for each lineage, sample 1-- #########
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_subSampling_max100/")

cells.stromal <- proj$cellNames[proj$category1 == 'stromal']
proj.stromal <- subsetArchRProject(
  ArchRProj = proj,
  cells = cells.stromal,
  outputDirectory = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_subSampling_max100_topic_stromal/",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = T
)

cells.epithelial <- proj$cellNames[proj$category1 == 'epithelial']
proj.epithelial <- subsetArchRProject(
  ArchRProj = proj,
  cells = cells.epithelial,
  outputDirectory = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_subSampling_max100_topic_epithelial/",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = T
)

cells.immune <- proj$cellNames[proj$category1 == 'immune']
proj.immune <- subsetArchRProject(
  ArchRProj = proj,
  cells = cells.immune,
  outputDirectory = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_subSampling_max100_topic_immune/",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = T
)




