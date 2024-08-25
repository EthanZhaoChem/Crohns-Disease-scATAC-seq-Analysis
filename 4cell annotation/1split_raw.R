dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/work/final_GCArna/scripts/gca_colors.R')
source('~/yuzhao1/work/final_GCArna/scripts/gca_markers.R')

addArchRThreads(12)
# read files
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2/projdir/")
out.dir <- '~/yuzhao1/work/atac_gca2024/3annotation/plots/1split_raw/'
dir.create(out.dir, showWarnings = F, recursive = T)


# ############################# Section0: subset ################################
df_annotation_res0.5 <- list(
  'C1' ='epithelial',
  'C2' ='epithelial',
  'C3' ='epithelial',
  'C4' ='LQ',
  'C5' ='LQ',
  'C6' ='epithelial',
  'C7' ='epithelial',
  'C8' ='epithelial',
  'C9' ='epithelial',
  'C10' ='epithelial',
  'C11' ='epithelial',
  'C12' ='epithelial',
  'C13' ='epithelial',
  'C14' ='epithelial',
  'C15' ='epithelial',
  'C16' ='epithelial',
  'C17' ='epithelial',
  'C18' ='epithelial',
  'C19' ='immune',
  'C20' ='LQ',
  'C21' ='immune',
  'C22' ='LQ',
  'C23' ='immune',
  'C24' ='immune',
  'C25' ='immune',
  'C26' ='immune',
  'C27' ='immune',
  'C28' ='immune',
  'C29' ='immune',
  'C30' ='immune',
  'C31' ='immune',
  'C32' ='stromal',
  'C33' ='stromal',
  'C34' ='stromal',
  'C35' ='stromal',
  'C36' ='stromal',
  'C37' ='stromal',
  'C38' ='stromal'
)


proj$category1 <- unlist(mapvalues(proj$LSI_Clusters_res0.5, names(df_annotation_res0.5), df_annotation_res0.5))
saveArchRProject(ArchRProj = proj, load = T)


######################################################################
# epithelial
idxPass <- which(proj$category1 %in% c("epithelial"))
cellsPass.epithelial <- proj$cellNames[idxPass]
proj.epithelial <- subsetArchRProject(
  ArchRProj = proj,
  cells = cellsPass.epithelial,
  outputDirectory = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_epithelial1",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = T
)
proj.epithelial <- addImputeWeights(proj.epithelial)


######################################################################
# immune
idxPass <- which(proj$category1 %in% c("immune"))
cellsPass.immune <- proj$cellNames[idxPass]
proj.immune <- subsetArchRProject(
  ArchRProj = proj,
  cells = cellsPass.immune,
  outputDirectory = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_immune1",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = T
)
proj.immune <- addImputeWeights(proj.immune)

######################################################################
# stromal
idxPass <- which(proj$category1 %in% c("stromal"))
cellsPass.stromal <- proj$cellNames[idxPass]
proj.stromal <- subsetArchRProject(
  ArchRProj = proj,
  cells = cellsPass.stromal,
  outputDirectory = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_stromal1",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = T
)
proj.stromal <- addImputeWeights(proj.stromal)



























