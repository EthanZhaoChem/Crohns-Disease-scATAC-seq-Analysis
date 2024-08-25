# this step is not necessary if the peak calling was done without ceiling

dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)
addArchRThreads(24)

# # read files
proj_union <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1/")
proj_healthy <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_healthy/")
proj_sub100 <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_subSampling_max100/")
proj_epithelial <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_subSampling_max100_topic_epithelial/")
proj_immune <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_subSampling_max100_topic_immune/")
proj_stromal <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_subSampling_max100_topic_stromal/")


temp_lineage <- 'union'
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)
proj <- addPeakMatrix(proj, force = T, ceiling=10^9)
saveArchRProject(proj)

temp_lineage <- 'healthy'
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)
proj <- addPeakMatrix(proj, force = T, ceiling=10^9)
saveArchRProject(proj)

temp_lineage <- 'sub100'
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)
proj <- addPeakMatrix(proj, force = T, ceiling=10^9)
saveArchRProject(proj)

temp_lineage <- 'epithelial'
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)
proj <- addPeakMatrix(proj, force = T, ceiling=10^9)
saveArchRProject(proj)

temp_lineage <- 'immune'
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)
proj <- addPeakMatrix(proj, force = T, ceiling=10^9)
saveArchRProject(proj)

temp_lineage <- 'stromal'
proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)
proj <- addPeakMatrix(proj, force = T, ceiling=10^9)
saveArchRProject(proj)




















