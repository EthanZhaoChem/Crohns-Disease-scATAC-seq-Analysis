dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)

addArchRThreads(24)


######### healthy #########
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_healthy")
proj <- addMotifAnnotations(ArchRProj = proj,  motifSet = "cisbp", force = TRUE)
proj <- addBgdPeaks(proj, force = TRUE)
saveArchRProject(ArchRProj = proj, load = F)

######### sample 100 #########
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_subSampling_max100")
proj <- addMotifAnnotations(ArchRProj = proj,  motifSet = "cisbp", force = TRUE)
proj <- addBgdPeaks(proj, force = TRUE)
saveArchRProject(ArchRProj = proj, load = F)


proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_subSampling_max100_topic_epithelial/")
proj <- addMotifAnnotations(ArchRProj = proj,  motifSet = "cisbp", force = TRUE)
proj <- addBgdPeaks(proj, force = TRUE)
saveArchRProject(ArchRProj = proj, load = F)


proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_subSampling_max100_topic_immune/")
proj <- addMotifAnnotations(ArchRProj = proj,  motifSet = "cisbp", force = TRUE)
proj <- addBgdPeaks(proj, force = TRUE)
saveArchRProject(ArchRProj = proj, load = F)


proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_subSampling_max100_topic_stromal/")
proj <- addMotifAnnotations(ArchRProj = proj,  motifSet = "cisbp", force = TRUE)
proj <- addBgdPeaks(proj, force = TRUE)
saveArchRProject(ArchRProj = proj, load = F)



