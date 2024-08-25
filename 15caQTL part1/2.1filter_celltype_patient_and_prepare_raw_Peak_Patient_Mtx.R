dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/work/final_GCArna/scripts/gca_colors.R')
source('~/yuzhao1/work/final_GCArna/scripts/gca_markers.R')

out.dir <- '~/yuzhao1/work/atac_gca2024/19rasqual/2.1mtx_peak_patient_raw/'

addArchRThreads(1)
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1/")
peak.mtx <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/peakMat_unbinarized.rds'))
col.cellnames <- colnames(peak.mtx@assays@data$PeakMatrix)
row.peaknames <- paste0(seqnames(peak.mtx@rowRanges), '_', start(peak.mtx@rowRanges), '_', end(peak.mtx@rowRanges))
peak.mtx <- peak.mtx@assays@data$PeakMatrix
rownames(peak.mtx) <- row.peaknames


## select cell types to analyze
# check which cell types have more than 20 patients, with more than 20 cells per patient per cell type
cells_metadata <- as.data.frame(proj@cellColData)
df <- cells_metadata[, c('Patient_ID_masked', 'anno1')]
cellMtx_patientByCelltype <- data.frame(unclass(table(df)))
# apply cell number cutoff, colSums is nPatients with cells more than this cutoff
nPatients_perCelltype <- colSums(cellMtx_patientByCelltype > 20)
# apply patient number cutoff, get cell types to analyze
celltypes_enoughPatients <- names(nPatients_perCelltype)[nPatients_perCelltype > 20]
# write to main folder
writeLines(celltypes_enoughPatients, '~/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt')
for (temp.celltype in celltypes_enoughPatients){
  # nCells per patient
  temp.cellMtx_patient <- cellMtx_patientByCelltype[, temp.celltype]
  # patients who passed the cutoff
  temp.patients <- rownames(cellMtx_patientByCelltype)[temp.cellMtx_patient > 20]
  temp.patients <- mixedsort(temp.patients)
  writeLines(temp.patients, paste0('~/yuzhao1/work/atac_gca2024/19rasqual/0patientID_perCelltype/', temp.celltype, '.txt'))
}


## for each cell type, create a peak-patient matrix
list1 <- list()
for (temp.celltype in celltypes_enoughPatients){
  # nCells per patient
  temp.cellMtx_patient <- cellMtx_patientByCelltype[, temp.celltype]
  # patients who passed the cutoff
  temp.patients <- rownames(cellMtx_patientByCelltype)[temp.cellMtx_patient > 20]
  temp.patients <- mixedsort(temp.patients)
  
  # build peak-sample count matrix
  temp.peakMtx_peakByPatient <- matrix(0, ncol = length(temp.patients), nrow = nrow(peak.mtx))
  colnames(temp.peakMtx_peakByPatient) <- temp.patients
  rownames(temp.peakMtx_peakByPatient) <- rownames(peak.mtx)
  
  # cells from these patients and this celltype
  for (patient2 in temp.patients) {
    cells2 <- proj$cellNames[proj$Patient_ID_masked==patient2 & proj$anno1==temp.celltype]
    nPeaks_perCelltypePatient <- rowSums(peak.mtx[, cells2])
    temp.peakMtx_peakByPatient[, patient2] <- nPeaks_perCelltypePatient
  }
  list1[[temp.celltype]] <- temp.peakMtx_peakByPatient
}


for (temp.celltype in celltypes_enoughPatients){
  write.table(list1[[temp.celltype]], paste0(out.dir, temp.celltype, '.txt'), sep = '\t', quote = FALSE)
}



