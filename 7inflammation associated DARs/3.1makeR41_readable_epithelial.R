# dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(SingleCellExperiment)
library(dreamlet)
library(muscat)
library(zenith)
library(ExperimentHub)
library(scater)
library(variancePartition)
source('~/yuzhao1/scripts/helper_archr.R')
addArchRThreads(1)

# ############################################################
# # provide the rowID reference table to peaks
# proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1")
# peakMtx <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/peakMat_unbinarized.rds'))
# 
# #this is different from peakset in archr
# saveRDS(peakMtx@rowRanges, paste0(proj@projectMetadata$outputDirectory, '/', 'peakMtx_unbinarized_rowranges.rds'))

############################################################
# union inf vs control
coef <- 'inflammation_statusinf'
filename <- '~/yuzhao1/work/atac_gca2024/7dreamlet/differential_test/dl_union_Control_inf'
res.dl <- readRDS(filename)
celltypes_all <- names(res.dl)

statistics_singleContrast <- list()
for (j in 1:length(celltypes_all)) {
  j.celltype <- celltypes_all[[j]]
  
  # if this coef is not tested after sample filtering
  if(!coef %in% colnames(res.dl[[j.celltype]]$coefficients)){
    next
  }
  j.statistics <- variancePartition::topTable(res.dl[[j.celltype]],
                                              coef = coef,
                                              number = 10000000)
  statistics_singleContrast[[j.celltype]] <- j.statistics
}

saveRDS(statistics_singleContrast, '~/yuzhao1/work/atac_gca2024/7dreamlet/differential_test/statistics_inf_vs_control.rds')

############################################################
# union inf vs control
coef <- 'inflammation_statusControl'
filename <- '~/yuzhao1/work/atac_gca2024/7dreamlet/differential_test/dl_union_Control_adjInf'
res.dl <- readRDS(filename)
celltypes_all <- names(res.dl)

statistics_singleContrast <- list()
for (j in 1:length(celltypes_all)) {
  j.celltype <- celltypes_all[[j]]
  
  # if this coef is not tested after sample filtering
  if(!coef %in% colnames(res.dl[[j.celltype]]$coefficients)){
    next
  }
  j.statistics <- variancePartition::topTable(res.dl[[j.celltype]],
                                              coef = coef,
                                              number = 10000000)
  statistics_singleContrast[[j.celltype]] <- j.statistics
}

saveRDS(statistics_singleContrast, '~/yuzhao1/work/atac_gca2024/7dreamlet/differential_test/statistics_control_vs_adjInf.rds')


############################################################
# union inf vs control
coef <- 'inflammation_statusnonInf'
filename <- '~/yuzhao1/work/atac_gca2024/7dreamlet/differential_test/dl_union_Control_nonInf'
res.dl <- readRDS(filename)
celltypes_all <- names(res.dl)

statistics_singleContrast <- list()
for (j in 1:length(celltypes_all)) {
  j.celltype <- celltypes_all[[j]]
  
  # if this coef is not tested after sample filtering
  if(!coef %in% colnames(res.dl[[j.celltype]]$coefficients)){
    next
  }
  j.statistics <- variancePartition::topTable(res.dl[[j.celltype]],
                                              coef = coef,
                                              number = 10000000)
  statistics_singleContrast[[j.celltype]] <- j.statistics
}

saveRDS(statistics_singleContrast, '~/yuzhao1/work/atac_gca2024/7dreamlet/differential_test/statistics_nonInf_vs_control.rds')




















