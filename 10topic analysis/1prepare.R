library(fastTopics)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(plyr)
library(dplyr)
library(stringr)
library(ArchR)

source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/scripts/helper_archr.R')

plot.dir <- '~/yuzhao1/work/atac_gca2024/13fasttopic/plots/'
data.dir <- '~/yuzhao1/work/atac_gca2024/13fasttopic/rds/'


# ################################## union ##############################################
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_subSampling_max100/")
union_metadata <- as.data.frame(proj@cellColData)

peak.mtx <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/peakMat_binarized.rds'))
counts <- peak.mtx@assays@data$PeakMatrix
counts <- t(counts)
colnames(counts) <- as.character(peak.mtx@rowRanges, ignore.strand=FALSE) %>% gsub(':', '_', .) %>% gsub('-', '_', .)

save(list = c("counts"), 
     file = file.path(data.dir, "union_raw_sub100.RData"))


# ################################## epithelial ##############################################
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_subSampling_max100_topic_epithelial/")
epithelial_metadata <- as.data.frame(proj@cellColData)

peak.mtx <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/peakMat_binarized.rds'))
counts <- peak.mtx@assays@data$PeakMatrix
counts <- t(counts)
colnames(counts) <- as.character(peak.mtx@rowRanges, ignore.strand=FALSE) %>% gsub(':', '_', .) %>% gsub('-', '_', .)

save(list = c("counts"), 
     file = file.path(data.dir, "epithelial_raw_sub100.RData"))


# ################################## immune ##############################################
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_subSampling_max100_topic_immune/")
immune_metadata <- as.data.frame(proj@cellColData)

peak.mtx <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/peakMat_binarized.rds'))
counts <- peak.mtx@assays@data$PeakMatrix
counts <- t(counts)
colnames(counts) <- as.character(peak.mtx@rowRanges, ignore.strand=FALSE) %>% gsub(':', '_', .) %>% gsub('-', '_', .)

save(list = c("counts"), 
     file = file.path(data.dir, "immune_raw_sub100.RData"))


# ################################## stromal##############################################
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_subSampling_max100_topic_stromal/")
stromal_metadata <- as.data.frame(proj@cellColData)

peak.mtx <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/peakMat_binarized.rds'))
counts <- peak.mtx@assays@data$PeakMatrix
counts <- t(counts)
colnames(counts) <- as.character(peak.mtx@rowRanges, ignore.strand=FALSE) %>% gsub(':', '_', .) %>% gsub('-', '_', .)

save(list = c("counts"), 
     file = file.path(data.dir, "stromal_raw_sub100.RData"))

