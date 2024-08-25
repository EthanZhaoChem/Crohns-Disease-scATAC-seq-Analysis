# read parameters
args <- commandArgs(trailingOnly = T)
file_counts <- as.character(args[1])
file_fasttopic <- as.character(args[2]) 
Dir_tmp <- as.character(args[3])
i  <- as.numeric(args[4])
nThreads <- as.numeric(args[5])
file_de <- paste0(Dir_tmp, '/', i, '.rds')

# preset
set.seed(6)
library(fastTopics)
library(tidyverse)
library(plyr)
library(dplyr)
library(stringr)
load(file_counts)
fit <- readRDS(file_fasttopic)
fit <- poisson2multinom(fit)

# locate the subset peaks to calculate
peakIDs <- 1:ncol(counts)
chunklength <- ceiling(length(peakIDs)/nThreads)
group.list <- split(peakIDs, ceiling(peakIDs / chunklength))
peaks_sub <- group.list[[i]]

# calculate and save
fit2   <- fit
fit2$F <- fit$F[peaks_sub, ]
fit2$Fn <- fit$Fn[peaks_sub, ]
fit2$Fy <- fit$Fy[peaks_sub, ]
counts2 <- counts[, peaks_sub]
# rm(counts, fit)

de <- de_analysis(fit2, 
                  counts2,
                  shrink.method = "none", 
                  lfc.stat = 'vsnull',
                  control = list(ns = 1e4, nc = 12))

saveRDS(de, file_de)






