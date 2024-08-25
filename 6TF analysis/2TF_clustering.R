dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(ggplot2)
library(dplyr)
library(plyr)
library(stringr)
library(harmony)
library(Seurat)
library(ArchR)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggseqlogo)
library(seqLogo)
library(motifStack)
library(ade4)
library(TFBSTools)
source('~/yuzhao1/scripts/plot.R')

##################  ################ ##################  ################
### prepare pfm list
out.dir <- '~/yuzhao1/work/atac_gca2024/5TF/2TF_clustering/'

# reduce redundancy with cisbp logos
pwms <- chromVARmotifs::human_pwms_v2

# have checked all pseudo counts are zeros
pfms <- list()
for (i in 1:length(pwms)){
  pwm <- pwms[[i]]
  TF.name <- pwm@name
  pfm <- exp(pwm@profileMatrix)*0.25 
  pfms[[TF.name]] <- pfm
}

# pfm to list
motifs <- mapply(names(pfms), 
               pfms, 
               FUN=function(.ele, .pfm){new("pfm",mat=.pfm, name=.ele)},
               SIMPLIFY = FALSE)

##################  ################ ##################  ################
## prepare clustering
# hc <- clusterMotifs(motifs)
# saveRDS(hc, paste0(out.dir, 'hc.rds'))
hc <- readRDS(paste0(out.dir, 'hc.rds'))

num_clusters <- 300
clusters <- cutree(hc, k = num_clusters)
cluster_list <- vector("list", num_clusters)
for (i in 1:num_clusters) {
  cluster_list[[i]] <- hc$labels[which(clusters == i)]
}

# Write the data frame to a CSV file
df <- data.frame(Cluster = sapply(cluster_list, function(x) paste(x, collapse = ",")))
write.csv(df, "~/yuzhao1/work/atac_gca2024/5TF/2TF_clustering/300clusters.csv", row.names = FALSE, quote = FALSE)
saveRDS(clusters, '~/yuzhao1/work/atac_gca2024/5TF/2TF_clustering/300clusters.rds')

#################  ################ ##################  ################
# plot logo
motifs <- pfms
motifs <- mapply(names(motifs),
               motifs,
               FUN=function(.ele, .pfm){new("pfm",mat=.pfm, name=.ele)},
               SIMPLIFY = FALSE)

pdf(paste0(out.dir, 'singleTF.pdf'), width = 10, height = 3)
for (kk in motifs){
  print(motifStack::plot(kk))
}
dev.off()







