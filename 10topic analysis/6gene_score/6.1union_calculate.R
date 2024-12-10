dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
args <- commandArgs(trailingOnly = T)
topic_id_pathway  <- as.numeric(args[1])

library(fastTopics)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(plyr)
library(dplyr)
library(stringr)
library(ArchR)
library(Seurat)
library(fitdistrplus)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/scripts/deg_dep_utils.R')
source('~/yuzhao1/work/final_GCArna/scripts/gca_colors.R')
source('~/yuzhao1/scripts/helper_archr.R')
source('~/yuzhao1/work/atac_gca2024/13fasttopic/6gene_score/helper_gene_annotation.R')
source('~/yuzhao1/work/atac_gca2024/13fasttopic/6gene_score/helper_gene_scores.R')
source('~/yuzhao1/work/atac_gca2024/13fasttopic/6gene_score/helper_gene_gsea.R')

dir.pathway <- '~/yuzhao1/work/atac_gca2024/13fasttopic/6gene_score/pathways/'
############################# gene reference ###############################
# proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1/")
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# OrgDb <- org.Hs.eg.db
# keytypes(org.Hs.eg.db)
# genes_all <- getGeneAnnotation(proj)
# genes_all <- genes_all$genes
# genes_all <- genes_all[!is.na(genes_all$symbol)]
# genes_annoType <- get_gene_annotations(txdb, OrgDb, columns_extract = c("SYMBOL", "ENSEMBL", "GENETYPE"), single.strand.genes.only = F)
# genes_all$GENETYPE <- mapvalues(genes_all$symbol, genes_annoType$SYMBOL, genes_annoType$GENETYPE, warn_missing = F)
# 
# genes_df <- data.frame(chr = seqnames(genes_all),
#                        start = start(genes_all),
#                        end = end(genes_all), 
#                        strand = strand(genes_all),
#                        GeneID = genes_all$symbol,
#                        GENETYPE = genes_all$GENETYPE)
# genes_df <- genes_df[genes_df$GENETYPE=='protein-coding', ]
# 
# write.csv(genes_df, '/project/gca/yuzhao1/work/atac_gca2024/13fasttopic/6gene_score/protein_coding_genes.csv')
genes_df <- read.csv('/project/gca/yuzhao1/work/atac_gca2024/13fasttopic/6gene_score/protein_coding_genes.csv', row.names = 1)

############################# gene score ######################################
nTopics <- 45
parallel_results_dir <- '~/yuzhao1/work/atac_gca2024/13fasttopic/rds/fit_union_sub100_k45_converged_de_vsnull/'
nThreads <- 20
parallel_results <- list()
z_all <- matrix(0, nrow = 0, ncol = nTopics) # z-scores for posterior mean LFC estimates.
for (i in 1:nThreads) {
  parallel_results[[i]] <- readRDS(paste0(parallel_results_dir, i, '.rds'))
  z_all <- rbind(z_all, parallel_results[[i]]$z)
}

ATAC.regions <- data.frame(chr = rownames(z_all) %>% strsplit(split = '_', fixed=T) %>% sapply(.,`[[`,1),
                           start = rownames(z_all) %>% strsplit(split = '_', fixed=T) %>% sapply(.,`[[`,2) ,
                           end = rownames(z_all) %>% strsplit(split = '_', fixed=T) %>% sapply(.,`[[`,3) )

df_gene_scores <- compute_gene_scores_tss_model(
  Z = z_all,
  ATAC.regions = ATAC.regions,
  genes = genes_df,
  transform = "none", # c("abs", "squared", "none")
  normalization = "l2", # c( "sum", "l2", "none")
)


df_topGenes <- data.frame(matrix(0, nrow = 250, ncol = nTopics))
colnames(df_topGenes) <- paste0('k', 1:nTopics)
for(topic_id in 1:nTopics){
  df_topGenes[, topic_id] <- rownames(df_gene_scores)[order(df_gene_scores[[paste0('k', topic_id)]], decreasing = T)[1:250]]
}
  
# write.csv(df_topGenes, '~/yuzhao1/work/atac_gca2024/13fasttopic/6gene_score/genes/top250_proteinCoding.csv')

# for supplementary table
df_topGenes <- data.frame(matrix(0, nrow = 1000, ncol = nTopics * 2))
for(topic_id in 1:nTopics){
  colnames(df_topGenes)[[2*topic_id-1]] <- paste0('Topic', topic_id, '_Gene')
  colnames(df_topGenes)[[2*topic_id]] <- paste0('Topic', topic_id, '_GeneScore')
  xx <- df_gene_scores[, paste0('k', topic_id)]
  df_topGenes[, 2*topic_id-1] <- rownames(df_gene_scores)[order(xx, decreasing = T)[1:1000]]
  df_topGenes[, 2*topic_id] <- xx[order(xx, decreasing = T)[1:1000]]
}
write.csv(df_topGenes, '~/yuzhao1/work/atac_gca2024/13fasttopic/6gene_score/genes/top1k_proteinCoding.csv')


############################### pathway ###################################
cat(paste0('working on ', topic_id_pathway, '\n'))
gvec <- df_gene_scores[[paste0('k', topic_id_pathway)]]
ddf <- data.frame(gene=rownames(df_gene_scores),
                  z_logFC=gvec)
gvec <- GetGenelist(ddf, by = 'z_logFC')
reslist <- myEnrichGSEA_msigdb(gvec)
resdf <- bind_rows(lapply(reslist, function(x) x@result))
saveRDS(resdf, paste0(dir.pathway, 'separated/k', topic_id_pathway, '.rds'))

library(openxlsx)
write.csv(resdf, paste0(dir.pathway, 'separated/k', topic_id_pathway, '.csv'))


########### summarize and save ##########
list1 <- list()
for (i in 1:nTopics){
  xx <- readRDS(paste0(dir.pathway, 'separated/k', i, '.rds'))
  list1[[paste0('k', i)]] <- xx
}

library(openxlsx)
saveRDS(list1, paste0(dir.pathway, 'all_topics.rds'))
write.xlsx(list1, paste0(dir.pathway, 'all_topics.xlsx'))






                            
                                