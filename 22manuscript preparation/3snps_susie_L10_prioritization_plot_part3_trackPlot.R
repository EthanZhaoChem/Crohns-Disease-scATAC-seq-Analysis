library(GenomicRanges) # required for all Gviz input
library(tidyverse) # manipulating tibbles
library(rtracklayer) # loading bigwigs/bed files
library(bigsnpr) # loading genotype data from 1000Genomes for LD calculation
library(Gviz) # main plotting
library(GenomicInteractions) # hic plots
library(data.table)
library(dplyr)
library(plyr)
library(mapgen)
library(stringr)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/work/final_GCArna/scripts/gca_colors.R')
source('~/yuzhao1/scripts/helper_archr.R')
# source('R/analysis_utils.R')
plot.dir <- '~/yuzhao1/work/atac_gca2024/0manu/plots/3trackplots/'

##################  ################ ##################  ################
# # get a rsID-hg38_pos correspondance
# GWAS <- as.data.frame(
#   fread(
#     paste0("~/yuzhao1/work/final_GCAatac/0gwas/b151_hg38/cd_40266_20161107_hg38_b151_gwas_annotated_notFilterFor1KG.txt")
#   )
# )
# GWAS <- GWAS[, c(1:9)]
# colnames(GWAS) <- c('snp', 'chr', 'pos', 'beta', 'se', 'pval', 'a0', 'a1', 'zscore')
# df_rsID_hg38pos <- GWAS[, c('snp', 'pos')] 
# df_rsID_hg38pos <- df_rsID_hg38pos[!duplicated(df_rsID_hg38pos$snp),]
# rownames(df_rsID_hg38pos) <- df_rsID_hg38pos$snp
# 
# gwas <- readRDS('~/yuzhao1/work/final_GCAatac/8susie/results/gwas_processed.rds')
# gwas$pos <- df_rsID_hg38pos[gwas$snp, 'pos']
# saveRDS(gwas, '~/yuzhao1/work/final_GCAatac/8susie/results/gwas_processed_hg38.rds')

##################  ################ ##################  ################
## read data
celltype_ideal_order <- names(gca_colors_atac_union_anno1)
palette <- gca_colors_atac_union_anno1
GWAS <- readRDS('~/yuzhao1/work/final_GCAatac/8susie/results/gwas_processed_hg38.rds') # this is the same gwas with hg19, but pos is hg38
GWAS <- GWAS[!is.na(GWAS$pos), ] # some pos are NA because they can't be mapped to hg38
big.snp <- snp_attach(rdsfile = '~/yuzhao1/resource/magma/g1000_eur.rds')
gene.annots <- readRDS('~/yuzhao1/work/0archive/alan_heart_project/genomic_annotations/hg38_gtf_genomic_annots.gr.rds')


finemapped_results <- read.csv('~/yuzhao1/work/final_GCAatac/8susie/results/cd_finemapping_unifprior_95loci_L10.csv', row.names = 1)
finemapped_results <- finemapped_results[, c('snp', 'susie_pip')]
GWAS_finemapped <- left_join(GWAS, finemapped_results)
GWAS_finemapped$susie_pip[is.na(GWAS_finemapped$susie_pip)] <- 0

##################  ################ ##################  ################
col_vec <- c("0-0.2"="#463699", "0.2-0.4"="#26bce1", "0.4-0.6"="#6efe68", "0.6-0.8"="#f8c32a", "0.8-1"="#db3d11")

## specify the region and snp to plot


## No.2
highlight.width <- 1000
plot.name <- 'LRRC32_rs11236797'
gr.start <-76580001
gr.end <-  76690000
snp.loc <- 76588605 # rs11236797
chr.id <- 'chr11'

## No.3
highlight.width <- 500
plot.name <- 'PLEKHG6_rs28999107'
gr.start <-6305001
gr.end <-  6395000
snp.loc <- 6383934 # rs28999107
chr.id <- 'chr12'

## No.4
highlight.width <- 500
plot.name <- 'GPR65_rs4462528'
gr.start <-87955001
gr.end <-  88013000
snp.loc <- 87959953 # rs4462528
chr.id <- 'chr14'

## No.5
highlight.width <- 500
plot.name <- 'SNX20_rs146528649'
gr.start <-50620001
gr.end <-  50690000
snp.loc <- 50627053 # rs146528649
chr.id <- 'chr16'

## No.6
highlight.width <- 2000
plot.name <- 'MIDN_rs4807569'
gr.start <-1117001
gr.end <-  1260000
snp.loc <- 1123379 # rs4807569
chr.id <- 'chr19'

## No.8
highlight.width <- 500
plot.name <- 'RSPO3_rs7741021'
gr.start <-127112001
gr.end <-  127157000
snp.loc <- 127147129 # rs7741021
chr.id <- 'chr6'

## No.9
highlight.width <- 1000
plot.name <- 'P4HA2_rs2188962'
gr.start <-132278001
gr.end <-  132448000
snp.loc <- 132435113 # rs2188962
chr.id <- 'chr5'

## No.10
highlight.width <- 100
plot.name <- 'SP110_rs13397985'
gr.start <-230223001
gr.end <-  230230000
snp.loc <- 230226508 # rs13397985 
chr.id <- 'chr2'

## No.11
highlight.width <- 50
plot.name <- 'HLA-DMB_rs111369841'
gr.start <-32950001
gr.end <-  32956000
snp.loc <- 32953327 # rs111369841
chr.id <- 'chr6'

## No.12
highlight.width <- 500
plot.name <- 'HLA-DRA_rs116980366'
gr.start <-32432001
gr.end <-  32476000
snp.loc <- 32473142 # rs116980366
chr.id <- 'chr6'

## No.13
highlight.width <- 100
plot.name <- 'PSMB8-AS1_rs74365910'
gr.start <-32837000
gr.end <-  32847000
snp.loc <- 32838776 # rs74365910
chr.id <- 'chr6'

## No.14
highlight.width <- 2000
plot.name <- 'HLA-DQA1_rs73739311'
gr.start <-32623001
gr.end <-  32818000
snp.loc <- 32796986 # rs73739311
chr.id <- 'chr6'

###### double snp
## No.1
highlight.width <- 500
plot.name <- 'IKZF1_rs1456896_rs9656588'
gr.start <- 50260001
gr.end <- 50310000
snp.loc1 <- 50264865 # rs1456896
snp.loc2 <- 50267184 # rs9656588
chr.id <- 'chr7'

## No.2
highlight.width <- 100
plot.name <- 'NCF4_rs4821544_rs760517'
gr.start <- 36859001
gr.end <- 36870000
snp.loc1 <- 36862461 # rs4821544
snp.loc2 <- 36862944 # rs760517
chr.id <- 'chr22'

## No.3
highlight.width <- 100
plot.name <- 'CIITA_rs7194862_rs6416647'
gr.start <- 10870001
gr.end <- 10880000
snp.loc1 <- 10871258
snp.loc2 <- 10871740
chr.id <- 'chr16'


##################  ################ ##################  ################
region.start <- gr.start
region.end <- gr.end
snp.chr <- chr.id
snp.p <- snp.loc # only used for highlight

## prepare input
ss <- region.start
ee <- region.end
pip.df <- GWAS_finemapped[paste0('chr', GWAS_finemapped$chr) == snp.chr,]
pip.df$pval <- -log10(pip.df$pval)
pip.df <- pip.df[pip.df$pos > ss & pip.df$pos < ee,]
curr.locus.gr <- GRanges(seqnames = snp.chr, 
                         IRanges(start = ss, 
                                 end = ee ))
# prepare P-value track
## Add LD information
pip.df$isAnnotate <- pip.df$susie_pip > 0.1
top.snp <- pip.df$bigSNP_index[which.max(pip.df$pval)]
top.snp.G <- big.snp$genotypes[, top.snp]
G.mat <- big.snp$genotypes[, pip.df$bigSNP_index]
r2.vals <- as.vector(cor(top.snp.G, G.mat))^2
r2.brackets <- cut(r2.vals, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                   labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1"))
pip.df$r2 <- r2.brackets


##################  ################ ##################  ################
## prepare tracks (pval track)
## actually make the track with Gviz
pval.df <- pip.df[,c("chr","pos", "pval","r2")] %>%
  mutate(start = pos, end = pos) %>% 
  dplyr::select(-pos) %>%
  pivot_wider(names_from = r2, values_from = "pval") 

pval.df.gr <- makeGRangesFromDataFrame(pval.df, keep.extra.columns = T)
seqlevelsStyle(pval.df.gr) <- "UCSC"
pval.track <- DataTrack(range = pval.df.gr,  
                        genome = "hg38", 
                        groups = factor(names(mcols(pval.df.gr)), levels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1")), 
                        col = col_vec, 
                        name = "-log10 pvalue")

# axis track - so we know where we are in the genome 
axisTrack <- GenomeAxisTrack()

#PIP TRACK
pip.track <- DataTrack(data = pip.df$susie_pip, chromosome = pip.df$chr, start = pip.df$pos, end = pip.df$pos, genome = "hg38", name = "PIP", col= "#463699")


# put all tracks into a list
list.of.tracks <- c(pval.track, pip.track,  axisTrack)

# # highlight a particular SNP
# ht <- HighlightTrack(trackList = list.of.tracks, 
#                      start = c(snp.p-highlight.width/2), width = highlight.width, 
#                      chromosome = as.character(seqnames(curr.locus.gr)), 
#                      col='white', fill = '#f5d3d3')

# highlight two particular SNPs
ht <- HighlightTrack(trackList = list.of.tracks, 
                      start = c(snp.loc1-highlight.width/2, snp.loc2-highlight.width/2), width = highlight.width, 
                      chromosome = as.character(seqnames(curr.locus.gr)), 
                      col='white', fill = '#f5d3d3')


# plot
pdf(paste0('~/yuzhao1/work/atac_gca2024/0manu/plots/3trackplots/', plot.name, '.pdf'), 
    width=12, 
    height=8)

plotTracks(ht, 
           chromosome = as.character(seqnames(curr.locus.gr)), 
           transcriptAnnotation = "symbol", 
           collapseTranscripts= 'longest', 
           from = start(curr.locus.gr), 
           to = end(curr.locus.gr),
           sizes = c(1, 0.6, 0.3),
           panel.only = F, 
           background.title = "white",
           col.title = "black",
           col.axis = "black")

dev.off()


