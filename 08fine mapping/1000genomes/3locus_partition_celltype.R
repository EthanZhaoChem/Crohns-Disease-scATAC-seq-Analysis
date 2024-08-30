library(data.table)
library(dplyr)
library(stringr)
library(GenomicRanges)
library(bigsnpr)
library(mapgen)
library(susieR)
library(data.table)
library(Repitools)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/work/final_GCArna/scripts/gca_colors.R')
source('~/yuzhao1/scripts/helper_archr.R')


plot.dir <- '~/yuzhao1/work/final_GCAatac/8susie/3locus_partition_celltype/'

############ 1.1 read fine mapping results ############
# hg19 susie
merge_susie_sumstats_ethan <- function (susie_results, sumstats) 
{
  sumstats$susie_pip <- 0
  sumstats$CS <- '0'
  loci <- names(susie_results)
  for (l in loci) {
    n.snps <- length(susie_results[[l]]$pip)
    sumstats[sumstats$locus == as.numeric(l), "susie_pip"] <- susie_results[[l]]$pip
    snps.in.cs <- rep('0', n.snps)
    if (!is.null(susie_results[[l]]$sets$cs)) {
      csNames <- names(susie_results[[l]]$sets$cs)
      for (csName in csNames) {
        snps.in.cs[unlist(susie_results[[l]]$sets$cs[[csName]])] <- csName
      }
    }
    sumstats[sumstats$locus == as.numeric(l), "CS"] <- snps.in.cs
  }
  return(sumstats)
}

gwas.sumstats.sigloci <- readRDS('~/yuzhao1/work/final_GCAatac/8susie/results/gwas.sumstats.sigloci_95loci.rds')
susie_raw_L10 <- readRDS('~/yuzhao1/work/final_GCAatac/8susie/results/cd_finemapping_unifprior_95loci_L10.rds')
susie_merged_L10 <- merge_susie_sumstats_ethan(susie_results = susie_raw_L10, sumstats = gwas.sumstats.sigloci)
susie_clean_L10 <- susie_merged_L10[susie_merged_L10$CS!='0',]
susie_clean_L10 <- as.data.frame(susie_clean_L10)
susie_clean_L10 <- susie_clean_L10[susie_clean_L10$snp!='.', ]
susie_clean_L10$signal_ID <- paste0(susie_clean_L10$locus, '_', susie_clean_L10$CS)
rownames(susie_clean_L10) <- susie_clean_L10$snp
risk_loci_IDs <- unique(susie_clean_L10$locus)
df_RiskSNPs <- susie_clean_L10
df_RiskSNPs$pip <- df_RiskSNPs$susie_pip
df_RiskSNPs$chr <- paste0('chr', df_RiskSNPs$chr)

########################## 1.2 read published results (nature 2017) ####################
credibleSets <- read.table('~/yuzhao1/work/final_GCAatac/12snps/1nature2017cs/nature2017cs.csv', sep = ',', header = T)
snps_pip <- read.table('~/yuzhao1/work/final_GCAatac/12snps/1nature2017cs/snps.csv', sep = ',', header = T)
snps_pip <- snps_pip[!is.na(snps_pip$chr),] # remove NA rows based on chr

# there is an error in the snps_pip (rs3838334 should be in TNFSF8, remove the wrong one)
snps_pip <- snps_pip[!(snps_pip$variant == 'rs3838334' & snps_pip$Gene.refGene != 'TNFSF8'),]
rownames(snps_pip) <- snps_pip$variant

df_RiskSNPs <- snps_pip
df_RiskSNPs$pos <-  df_RiskSNPs$position
df_RiskSNPs$pip <- df_RiskSNPs$P_mean_95

############ 2. read locus information ############
# hg19 locus
locus_all <- readRDS('~/yuzhao1/work/final_GCAatac/8susie/LD_blocks_mapgen/LD_Blocks.rds') %>% as.data.frame(.)
colnames(locus_all) <- c('chr', 'start', 'end', 'locus_id')
locus_all$chr <- paste0('chr', locus_all$chr)
loci_risk <- locus_all[locus_all$locus_id %in% risk_loci_IDs, ]

############ 3. read peakset bed files ############
# hg19 non overlapping peak sets
beds_dir <- '~/yuzhao1/work/final_GCAatac/14ldsc/results/healthy_union_anno1Group_higherResolution_peaks/FDR0_1FC0_5/bed_lifted/'
# beds_dir <- '~/yuzhao1/work/final_GCAatac/14ldsc/results/fasttopic_TI_immune_sub100_k20_p0.95/bed_lifted/'
peaksets <- list()
peakset_names <- list.files(beds_dir, full.names = F)
for (xx in peakset_names) {
  xx.path <- paste0(beds_dir, xx)
  xx.df <- read.table(xx.path)
  colnames(xx.df) <- c('chr', 'start', 'end')
  peaksets[[xx]] <- xx.df
}

############ 4. assign a score to each cell type in each locus ############
score_celltype_locus <- data.frame(matrix(0, nrow = nrow(loci_risk), ncol = length(peakset_names)))
rownames(score_celltype_locus) <- paste0('locus_', loci_risk$locus_id)
colnames(score_celltype_locus) <- peakset_names

for (i in 1:nrow(df_RiskSNPs)) {
  snp_chr <- df_RiskSNPs[i, 'chr']
  snp_pos <- df_RiskSNPs[i, 'pos']
  snp_pip <- df_RiskSNPs[i, 'pip']
  gr_snp <- GRanges(seqnames = snp_chr, ranges = IRanges(start = snp_pos, end = snp_pos))
  cat(paste0('working on snp: ', i, ' out of ', nrow(df_RiskSNPs), '\n'))
  
  for (j in 1:nrow(score_celltype_locus)) {
    risk_locus_chr <- loci_risk[j, 'chr']
    risk_locus_start <- loci_risk[j, 'start']
    risk_locus_end <- loci_risk[j, 'end']
    gr_locus <- GRanges(seqnames = risk_locus_chr[[1]], ranges = IRanges(start = risk_locus_start[[1]], end = risk_locus_end[[1]]))

    if(any(findOverlaps(gr_snp, gr_locus)@from)){
      for (k in 1:ncol(score_celltype_locus)) {
        peaks_set_temp <- peaksets[[k]]
        gr_peakset <- GRanges(seqnames = peaks_set_temp$chr, ranges = IRanges(start = peaks_set_temp$start, end = peaks_set_temp$end))
        
        if(any(findOverlaps(gr_snp, gr_peakset)@from)){
          # update score
          cat(paste0('updating score: ', j, ' * ', k, '\n'))
          score_celltype_locus[j, k] <- score_celltype_locus[j, k] + snp_pip
        }
      }
    }
  }
}

saveRDS(score_celltype_locus, '~/yuzhao1/work/final_GCAatac/8susie/3locus_partition_celltype/score_celltype_locus_raw.rds')

############ normalize ############
# for locus with more than 1 credible set, should we divide the scores by N credible sets?
# or simply normalize the scores to a sum of 1 in one locus?, before doing that, make sure you removed all peaks bed column
df <- score_celltype_locus
df$all_cellTypeSpecific.bed <- NULL
df <- t(df[rowSums(df) > 0.25, ])
df_normalized <- apply(df, 2, function(x) x / sum(x))
xx <- orders_anno1_group_highRes[paste0(orders_anno1_group_highRes, '.bed')%in%rownames(df)]
df_normalized <- df_normalized[paste0(xx, '.bed'), ]


## heatmap plot
library(circlize)
col_fun = colorRamp2(c(0, 2), c("white", "#d6604d"))
col_fun = c("white", "#d6604d")
# column_ha <- HeatmapAnnotation(foo = anno_text(colnames(df_normalized), location = 1, rot = 45, 
#                                                gp = gpar(fontsize = 10, fontface='bold')))

pdf(paste0(plot.dir, 'locus_partition_celltype.pdf'), height = 10, width = 30)
ht <- Heatmap(df_normalized, 
              name = 'pip_proportion',
              # bottom_annotation = column_ha,
              col = col_fun,
              rect_gp = gpar(col = "black", lwd = 2),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", df_normalized[i, j]), x, y, gp = gpar(fontsize = 10))
              },
              cluster_columns = T, cluster_rows = F,
              show_row_dend = F, show_column_dend = F, 
              show_row_names = T, show_column_names = T,
              row_names_side = "left",
              border = T,
              use_raster = F,
              show_heatmap_legend = T)
draw(ht, padding = unit(c(.1, .1, .1, .1), "npc"))
dev.off()












