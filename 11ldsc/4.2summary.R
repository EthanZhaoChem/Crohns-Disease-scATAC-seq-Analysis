fdr_threshold <- 0.1
nTopics <- 45 # set this for topic analysis, otherwise, ignore it
root_dir <- '/home/yuzhao1/yuzhao1/work/atac_gca2024/14ldsc/results/union_sub100_k45_daPeaks_positive_flexibleLpval30k_vsnull/'

# results_cell <- paste0(root_dir, '/output')
# post_process_cell <- paste0(root_dir, "/", 'post_process')

results_cell <- paste0(root_dir, '/output_conditioned')
post_process_cell <- paste0(root_dir, "/", 'post_process_conditioned')

# create plot folder
post_process_plot_cell <- paste0(post_process_cell, "/", 'plots')
dir.create(post_process_plot_cell, showWarnings = F)

library(ggplot2)
library(cowplot)
library(tidyverse)
library(plyr)
library(dplyr)
library(stringr)
library(rmeta)
library(data.table)
library(ArchR)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/work/final_GCArna/scripts/gca_colors.R')

########################## part 1, read files ##############################
sumstats_taskfile='~/yuzhao1/work/atac_gca2024/14ldsc/3sumstats_tasks'
traits <- readLines(sumstats_taskfile, warn = F)
annot_cell <- paste0(root_dir, '/annots')
annot_names <- list.dirs(annot_cell, recursive = F, full.names = F)
annot_names <- setdiff(annot_names, 'all_cellTypeSpecific')

# # order rownames of plot
# order_rownames <- setdiff(orders_anno1_group_highRes, 'all_cellTypeSpecific')
order_rownames <- paste0('daPeaks_Topic', 1:nTopics)
# order_rownames <- names(gca_colors_atac_union_anno1)


########################### part 2, initialization #############################
temp.mtx <- matrix(0, nrow = length(annot_names), ncol = length(traits))
rownames(temp.mtx) <- annot_names
colnames(temp.mtx) <- traits

df_tauStar_mlog10p <- temp.mtx
df_E_mlog10p <- temp.mtx
df_significant <- temp.mtx

# currently not useful 
df_tauStar <- temp.mtx
df_tauStar_p <- temp.mtx
df_E <- temp.mtx
df_tauStar_se <- temp.mtx
df_E_se <- temp.mtx
df_E_p <- temp.mtx

#################### part 3: identify pairs with positive tau star ############
pairs_positive_tau_star <- list()
pairs_names <- c()
p_tau_positive_tau_star <- c()
p_enrichment_positive_tau_star <- c()
count <- 1

for(trait in traits){
  df <- read.table(paste0(post_process_cell, '/', trait, "_ldsc_postprocess.txt"), sep = "\t")
  for(annot_name in annot_names){
    df_tauStar[annot_name, trait] <- df[annot_name, 'tau.star']
    df_tauStar_mlog10p[annot_name, trait] <- -log10(df[annot_name, 'p.tau.star'])
    df_E_mlog10p[annot_name, trait] <- -log10(df[annot_name, 'p.E.'])
    df_E[annot_name, trait] <- df[annot_name, 'E']
    
    if(df[annot_name, 'tau.star'] > 0){
      pairs_positive_tau_star[[count]] <- c(trait, annot_name)
      pairs_names <- c(pairs_names, paste0(trait, '-', annot_name))
      p_tau_positive_tau_star <- c(p_tau_positive_tau_star, df[annot_name, 'p.tau.star']) # this is p tau
      p_enrichment_positive_tau_star <- c(p_enrichment_positive_tau_star, df[annot_name, 'p.E.']) # this is p enr
      count <- count + 1
    }
  }
}
names(pairs_positive_tau_star) <- pairs_names



#################### part 4.1 plot: tau* ############
df <- df_tauStar
df <- df[order_rownames,]

library(circlize)
color = ArchRPalettes$solarExtra
breaks <- seq(-2, 2, length.out = length(color))
col_fun <- circlize::colorRamp2(breaks, color)

column_ha <- HeatmapAnnotation(foo = anno_text(traits, location = 1, rot = 30, 
                                               gp = gpar(fontsize = 10, fontface='bold')))

pdf(paste0(post_process_plot_cell, '/tau_star.pdf'), height = 13, width = 9)
ht <- Heatmap(df, 
              name = paste0('tau*'),
              bottom_annotation = column_ha,
              col = col_fun,
              rect_gp = gpar(col = "black", lwd = 2),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", df[i, j]), x, y, gp = gpar(fontsize = 10))
              },
              cluster_columns = F, cluster_rows = F,
              show_row_dend = F, show_column_dend = F, 
              show_row_names = T, show_column_names = F,
              row_names_side = "left",
              border = T,
              use_raster = F,
              show_heatmap_legend = T)
draw(ht, padding = unit(c(.1, .1, .1, .1), "npc"))
dev.off()



#################### part 4.2 plot: p tau* ############
df <- df_tauStar_mlog10p
df <- df[order_rownames,]

library(circlize)
color = c('white', 'red')
breaks <- seq(0, 10, length.out = length(color))
col_fun <- circlize::colorRamp2(breaks, color)

column_ha <- HeatmapAnnotation(foo = anno_text(traits, location = 1, rot = 30, 
                                               gp = gpar(fontsize = 10, fontface='bold')))

pdf(paste0(post_process_plot_cell, '/tau_star_mlog10p.pdf'), height = 13, width = 9)
ht <- Heatmap(df, 
              name = paste0('tauStar_mlog10p'),
              bottom_annotation = column_ha,
              col = col_fun,
              rect_gp = gpar(col = "black", lwd = 2),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", df[i, j]), x, y, gp = gpar(fontsize = 10))
              },
              cluster_columns = F, cluster_rows = F,
              show_row_dend = F, show_column_dend = F, 
              show_row_names = T, show_column_names = F,
              row_names_side = "left",
              border = T,
              use_raster = F,
              show_heatmap_legend = T)
draw(ht, padding = unit(c(.1, .1, .1, .1), "npc"))
dev.off()



#################### part 4.3 plot: E ############
df <- df_E
df <- df[order_rownames,]

library(circlize)
color = ArchRPalettes$solarExtra
breaks <- seq(-10, 10, length.out = length(color))
col_fun <- circlize::colorRamp2(breaks, color)

column_ha <- HeatmapAnnotation(foo = anno_text(traits, location = 1, rot = 30, 
                                               gp = gpar(fontsize = 10, fontface='bold')))

pdf(paste0(post_process_plot_cell, '/enrichment.pdf'), height = 13, width = 9)
ht <- Heatmap(df, 
              name = paste0('Enrichment'),
              bottom_annotation = column_ha,
              col = col_fun,
              rect_gp = gpar(col = "black", lwd = 2),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", df[i, j]), x, y, gp = gpar(fontsize = 10))
              },
              cluster_columns = F, cluster_rows = F,
              show_row_dend = F, show_column_dend = F, 
              show_row_names = T, show_column_names = F,
              row_names_side = "left",
              border = T,
              use_raster = F,
              show_heatmap_legend = T)
draw(ht, padding = unit(c(.1, .1, .1, .1), "npc"))
dev.off()


#################### part 4.4 plot: p E ############
df <- df_E_mlog10p
df <- df[order_rownames,]

library(circlize)
color = c('white', 'red')
breaks <- seq(0, 10, length.out = length(color))
col_fun <- circlize::colorRamp2(breaks, color)

column_ha <- HeatmapAnnotation(foo = anno_text(traits, location = 1, rot = 30, 
                                               gp = gpar(fontsize = 10, fontface='bold')))

pdf(paste0(post_process_plot_cell, '/enrichment_mlog10p.pdf'), height = 13, width = 9)
ht <- Heatmap(df, 
              name = paste0('mlog10p Enrichment'),
              bottom_annotation = column_ha,
              col = col_fun,
              rect_gp = gpar(col = "black", lwd = 2),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", df[i, j]), x, y, gp = gpar(fontsize = 10))
              },
              cluster_columns = F, cluster_rows = F,
              show_row_dend = F, show_column_dend = F, 
              show_row_names = T, show_column_names = F,
              row_names_side = "left",
              border = T,
              use_raster = F,
              show_heatmap_legend = T)
draw(ht, padding = unit(c(.1, .1, .1, .1), "npc"))
dev.off()





############################ part 5, significant p tau* ##############################
# # filter pairs by fdr or not
if(fdr_threshold == 'NA') {
  pairs_significant <- pairs_positive_tau_star
} else if(is.numeric(fdr_threshold)){
  fdr_positive_tau_star <- p.adjust(p_tau_positive_tau_star, method = 'fdr') # fdr is BH by default
  pairs_significant <- pairs_positive_tau_star[fdr_positive_tau_star < fdr_threshold]
  pairs_significant_fdrs <- fdr_positive_tau_star[fdr_positive_tau_star < fdr_threshold]
} else (stop("ERROR: fdr_threshold should be either a numer of a 'NA' character"))


# initialize mtx
df_tauStar_mlog10p <- temp.mtx
df_significant <- temp.mtx
df_significant_fdr <- temp.mtx + 1 # thus it is non significant by default

for(i in 1:length(pairs_significant)){
  tmp.pair <- pairs_significant[[i]]
  trait <- tmp.pair[[1]]
  annot_name <- tmp.pair[[2]]
  df <- read.table(paste0(post_process_cell, '/', trait, "_ldsc_postprocess.txt"), sep = "\t")
  
  df_tauStar_mlog10p[annot_name, trait] <- -log10(df[annot_name, 'p.tau.star'])
  df_significant[annot_name, trait] <- 1
  df_significant_fdr[annot_name, trait] <- pairs_significant_fdrs[[i]]
}

######################### part 5: plot #########################
# order rows
df <- df_tauStar_mlog10p
df <- df[order_rownames,]
df_significant_fdr <- df_significant_fdr[order_rownames,]

library(circlize)
color = c('white', 'red')
breaks <- seq(0, 10, length.out = length(color))
col_fun <- circlize::colorRamp2(breaks, color)

column_ha <- HeatmapAnnotation(foo = anno_text(traits, location = 1, rot = 30, 
                                               gp = gpar(fontsize = 10, fontface='bold')))

pdf(paste0(post_process_plot_cell, '/tau_star_mlog10p', '_FDR_', fdr_threshold, '.pdf'), height = 13, width = 9)
ht <- Heatmap(df, 
              name = paste0('-log10(tau* p-val), FDR < ', fdr_threshold),
              bottom_annotation = column_ha,
              col = col_fun,
              rect_gp = gpar(col = "black", lwd = 2),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", df[i, j]), x, y, gp = gpar(fontsize = 10))
              },
              cluster_columns = F, cluster_rows = F,
              show_row_dend = F, show_column_dend = F, 
              show_row_names = T, show_column_names = F,
              row_names_side = "left",
              border = T,
              use_raster = F,
              show_heatmap_legend = T)
draw(ht, padding = unit(c(.1, .1, .1, .1), "npc"))
dev.off()


pdf(paste0(post_process_plot_cell, '/tau_star_mlog10p', '_FDR_', fdr_threshold, '_colorOnly.pdf'), height = 13, width = 9)
ht <- Heatmap(df, 
              name = paste0('-log10(tau* p-val), FDR < ', fdr_threshold),
              bottom_annotation = column_ha,
              col = col_fun,
              rect_gp = gpar(col = "black", lwd = 0.5),
              cluster_columns = F, cluster_rows = F,
              show_row_dend = F, show_column_dend = F, 
              show_row_names = T, show_column_names = F,
              row_names_side = "left",
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(df_significant_fdr[i, j] < 0.001) {
                  grid.text("***", x, y)
                } 
                else if(df_significant_fdr[i, j] < 0.01) {
                  grid.text("**", x, y)
                }
                else if(df_significant_fdr[i, j] < 0.1) {
                  grid.text("*", x, y)
                }
                },
              border = T,
              use_raster = F,
              show_heatmap_legend = T)
draw(ht, padding = unit(c(.1, .1, .1, .1), "npc"))
dev.off()


            
########################### part 6: significant p E ##############################
# # filter pairs by fdr or not
if(fdr_threshold == 'NA') {
  pairs_significant <- pairs_positive_tau_star
} else if(is.numeric(fdr_threshold)){
  fdr_positive_tau_star <- p.adjust(p_enrichment_positive_tau_star, method = 'fdr') # fdr is BH by default
  pairs_significant <- pairs_positive_tau_star[fdr_positive_tau_star < fdr_threshold]
} else (stop("ERROR: fdr_threshold should be either a numer of a 'NA' character"))

# initialize mtx
df_E_mlog10p <- temp.mtx
df_significant <- temp.mtx

for(tmp.pair in pairs_significant){
  trait <- tmp.pair[[1]]
  annot_name <- tmp.pair[[2]]
  df <- read.table(paste0(post_process_cell, '/', trait, "_ldsc_postprocess.txt"), sep = "\t")
  
  df_E_mlog10p[annot_name, trait] <- -log10(df[annot_name, 'p.E.'])
  df_significant [annot_name, trait] <- 1
}


########################### part 6: plot ##############################
# order rows
df <- df_E_mlog10p
df <- df[order_rownames,]

## heatmap plot
library(circlize)
color = c('white', 'red')
breaks <- seq(0, 10, length.out = length(color))
col_fun <- circlize::colorRamp2(breaks, color)

column_ha <- HeatmapAnnotation(foo = anno_text(traits, location = 1, rot = 30, 
                                               gp = gpar(fontsize = 10, fontface='bold')))

pdf(paste0(post_process_plot_cell, '/enrichment_mlog10p', '_FDR_', fdr_threshold, '.pdf'), height = 13, width = 9)
ht <- Heatmap(df, 
              name = paste0('-log10(p-val), FDR < ', fdr_threshold),
              bottom_annotation = column_ha,
              col = col_fun,
              rect_gp = gpar(col = "black", lwd = 2),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", df[i, j]), x, y, gp = gpar(fontsize = 10))
              },
              cluster_columns = F, cluster_rows = F,
              show_row_dend = F, show_column_dend = F, 
              show_row_names = T, show_column_names = F,
              row_names_side = "left",
              border = T,
              use_raster = F,
              show_heatmap_legend = T)
draw(ht, padding = unit(c(.1, .1, .1, .1), "npc"))
dev.off()
                                