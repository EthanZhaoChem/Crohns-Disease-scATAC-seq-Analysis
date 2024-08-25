library(stringr)
library(ArchR)
library(ComplexHeatmap)
metadata <- read.table('~/yuzhao1/work/atac_gca2024/0metadata/meta_Ethan_curated_20240211.csv', header = T, sep = ',')
sample_pairs_all <- readLines('~/yuzhao1/work/final_GCAatac/18glimpse/4.0prepare_samplePair_all_chr1.txt')



sample_names <- unique(metadata$sample)
df <- data.frame(matrix(1, nrow = length(sample_names), ncol = length(sample_names)))
rownames(df) <- sample_names
colnames(df) <- sample_names

for (sample_pair in sample_pairs_all) {
  sample1 <- strsplit(sample_pair, ' ')[[1]][[1]]
  sample2 <- strsplit(sample_pair, ' ')[[1]][[2]]
  resultPath <- paste0("~/yuzhao1/work/final_GCAatac/18glimpse/4.3sample_concordance_all_chr1/", sample1, '_', sample2, ".rds")
  result <- readRDS(resultPath)
  df[sample1, sample2] <- result$nonreference_concordance
  df[sample2, sample1] <- result$nonreference_concordance
}

### plot ###
library(circlize)
color = c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac') %>% rev()
breaks <- seq(0, 1, length.out = length(color))
col_fun <- circlize::colorRamp2(breaks, color)

column_ha <- HeatmapAnnotation(foo = anno_text(sample_names, location = 1, 
                                               rot = 90, 
                                               gp = gpar(fontsize = 5,
                                                         fontface='plain')))

pdf('~/yuzhao1/work/final_GCAatac/18glimpse/4.3sample_concordance_all_chr1/00heatmap.pdf', height = 30, width = 40)
ht <- Heatmap(df, 
              name = paste0('NRC'),
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



pdf('~/yuzhao1/work/final_GCAatac/18glimpse/4.3sample_concordance_all_chr1/00heatmap_clean.pdf', height = 7, width = 8)
ht <- Heatmap(df, 
              name = paste0('NRC'),
              bottom_annotation = column_ha,
              col = col_fun,
              column_names_gp = grid::gpar(fontsize = 5),
              row_names_gp = grid::gpar(fontsize = 5),
              rect_gp = gpar(col = "black", lwd = 0.5),
              cluster_columns = F, cluster_rows = F,
              show_row_dend = F, show_column_dend = F, 
              show_row_names = T, show_column_names = F,
              row_names_side = "left",
              border = T,
              use_raster = F,
              show_heatmap_legend = T)
draw(ht, padding = unit(c(.1, .1, .1, .1), "npc"))
dev.off()







