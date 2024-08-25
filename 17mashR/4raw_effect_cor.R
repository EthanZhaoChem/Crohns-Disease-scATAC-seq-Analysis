dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(dplyr)
library(tidyverse)
library(data.table)
library(gtools)
library(limma)
library(ArchR)
library(ashr)
library(mashr)
library(ComplexHeatmap)
source('~/yuzhao1/work/atac_gca2024/scripts/gca_colors.R')
source('~/yuzhao1/work/final_GCArna/scripts/gca_markers.R')

"%&%" <- function(a, b) paste0(a, b)
##################  ################ ##################  ################
## presets
celltypes <- readLines('~/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt')
out.dir <- '~/yuzhao1/work/atac_gca2024/26mash3/4raw_effect_cor/'
prep.dir <- '~/yuzhao1/work/atac_gca2024/26mash3/1select_qtl/'

df_lfsr <- readRDS('~/yuzhao1/work/atac_gca2024/26mash3/2mash/df_lfsr.rds')
df_effect <- readRDS(paste0(prep.dir, 'df_effect_query.rds'))
df_effect$Feature_ID <- df_effect %>% rownames() %>% strsplit(split = '@', fixed=T) %>% sapply(.,`[[`,1) %>% unlist()
df_effect$rsID <- df_effect %>% rownames() %>% strsplit(split = '@', fixed=T) %>% sapply(.,`[[`,2) %>% unlist()

##################  ################ ##################  ################
qtls_ctList <- list()
for (ct in celltypes) {
  qtls_ctList[[ct]] <- rownames(df_lfsr)[df_lfsr[, ct] <= 0.05]
}

##################  ################ ##################  ################
df <- data.frame(matrix(0, nrow = length(celltypes), ncol = length(celltypes)))
rownames(df) <- celltypes
colnames(df) <- celltypes
df_cor <- df
df_p <- df

for (ct1 in celltypes) {
  col_id_tested <- c() # this is to help adjust p value (based on p values that are calculated)
  for (ct2 in celltypes) {
    if(ct2 == ct1){
      next
    }
    col_id_tested <- c(col_id_tested, ct2)
    
    qtls1 <- qtls_ctList[[ct1]]
    qtls2 <- qtls_ctList[[ct2]]

    qtls.test <- union(qtls1, qtls2) %>% unique()
    
    effects1 <- df_effect[qtls.test, ct1]
    effects2 <- df_effect[qtls.test, ct2]
    
    xx <- cor.test(effects1, effects2, method = 'spearman')
    
    df_cor[ct1, ct2] <- xx$estimate %>% as.numeric()
    df_p[ct1, ct2] <- xx$p.value %>% as.numeric()
  }
  p_values_original <- df_p[ct1, col_id_tested]
  p_values_adjusted <- p.adjust(p_values_original, method = 'BH')
  df_p[ct1, col_id_tested] <- p_values_adjusted
}

# ##
# saveRDS(df_cor, paste0(out.dir, 'df_cor.rds'))
# saveRDS(df_p, paste0(out.dir, 'df_p.rds'))

df_cor <- df_cor %>% na.replace(0)
df_p <- df_p %>% na.replace(0)

# change non-significant cor values to 0
df_cor_plot <- df_cor
adjp_threshold <- 0.05
for (ct1 in celltypes){
  for (ct2 in celltypes){
    if(df_p[ct1, ct2] > adjp_threshold){
      df_cor_plot[ct1, ct2] <- 0
    }
  }
}

##################  ################ ##################  ################
# plot
color = c('#fff7ec', '#fc8d59', '#d7301f')
color = c('white', '#fff7ec','#fee8c8','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#b30000','#7f0000')
breaks <- seq(0.2, 1, length.out = length(color))
col_fun <- circlize::colorRamp2(breaks, color)

library(circlize)
celltype_order <- names(gca_colors_atac_union_anno1)[names(gca_colors_atac_union_anno1) %in% celltypes]
celltype_order <- c("TI_Stem","AC_Stem","Early_Enterocyte","Enterocyte","Early_Colonocyte","Colonocyte","TI_Goblet","AC_Goblet","BEST4",
                    "CD4T","CD8T","gdT","MemoryB","Plasma","NaiveB","Macrophage","Neutrophil","Fibroblast","Endothelium")

corr_mtx <- df_cor_plot[celltype_order, celltype_order]

column_ha <- HeatmapAnnotation(foo = anno_text(celltype_order, location = 1, rot = 45, 
                                               gp = gpar(fontsize = 10)))


pdf(paste0(out.dir, 'raw_effect_corr_FromEitherCT_QTL.pdf'), height = 7.5, width = 10)
ht <- Heatmap(corr_mtx, 
              name = paste0('spearman cor'),
              cell_fun = function(j, i, x, y, w, h, fill) {
                if(i < j ) {
                  grid.rect(x, y, w, h, gp = gpar(fill = 'white', col = 'white', border = T))
                }
              },
              bottom_annotation = column_ha,
              col = col_fun,
              rect_gp = gpar(col = "black", lwd = 0.5),
              cluster_columns = F, cluster_rows = F,
              show_row_dend = F, show_column_dend = F, 
              show_row_names = T, show_column_names = F,
              row_names_side = "left",
              border = F,
              use_raster = F,
              show_heatmap_legend = T)
draw(ht, padding = unit(c(.1, .1, .1, .1), "npc"))
dev.off()



















