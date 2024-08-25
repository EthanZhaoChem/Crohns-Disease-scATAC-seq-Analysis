# make sure this script is shared across three folders: atac_gca2024, final_GCArna, final_GCAatac
# original file is in final_GCArna, other folders have soft links
library(ggsci)
library(ggpubr)
library(colorspace)
library(paletteer) 
library(wesanderson)
library(ggrepel)
library(ggrastr)
library(ggforce)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(scales)

############### archived colors, the new ones are in plot.R ##############
gca_colors_peakType <- c('Promoter' = '#fc8d62', 'Intronic' = '#66c2a5', 'Exonic' = '#8da0cb', 'Distal' = '#e78ac3')
gca_colors_snpType <- c('Promoter' = '#fc8d62', 
                        'Exonic' = '#8da0cb',
                        'Intronic' = '#66c2a5', 
                        'UTR' = '#a6d854',
                        'Intergenic Proximal' = '#ffd92f',
                        'Intergenic Distal' = '#e78ac3')

gca_colors_location <- c(
  'TI'='#8dd3c7', 'AC'='#80b1d3'
)
gca_colors_disease <- c('Control' = '#33a02c', 'CD' = '#e7298a')

gca_colors_inflammation <- c('Control' = '#b3de69',
                             'nonInf' = '#bebada',
                             'adjInf' = '#fdb462',
                             'inf' = '#fb8072')


############ 1, atac umap ##############
gca_colors_atac_lineage <- c('epithelial' = '#72b9bf',
                             'immune' = '#82c07f',
                             'stromal' = '#efb788')


gca_colors_atac_epithelial_anno1 <- c(
  "TI_Stem" = '#fccde5',
  "AC_Stem" = '#ff7f00',
  "Early_Enterocyte" = '#9ecae1',
  "Enterocyte" = '#4292c6',
  "Early_Colonocyte" = '#80cdc1',
  "Colonocyte" = '#35978f',
  "TI_Goblet" = '#fbb4ae',
  "AC_Goblet" = '#fed9a6',
  "BEST4" = '#b3de69',
  "Tuft" = '#ffed6f',
  "EEC" = '#fb8072',
  "Paneth" = '#cd62af')


gca_colors_atac_immune_anno1 <- c("CD4T" = '#b2df8a', # green
                            "CD8T" = '#fc4e2a', # orange
                            "gdT" = '#bc80bd', # purple 
                            "NK" = '#de77ae', # red purple
                            "ILCs" = '#8c510a', # brown
                            "GC_B" = '#ccebc5', # blue green
                            "NaiveB" = '#b3de69', 
                            "MemoryB" = '#80b1d3', 
                            "Plasma" = '#8dd3c7', 
                            "Mast" = '#cfa061', 
                            "Macrophage" = '#4393c3', 
                            "DC" = '#542788', # purple
                            "Neutrophil" = '#db433a') # cute red


gca_colors_atac_stromal_anno1 <- c("Fibroblast" = '#bf812d',
                             "Pericyte" = '#fccde5',
                             "Endothelium" = '#ffed6f',
                             "Glial" = '#ec5f60')


gca_colors_atac_union_anno1 <- c(gca_colors_atac_epithelial_anno1,
                                 gca_colors_atac_immune_anno1,
                                 gca_colors_atac_stromal_anno1)

############ heatmap ##############
gca_heatmap_colors_gradient1 <-c(
  '#08306b',
  # '#08519c',
  '#2171b5',
  # '#4292c6',
  # '#6baed6',
  '#9ecae1',
  '#c6dbef',
  # '#deebf7',
  '#f7fbff',
  '#ffeda0',
  # '#fed976',
  '#feb24c',
  '#fd8d3c',
  '#fc4e2a',
  '#e31a1c',
  '#bd0026',
  '#800026'
)

orders_anno1 <- c(
  names(gca_colors_atac_epithelial_anno1),
  names(gca_colors_atac_immune_anno1),
  names(gca_colors_atac_stromal_anno1)
)
