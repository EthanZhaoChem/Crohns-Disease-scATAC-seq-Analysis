# 12 quantitative colors
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
library(plyr)
library(scales)

solarExtra_continuous <- readRDS( '~/yuzhao1/scripts/colors/solarExtra.rds')

manual_colors <- c('#fccde5','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928','#8dd3c7','#80b1d3', '#ccebc5', '#d9d9d9', '#a6cee3', '#ffed6f', '#35978f', '#bf812d', '#c51b7d', '#1b7837', '#f4a582', '#d1e5f0', '#053061', '#878787')

manual_colors_rc2_location <- c('AC'='#80b1d3', 'TI'='#8dd3c7', 'POU'='#fdb462', 'PP'='#b3de69')
manual_colors_rc2_anno2_grouped <- setNames(
  c("#72b9bf", "#72b9bf","#72b9bf", "#c9d764", "#c9d764"), 
  c("EC TI", "EC PP", "EC POU Group1", "EC POU Group2", "EC AC"))

manual_colors_composition1 <- manual_colors

colors51 <- pal_igv("default")(51)
colors51 <- sample(colors51, 51)

################ Gradient colors #####################
# 1. deep blue to deep red, middle is white
manual_colors_gradient1 <-c(
  '#08306b',
  '#08519c',
  '#2171b5',
  '#4292c6',
  '#6baed6',
  '#9ecae1',
  '#c6dbef',
  '#deebf7',
  '#f7fbff',
  '#ffeda0',
  '#fed976',
  '#feb24c',
  '#fd8d3c',
  '#fc4e2a',
  '#e31a1c',
  '#bd0026',
  '#800026'
)

# show_col(manual_colors_gradient1)

# 2. from warm yellow to warm red
manual_colors_gradient2 <-c('#ffffcc',
                            '#ffeda0',
                            '#fc4e2a',
                            '#e31a1c',
                            '#bd0026',
                            '#800026')

# show_col(manual_colors_gradient2)

# 3. from white to deep cold red
manual_colors_gradient3 <- c(
  '#FFFFFF',
  '#fff7ec',
  '#fee8c8',
  '#fdd49e',
  '#fdbb84',
  '#fc8d59',
  '#ef6548',
  '#d7301f',
  '#b30000',
  '#7f0000'
)
# show_col(manual_colors_gradient3)

# 3. from grey to deep cold red
fun_color_range <- colorRampPalette(c('#f0f0f0', '#a50f15'))
manual_colors_gradient4 <- fun_color_range(10)    



###################### Seurat related ##############################
# 1.1 seurat umap
plot_seurat_UMAP_custom <- function(seurat.obj, label, reduction_name){
  p <- DimPlot(seurat.obj, group.by = label, reduction = reduction_name, label = T, shuffle=T, repel = T)+
    theme(
      panel.background = element_rect(fill = 'white', colour = 'black'),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 15, face = 'bold'),
      legend.key = element_blank(),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 6),
      plot.margin = margin(2, 2, 2, 2, "cm"),
      plot.title = element_blank()
    ) +
    scale_color_manual(values = c(manual_colors, colors51)) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    labs(title = "UMAP", x = "UMAP Dimension_1", y = "UMAP Dimension_2")
  
  return(p)
}

# 1.2 seurat feature plot
plot_seurat_feature_custom <- function(seurat.obj, feature_to_plot, reduction_name,
                                       combine = T) {
  FeaturePlot(seurat.obj,
              order = F,
              reduction = reduction_name,
              features = feature_to_plot,
              pt.size = 0.1,
              max.cutoff = "q95",
              raster = F)+
    theme_pubr() +
    theme(
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 15, color = "black"),
      legend.text = element_text(size = 10, color = "black"),
      legend.title = element_text(size = 15, face = "bold"),
      legend.position = 'bottom',
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    )+    
    labs(title = feature_to_plot, x = "UMAP Dimension_1", y = "UMAP Dimension_2")+
    scale_color_continuous_sequential(palette = "Reds 2") +
    guides(color = guide_colorbar(
      barwidth = unit(6, "cm"),
      barheight = unit(0.4, "cm"),
    ))
}

# 1.3 seurat QC plot
plot_QC_seurat <- function (seurat.obj) {
  seurat.obj$xx <- 1
  VlnPlot(
    seurat.obj,
    group.by = 'xx',
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    pt.size = 0,
    ncol = 3
  ) + NoLegend()
}

# 1.4 seurat calculate gene percentage in all or clusters
PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}

# 1.5 seurat dot plot
 
plot_seurat_dot_custom <- function(seurat.obj, group, features, levels.custom = 'default'){
  K <- DotPlot(seurat.obj, group.by = group, features = features)
  x.color <- K$data$avg.exp
  x.color.scaled <- K$data$avg.exp.scaled
  x.size <- K$data$pct.exp
  x.gene <- K$data$features.plot
  x.cluster <- K$data$id
  x.size.max <- max(x.size)
  
  if(levels.custom != 'default'){
    x.cluster <- factor(x.cluster, levels = levels.custom)
  }
  
  x.df <- data.frame(x.color2 = x.color,
                     x.color.scaled2 = x.color.scaled,
                     x.size2 = x.size/x.size.max,
                     x.gene2 <- x.gene,
                     x.cluster2 = x.cluster)
  
  ggplot(x.df, aes(x=x.cluster2, y=x.gene2)) +
    geom_point(aes(fill = x.color.scaled2, size = x.size2), shape=21, color = 'black')+
    scale_fill_gradientn(colours = manual_colors_gradient1)+
    theme(
      panel.background = element_rect(fill = 'white', colour = 'black'),
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 10,  colour = 'black'),
      axis.text.x = element_text(size = 10, angle = 45,  hjust=1,  colour = 'black'),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.key = element_blank(),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 6),
      plot.margin = margin(1, 1, 1, 1, "cm"),
      plot.title = element_blank()
    ) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    labs(title = "", x = "Clusters", y = "Genes")
  
  
}


############################### dataframe related ########################
# 2.1 df for umap with cluster_name, embedding1, embedding2

# if plot feature: cap and threshold the values before input by quantile: e.g. 0.01, 0.99
plot_df_umap_custom <- function(df, plot_level  = 'default',
                                feature_title = 'default',
                                title = 'UMAP', 
                                custom_colors = 'default', show.label = 'number', 
                                plot_feature = F,
                                feature_name = 'default',
                                axistitle.textsize = 15,
                                legend.textsize = 6,
                                cap1 = 0.01,
                                cap2 = 0.99,
                                cap_customized = 'percentile'){
  # use feature_to_plot to specifity the feature column in df
  

  
  if(plot_level  == 'default'){
    plot_level <- names(sort(table(df$cluster_name), decreasing = T))
  }

  if (custom_colors == 'default'){
    custom_colors <- c(manual_colors, colors51)
  } 
  else{
    custom_colors <- mapvalues(plot_level, from = names(custom_colors), to = custom_colors) %>% unlist()
  }
  
  if(plot_feature == T){
    if(cap_customized == 'percentile'){
      df$feature_to_plot[df$feature_to_plot < quantile(df$feature_to_plot, cap1)] <- quantile(df$feature_to_plot, cap1)
      df$feature_to_plot[df$feature_to_plot > quantile(df$feature_to_plot, cap2)] <- quantile(df$feature_to_plot, cap2)
    }
    else{
      df$feature_to_plot[df$feature_to_plot < cap_customized[[1]]] <- cap_customized[[1]]
      df$feature_to_plot[df$feature_to_plot > cap_customized[[2]]] <- cap_customized[[2]]
    }
    umap <- ggplot(df, aes(embedding1, embedding2)) +
      geom_point(aes(colour = feature_to_plot), size = 0.03) +
      theme(
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = axistitle.textsize, face = 'bold'),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = legend.textsize),
        plot.margin = margin(2, 2, 2, 2, "cm"),
        plot.title = element_text(size = axistitle.textsize, face = 'bold')
      ) +
      scale_color_gradientn(colours = c(manual_colors_gradient4))+
      labs(title = feature_name, x = paste0(title, " Dimension_1"), y = paste0(title, " Dimension_2"))
    return(umap)
  }

  # shuffle plot sequence
  df <- df[sample(nrow(df)),]
  class_avg <- df %>%
    group_by(cluster_name) %>%
    summarise_at(vars(embedding1, embedding2), list(median))
  
  # from big to small clusters, rank starting from 1 
  class_avg$ranking <- mapvalues(class_avg$cluster_name, plot_level, 1:length(plot_level))
  class_avg$legend <- paste0(class_avg$ranking, '-', class_avg$cluster_name)
  df$ranking <- mapvalues(df$cluster_name, from = class_avg$cluster_name, to = class_avg$ranking)
  df$legend <- mapvalues(df$cluster_name, from = class_avg$cluster_name, to = class_avg$legend)
  df$legend <- factor(df$legend, levels = paste0(1:length(plot_level), '-', plot_level))
  df$cluster_name <- factor(df$cluster_name, levels = plot_level)
  
  # default, show number on plot
  if(show.label == 'number'){
    umap <- ggplot(df, aes(embedding1, embedding2)) +
      rasterize(geom_point(aes(colour = legend), size = 0.001, alpha = 0.8), dpi = 300, scale = 0.6) +
      geom_label_repel(
        aes(label = ranking),
        data = class_avg,
        label.size = NA,
        # fontface = "bold",
        max.overlaps = Inf,
        fill = NA
      ) +
      theme(
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = axistitle.textsize, face = 'bold'),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = legend.textsize),
        plot.margin = margin(2, 2, 2, 2, "cm"),
        plot.title = element_blank()
      ) +
      scale_color_manual(values = custom_colors) +
      guides(colour = guide_legend(override.aes = list(size = legend.textsize, shape = 15))) +
      labs(title = title, x = paste0(title, " Dimension_1"), y = paste0(title, " Dimension_2"))
  }
  
  # show name on plot
  else if(show.label == 'name') 
  {
    umap <- ggplot(df, aes(embedding1, embedding2)) +
      rasterize(geom_point(aes(colour = cluster_name), size = 0.001, alpha = 0.8), dpi = 300, scale = 0.6) +
      geom_label_repel(
        aes(label = cluster_name),
        data = class_avg,
        label.size = NA,
        # fontface = "bold",
        max.overlaps = Inf,
        fill = NA
      ) +
      theme(
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = axistitle.textsize, face = 'bold'),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = legend.textsize),
        plot.margin = margin(2, 2, 2, 2, "cm"),
        plot.title = element_blank()
      ) +
      scale_color_manual(values = custom_colors) +
      guides(colour = guide_legend(override.aes = list(size = legend.textsize, shape = 15))) +
      labs(title = title, x = paste0(title, " Dimension_1"), y = paste0(title, " Dimension_2"))
  }
  
  # if user specify random parameter, nothing will show up
  else{ # if user specify random parameter, nothing will show up
    umap <- ggplot(df, aes(embedding1, embedding2)) +
      rasterize(geom_point(aes(colour = cluster_name), size = 0.001, alpha = 0.8), dpi = 300, scale = 0.6) +
      theme(
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = axistitle.textsize, face = 'bold'),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = legend.textsize),
        plot.margin = margin(2, 2, 2, 2, "cm"),
        plot.title = element_blank()
      ) +
      scale_color_manual(values = custom_colors) +
      guides(colour = guide_legend(override.aes = list(size = legend.textsize, shape = 15))) +
      labs(title = title, x = paste0(title, " Dimension_1"), y = paste0(title, " Dimension_2"))
  }

  return(umap)
}

# 2.2 df for stacked col plot
plot_stacked_cols <- function(df, x, y, fill){
  p <- ggplot(df, aes(x = x, y = y, fill = fill)) +
    geom_col(width = 0.6) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 10, face = 'bold', colour = 'black', angle = 90, vjust = 0.5, hjust=1),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.key = element_blank(),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      plot.margin = margin(2, 2, 2, 2, "cm"),
      plot.title = element_blank()
    )
  return(p)
}

























