# 1.
# convert markerTest of peaks to df_stat
archr_helper_markerPeaks_converter <- function(markertest){
  stat.colnames <- c("Log2FC", "Mean", "FDR", "Pval", "MeanDiff", "AUC", "MeanBGD")
  df_stat <- as.data.frame(matrix(nrow = nrow(markertest), ncol = length(stat.colnames)))
  colnames(df_stat) <- stat.colnames
  for (stat.name in stat.colnames) {
    df_stat[[stat.name]] <- markertest@assays@data[[stat.name]]$x
  }
  peakLoc <- paste0(markertest@elementMetadata$seqnames,  ':',
                    markertest@elementMetadata$start,'-',markertest@elementMetadata$end)
  rownames(df_stat) <- peakLoc
  df_stat <- df_stat[complete.cases(df_stat[,c("Log2FC","FDR")]),]
  return(df_stat)
}
# # example usage
# df_stat <- archr_helper_markerPeaks_converter(markersPeaks)
# archr_peakset_Loc <- paste0(paste0(as.character(seqnames(proj@peakSet)), ':',
#                                    start(proj@peakSet),'-',end(proj@peakSet)))
# idx_matched <- match(rownames(df_stat), archr_peakset_Loc)
# peaks_info <- as.data.frame(proj@peakSet[idx_matched], row.names = NULL)
# df_stat <- bind_cols(df_stat, peaks_info)
# df_stat_significant <- df_stat[which(df_stat$FDR<0.001 & df_stat$Log2FC > 1), ]

# 2. convert marker peaks to dfs
archr_helper_markerPeaks_converter_multiple <- function(markertest){
  stats.all <- list()
  stat.colnames <- c("Log2FC", "Mean", "FDR", "Pval", "MeanDiff", "AUC", "MeanBGD")
  
  for(cluster.ID in 1: length(markertest@colData@rownames)){
    cluster.name <- markertest@colData@rownames[[cluster.ID]]
    df_stat <- as.data.frame(matrix(nrow = nrow(markertest), ncol = length(stat.colnames)))
    colnames(df_stat) <- stat.colnames
    for (stat.name in stat.colnames) {
      df_stat[[stat.name]] <- markertest@assays@data[[stat.name]][, cluster.ID]
    }
    peakLoc <- paste0(markertest@elementMetadata$seqnames,  ':',
                      markertest@elementMetadata$start,'-',markertest@elementMetadata$end)
    rownames(df_stat) <- peakLoc
    df_stat <- df_stat[complete.cases(df_stat[,c("Log2FC","FDR")]),]
    
    stats.all[[cluster.name]] <- df_stat
  }
  return(stats.all)
}

# # example usage
# markersPeaks <- readRDS(paste0(proj@projectMetadata$outputDirectory, '/', 'MarkersPeaks_', temp_lineage, '.rds'))
# df_stat_ByCluster <- archr_helper_markerPeaks_converter_multiple(markersPeaks)
# 
# for (cluster.name in names(df_stat_ByCluster)){
#   df_stat <- df_stat_ByCluster[[cluster.name]]
#   df_stat_significant <- df_stat[which(df_stat$FDR<cutoff_FDR & df_stat$Log2FC > cutoff_Log2FC), ]
#   DARs_plot[[cluster.name]] <- rownames(df_stat_significant)
# }


# 3.
# convert markerTest of TF motifs to df_stat
archr_helper_markerTFs_converter <- function(markertest){
  stat.colnames <- c("mlog10Padj", "mlog10p", "Enrichment", "BackgroundProporition", "nBackground", "BackgroundFrequency",
                     "CompareProportion", "nCompare", "CompareFrequency", "feature" )
  df_stat <- as.data.frame(matrix(nrow = nrow(markertest), ncol = length(stat.colnames)))
  colnames(df_stat) <- stat.colnames
  for (stat.name in stat.colnames) {
    df_stat[[stat.name]] <- markertest@assays@data[[stat.name]][[1]]
  }
  return(df_stat)
}

# 4.
# convert markerTest of TF chromVAR to df_stat
archr_helper_markerTFs_chromvar_converter <- function(markertest){
  stat.colnames <- c('Mean', 'FDR', 'Pval', 'MeanDiff', 'AUC', 'MeanBGD')
  df_stat <- as.data.frame(matrix(nrow = nrow(markertest), ncol = length(stat.colnames)))
  rownames(df_stat) <- markertest@elementMetadata$name
  colnames(df_stat) <- stat.colnames
  for (stat.name in stat.colnames) {
    df_stat[[stat.name]] <- markertest@assays@data[[stat.name]][[1]]
  }
  df_stat$TF <- rownames(df_stat)
  return(df_stat)
}

# 5
# plot peak summary
plotPeakCallSummary <- function(
  ArchRProj = NULL, 
  pal = NULL
){
  
  peakDF <- metadata(ArchRProj@peakSet)$PeakCallSummary
  
  if(is.null(peakDF)){
    stop("Error no Peak Call Summary available are you sure these peaks were called with CreateReproduciblePeakSet?")
  }
  
  if(is.null(pal)){
    pal <- paletteDiscrete(values=peakDF$Var1)
  }
  pal <- c("Distal" = "#b3de69", "Exonic" = "#73C6FF", "Intronic" = "#abd9e9", "Promoter" = "#FFC554", pal) 
  pal <- pal[!duplicated(names(pal))]
  
  lengthMax <- split(peakDF$Freq, peakDF$Group) %>% lapply(sum) %>% unlist %>% max
  
  colnames(peakDF)[colnames(peakDF)=="Var1"] <- "PeakAnno"
  
  # sort based on #peak
  xlabel_levels <- split(peakDF$Freq, peakDF$Group) %>% lapply(sum) %>% unlist %>% sort(decreasing = T) %>% names()
  peakDF$Group <- factor(peakDF$Group, levels = xlabel_levels)
  
  p <- ggplot(peakDF, aes(x=Group, y=Freq, fill=PeakAnno)) + 
    geom_bar(stat = "identity") + 
    theme_ArchR(xText90 = TRUE) +
    ylab("Number of Peaks (x10^3)") +
    xlab("") + 
    theme(
      axis.text.x = element_text(angle=45, hjust = 1, vjust = 1),
      legend.position = "bottom", 
      legend.key = element_rect(size = 2), 
      legend.box.background = element_rect(color = NA)
    ) +
    scale_fill_manual(values=pal) +
    scale_y_continuous(
      breaks = seq(0, lengthMax * 2,50), 
      limits = c(0, lengthMax * 1.1), 
      expand = c(0,0)
    )
  attr(p, "ratioYX") <- 0.5
  return(p)
}

# # example usage
# for (temp_lineage in c('epithelial', 'tcell', 'bcell', 'myeloid', 'others')) {
#   proj <- paste0('proj_', temp_lineage) %>% as.name(.) %>% eval(.)
#   ideal_width <- metadata(proj@peakSet)$PeakCallSummary[['Group']] %>% unique() %>% length() * 0.3+2
#   
#   p <- plotPeakCallSummary(proj)
#   pdf(paste0(out.dir, temp_lineage, '_peaks', '.pdf'), width = ideal_width, height = 6, pointsize = 1)
#   print(p)
#   dev.off()
# }




# 5
# hyper-geometric enrichment of candidate peaks
# caution: the test to obtain candidate peaks should be done in union peakset (not subsetted)

archr_customized_motif_enrichment <- function(ArchRProj = NULL,
                                              peakAnnotation = "Motif",
                                              candidate_peaks_vector = NULL){
  
  candidate_peaks <- candidate_peaks_vector %>% gsub(':', '_', .) %>% gsub('-', '_', .) 
  matches <- getMatches(ArchRProj, peakAnnotation)
  bgdPeaks <- SummarizedExperiment::assay(getBgdPeaks(ArchRProj))
  
  
  r1 <- SummarizedExperiment::rowRanges(matches)
  pr1 <- paste(seqnames(r1),start(r1),end(r1),sep="_")
  mcols(r1) <- NULL
  
  # candidate peak index
  idx <- which(pr1 %in% candidate_peaks)
  compare <-idx
  
  # bgd peaks: candidate peaks and similar peaks (indexed by themselves)
  background = c(idx, as.vector(bgdPeaks[idx,]))
  
  # the motif annotation matrix
  matches <- matches@assays@data$matches
  
  #Compute Totals
  matchCompare <- matches[compare, ,drop=FALSE]
  matchBackground <- matches[background, ,drop=FALSE]
  matchCompareTotal <- Matrix::colSums(matchCompare)
  matchBackgroundTotal <- Matrix::colSums(matchBackground)
  
  #Create Summary DF
  pOut <- data.frame(
    feature = colnames(matches),
    CompareFrequency = matchCompareTotal,
    nCompare = nrow(matchCompare),
    CompareProportion = matchCompareTotal/nrow(matchCompare),
    BackgroundFrequency = matchBackgroundTotal,
    nBackground = nrow(matchBackground),
    BackgroundProporition = matchBackgroundTotal/nrow(matchBackground)
  )
  
  pOut$Enrichment <- pOut$CompareProportion / pOut$BackgroundProporition
  
  #Get P-Values with Hyper Geometric Test
  pOut$mlog10p <- lapply(seq_len(nrow(pOut)), function(x){
    p <- -phyper(pOut$CompareFrequency[x] - 1, # Number of Successes the -1 is due to cdf integration
                 pOut$BackgroundFrequency[x], # Number of all successes in background
                 pOut$nBackground[x] - pOut$BackgroundFrequency[x], # Number of non successes in background
                 pOut$nCompare[x], # Number that were drawn
                 lower.tail = FALSE, log.p = TRUE)# P[X > x] Returns LN must convert to log10
    return(p/log(10))
  }) %>% unlist %>% round(4)
  
  #Minus Log10 Padj
  pOut$mlog10Padj <- pmax(pOut$mlog10p - log10(ncol(pOut)), 0)
  pOut <- pOut[order(pOut$mlog10p, decreasing = TRUE), , drop = FALSE]
  
  # return results
  return(pOut)
}

# # example usage
# xx <- archr_customized_motif_enrichment(proj, peakAnnotation = "Motif", candidate_peaks_vector = DARs_plot$EC_ACvsTI)





















































