library(fastTopics)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(plyr)
library(dplyr)
library(stringr)
library(ArchR)
library(fitdistrplus)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/work/final_GCArna/scripts/gca_colors.R')
source('~/yuzhao1/scripts/helper_archr.R')

# fixed
lineage <- 'immune'
lineage.anno1 <- gca_colors_atac_immune_anno1
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1_subSampling_max100_topic_immune/")
data.dir <- '~/yuzhao1/work/atac_gca2024/13fasttopic/rds/'
load(paste0(data.dir, "immune_raw_sub100.RData"))
thrP <- 0.95 # 	Probability threshold to use as cutoff on the probability distribution when using GammaFit as method.
subID <- '100'

# change for each k
nTopics <- 20

# fixed
file_fasttopic <- paste0('~/yuzhao1/work/atac_gca2024/13fasttopic/rds/fit_immune_sub100_k',  nTopics, '_converged.rds')
file_fasttopic_de <- paste0('~/yuzhao1/work/atac_gca2024/13fasttopic/rds/fit_immune_sub100_k',  nTopics, '_converged_de_default.rds')
plot.dir <- paste0('~/yuzhao1/work/atac_gca2024/13fasttopic/plots/immune_sub100_k', nTopics, '/')
dir.create(plot.dir, showWarnings = T)


################################## read topic-peak matrix ##############################################
# fit$F and fit$L are the parameters of the multinomial topic model; specifically, fit$L[i,] gives the topic probabilities for sample or document i, and fit$F[,k] gives the term probabilities for topic k.
fit <- readRDS(file_fasttopic)
fit <- poisson2multinom(fit)

# make sure the cell names are in the correct order
metadata1 <- as.data.frame(proj@cellColData)
metadata1 <- metadata1[rownames(fit$L),]

# matrix export
peak_topic_mtx <- fit$F # all peaks in one topic sum to 1
cell_topic_mtx <- fit$L # all topics in one cell sum to 1

################################## topic-topic corr mtx ##############################################
corr_mtx <- matrix(0, nrow = nTopics, ncol = nTopics)
rownames(corr_mtx) <- paste0('k', 1:nTopics)
colnames(corr_mtx) <- paste0('k', 1:nTopics)
for (i in 1:nTopics) {
  for (j in 1:nTopics) {
    corr_mtx[i,j] <- cor(cell_topic_mtx[, i], cell_topic_mtx[, j], method = 'spearman')
  }
}

# heatmap plot
library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("#4393c3", "white", "#d6604d"))
col_fun(seq(-1, 1))

column_ha <- HeatmapAnnotation(foo = anno_text(colnames(corr_mtx), location = 1, rot = 30,
                                               gp = gpar(fontsize = 11, fontface='bold')))

pdf(paste0(plot.dir, lineage, '_sub', subID, '_k_',  nTopics, '_topic_corr.pdf'), height = 7, width = 11)
ht <- Heatmap(corr_mtx,
              bottom_annotation = column_ha,
              col = col_fun,
              rect_gp = gpar(col = "black", lwd = 2),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", corr_mtx[i, j]), x, y, gp = gpar(fontsize = 10))
              },
              cluster_columns = F, cluster_rows = F,
              show_row_dend = F, show_column_dend = F,
              show_row_names = T, show_column_names = F,
              row_names_side = "left",
              border = T,
              use_raster = F,
              show_heatmap_legend = FALSE)
lgd <- Legend(col_fun = col_fun, title = "Spearman correlation score")
draw(ht, heatmap_legend_list = lgd)
dev.off()

#################### structure plot: anno1 ##################
# pre sets
labels <- factor(metadata1$anno1, levels = names(lineage.anno1))
topic_colors <- colors51

# fitting progress: check convergence
pdf(paste0(plot.dir, 'progress_plot_', lineage, '_sub', subID, '_k_', nTopics, '_anno1.pdf'), height = 7, width = 7)
plot_progress(fit,x = "iter",add.point.every = 5,colors = "black") +
  theme_cowplot(font_size = 10)
dev.off()

# structure plot
pdf(paste0(plot.dir, 'structure_plot_', lineage, '_sub', subID, '_k_', nTopics, '_anno1.pdf'), height = 7, width = 26)
structure_plot(fit,
               grouping = labels,
               colors = topic_colors,
               # topics = 1:10,
               gap = 20, 
               perplexity = 50,
               verbose = FALSE)
dev.off()


#################### structure plot: partial topics ##################
# only need to customize the topics you want to show
partial_topics <- c(9, 10, 13)

# the last merged topic is white now
labels <- factor(metadata1$anno1, levels = names(lineage.anno1))
colors_pool <- c("gold","darkorange","dodgerblue", colors51)
topic_colors <- rep('white', 1+length(partial_topics))
topic_colors[1:length(topic_colors)-1] <- colors_pool[1:length(topic_colors)-1]


fit2 <- readRDS(file_fasttopic)
fit2 <- poisson2multinom(fit2)
fit2 <- merge_topics(fit2, paste0("k",setdiff(1:nTopics, partial_topics))) # combined topics become the last one
colnames(fit2$L) <- c(paste0('k', partial_topics), "other")


set.seed(2)
pdf(paste0(plot.dir, 'structure_plot_', lineage, '_sub', subID, '_k_', nTopics, '_anno1_partial.pdf'), height = 2, width = 8)
p <- structure_plot(fit2,
                    grouping = labels,
                    colors = topic_colors, # this order is the order of colnames(fit2$L), unrelated to the topic order below
                    topics = c(length(partial_topics)+1, 1:length(partial_topics)), # the first one is shown on top
                    gap = 50,
                    verbose = FALSE)
print(p)
dev.off()



################################## DA peaks analysis ##########################
## get topic specific peaks by DA analysis
DA_res <- readRDS(file_fasttopic_de)
lpval_threshold <- -log10(1)
postmean_threshold <- 0
nPeaks_selected <- 30000
summary(DA_res)
dim(DA_res$z)

# volcano:  lpval vs lfc
plots <- vector("list", nTopics)
names(plots) <- 1:nTopics
for (k in 1:nTopics) {
  data <- data.frame(lpval = DA_res$lpval[,k])
  x_vertical <- sort(data$lpval, decreasing = T)[[nPeaks_selected]]
  plots[[k]] <- ggplot(data, aes(x = lpval)) +
    geom_histogram(bins = 1000, fill = "blue", color = "black") +
    geom_vline(xintercept = x_vertical, color = "red", linetype = "dashed", size = 1) +
    xlim(0, 3) +
    ggtitle(paste("Topic", k))
}
png(paste0(plot.dir, 'da_plot_sub', subID, '_k_', nTopics, '_lpval.png'), height = 4000, width = 7200, res = 300)
do.call(plot_grid, c(plots, nrow = 5, ncol = 5))
dev.off()


# Number of regions selected at different lfsr cutoffs:
sig_regions <- matrix(NA, nrow = nTopics, ncol = 10)
colnames(sig_regions) <- c("lpval > -log10(1)",
                           "lpval > -log10(0.5)",
                           "lpval > -log10(0.4)",
                           "lpval > -log10(0.3)",
                           "lpval > -log10(0.2)", 
                           "lpval > -log10(0.1)", 
                           "lpval > -log10(0.05)", 
                           "lpval > -log10(0.01)", 
                           "lpval > -log10(0.005)", 
                           "lpval > -log10(0.001)")
rownames(sig_regions) <- paste("topic", 1:nrow(sig_regions))

for(k in 1:nTopics){
  lpval <- DA_res$lpval[,k]
  postmean <- DA_res$postmean[,k]
  sig_regions[k, ] <- c(length(which((lpval > -log10(1)) & (postmean > 0))),
                        length(which((lpval > -log10(0.5)) & (postmean > 0))),
                        length(which((lpval > -log10(0.4)) & (postmean > 0))),
                        length(which((lpval > -log10(0.3)) & (postmean > 0))),
                        length(which((lpval > -log10(0.2)) & (postmean > 0))),
                        length(which((lpval > -log10(0.1)) & (postmean > 0))),
                        length(which((lpval > -log10(0.05)) & (postmean > 0))),
                        length(which((lpval > -log10(0.01)) & (postmean > 0))),
                        length(which((lpval > -log10(0.005)) & (postmean > 0))),
                        length(which((lpval > -log10(0.001)) & (postmean > 0)))
  )
}

sig_regions


# postmean > 0 as criteria to select regions
da_peaks <- list()
for(k in 1:nTopics){
  lpval <- DA_res$lpval[,k]
  # here the lpval_threshold is for top 10k peaks sorted by lpval
  lpval_threshold <- sort(lpval, decreasing = T)[[nPeaks_selected]]
  postmean <- DA_res$postmean[, k]
  da_peaks[[k]] <- rownames(DA_res$lpval)[which((lpval > lpval_threshold) & (postmean > postmean_threshold))]
}

daPeaks.dir <- paste0(plot.dir, 'daPeaks_positive/')
dir.create(daPeaks.dir)
# save da peak to bed files
for (i in 1:nTopics){
  peaks <- da_peaks[[i]]
  df <- data.frame(chr = strsplit(peaks, '_') %>% sapply(.,`[[`,1),
                   start = strsplit(peaks, '_') %>% sapply(.,`[[`,2),
                   end = strsplit(peaks, '_') %>% sapply(.,`[[`,3))
  write.table(df, paste0(daPeaks.dir, 'daPeaks_Topic', i, '.bed'),
              row.names = F,col.names = F, sep="\t", quote=FALSE)
}

motifs_perTopic_daPeaks <- list()
for (i in 1:nTopics) {
  motifs_perTopic_daPeaks[[i]] <- archr_customized_motif_enrichment(proj, peakAnnotation = "Motif", da_peaks[[i]])
}

library(openxlsx)
names(motifs_perTopic_daPeaks) <- paste0('Topic', 1:nTopics, '_', lengths(da_peaks), 'peaks')
write.xlsx(motifs_perTopic_daPeaks, file = paste0(plot.dir, 'enriched_motif_from_daPeaks_positive_', 'sub', subID, '_k_', nTopics, '.xlsx'))






################################## DA peaks analysis: from parallel ##########################
lpval_threshold <- -log10(1)
postmean_threshold <- 0
nPeaks_selected <- 30000

parallel_results_dir <- '~/yuzhao1/work/atac_gca2024/13fasttopic/rds/fit_immune_sub100_k20_converged_de_vsnull/'
nThreads <- 20
parallel_results <- list()
pval_all <- matrix(0, nrow = 0, ncol = nTopics)
postmean_all <- matrix(0, nrow = 0, ncol = nTopics)
for (i in 1:nThreads) {
  parallel_results[[i]] <- readRDS(paste0(parallel_results_dir, i, '.rds'))
  pval_all <- rbind(pval_all, parallel_results[[i]]$lpval)
  postmean_all <- rbind(postmean_all, parallel_results[[i]]$postmean)
}

# volcano:  lpval vs lfc
plots <- vector("list", nTopics)
names(plots) <- 1:nTopics
for (k in 1:nTopics) {
  data <- data.frame(lpval = pval_all[,k])
  x_vertical <- sort(data$lpval, decreasing = T)[[nPeaks_selected]]
  plots[[k]] <- ggplot(data, aes(x = lpval)) +
    geom_histogram(bins = 1000, fill = "blue", color = "black") +
    geom_vline(xintercept = x_vertical, color = "red", linetype = "dashed", size = 1) +
    xlim(0, 50) +
    ggtitle(paste("Topic", k))
}
png(paste0(plot.dir, 'da_vsnull_plot_sub', subID, '_k_', nTopics, '_lpval.png'), height = 4000, width = 7200, res = 300)
do.call(plot_grid, c(plots, nrow = 5, ncol = 5))
dev.off()


# Number of regions selected at different lfsr cutoffs:
sig_regions <- matrix(NA, nrow = nTopics, ncol = 10)
colnames(sig_regions) <- c("lpval > -log10(1)",
                           "lpval > -log10(0.5)",
                           "lpval > -log10(0.4)",
                           "lpval > -log10(0.3)",
                           "lpval > -log10(0.2)", 
                           "lpval > -log10(0.1)", 
                           "lpval > -log10(0.05)", 
                           "lpval > -log10(0.01)", 
                           "lpval > -log10(0.005)", 
                           "lpval > -log10(0.001)")
rownames(sig_regions) <- paste("topic", 1:nrow(sig_regions))

for(k in 1:nTopics){
  lpval <- pval_all[,k]
  postmean <- postmean_all[,k]
  sig_regions[k, ] <- c(length(which((lpval > -log10(1)) & (postmean > 0))),
                        length(which((lpval > -log10(0.5)) & (postmean > 0))),
                        length(which((lpval > -log10(0.4)) & (postmean > 0))),
                        length(which((lpval > -log10(0.3)) & (postmean > 0))),
                        length(which((lpval > -log10(0.2)) & (postmean > 0))),
                        length(which((lpval > -log10(0.1)) & (postmean > 0))),
                        length(which((lpval > -log10(0.05)) & (postmean > 0))),
                        length(which((lpval > -log10(0.01)) & (postmean > 0))),
                        length(which((lpval > -log10(0.005)) & (postmean > 0))),
                        length(which((lpval > -log10(0.001)) & (postmean > 0)))
  )
}

sig_regions

# postmean > 0 as criteria to select regions
da_peaks <- list()
for(k in 1:nTopics){
  lpval <- pval_all[,k]
  # here the lpval_threshold is for top 10k peaks sorted by lpval
  lpval_threshold <- sort(lpval, decreasing = T)[[nPeaks_selected]]
  postmean <- postmean_all[, k]
  da_peaks[[k]] <- rownames(pval_all)[which((lpval > lpval_threshold) & (postmean > postmean_threshold))]
}

daPeaks.dir <- paste0(plot.dir, 'daPeaks_positive_flexibleLpval30k_vsnull/')
dir.create(daPeaks.dir)
# save da peak to bed files
for (i in 1:nTopics){
  peaks <- da_peaks[[i]]
  df <- data.frame(chr = strsplit(peaks, '_') %>% sapply(.,`[[`,1),
                   start = strsplit(peaks, '_') %>% sapply(.,`[[`,2),
                   end = strsplit(peaks, '_') %>% sapply(.,`[[`,3))
  write.table(df, paste0(daPeaks.dir, 'daPeaks_Topic', i, '.bed'),
              row.names = F,col.names = F, sep="\t", quote=FALSE)
}

motifs_perTopic_daPeaks <- list()
for (i in 1:nTopics) {
  motifs_perTopic_daPeaks[[i]] <- archr_customized_motif_enrichment(proj, peakAnnotation = "Motif", da_peaks[[i]])
}

library(openxlsx)
names(motifs_perTopic_daPeaks) <- paste0('Topic', 1:nTopics, '_', lengths(da_peaks), 'peaks')
write.xlsx(motifs_perTopic_daPeaks, file = paste0(plot.dir, 'enriched_motif_from_daPeaks_vsnull_positive_', 'sub', subID, '_k_', nTopics, '.xlsx'))





################### topic score vln distribution in different cell groups, condition, location  ##########################
df <- bind_cols(cell_topic_mtx, metadata1[rownames(cell_topic_mtx),])
df$anno1 <- factor(df$anno1, levels = names(gca_colors_atac_union_anno1))
df$inflammation_status<- factor(df$inflammation_status, levels = names(gca_colors_inflammation))

## disease
topic_plots <- list()
for (i in 1:nTopics) {
  pp <- ggviolin(df, x = 'anno1', y = paste0('k', i),   fill = 'disease_status', palette = gca_colors_disease,
                 scale = 'width', width=0.8, 
                 trim = T, position = position_dodge(0.7), draw_quantiles = 0.5)+
    labs( x = NULL) +
    theme_pubr()+
    theme(axis.text.y = element_text(size=15),
          axis.text.x = element_text(angle=45, hjust = 1, vjust = 1),
          axis.title = element_text(size=20),
          legend.title = element_blank(),
          legend.text = element_text(size=10),
          legend.position = "bottom",
          plot.title = element_text(size=20, hjust=0.5, face = 'bold'))+
    labs(title = '')+
    facet_wrap(~ biopsy_location, nrow = 2, strip.position = 'left') +
    theme(
      strip.background = element_rect(fill = "white", colour = "white"),
      strip.text = element_text(size = 15, face = 'bold'),
      strip.placement = "outside"
    )
  
  topic_plots[[i]] <- pp
}

pdf(paste0(plot.dir, lineage, '_sub', subID, '_k_',  nTopics, '_cellGroup_disease.pdf'), height = 7, width = 10)
print(topic_plots)
dev.off()

## inflammation
topic_plots <- list()
for (i in 1:nTopics) {
  pp <- ggviolin(df, x = 'anno1', y = paste0('k', i),   fill = 'inflammation_status', palette = gca_colors_inflammation,
                 scale = 'width', width=0.8, 
                 trim = T, position = position_dodge(0.7), draw_quantiles = 0.5)+
    labs( x = NULL) +
    theme_pubr()+
    theme(axis.text.y = element_text(size=15),
          axis.text.x = element_text(angle=45, hjust = 1, vjust = 1),
          axis.title = element_text(size=20),
          legend.title = element_blank(),
          legend.text = element_text(size=10),
          legend.position = "bottom",
          plot.title = element_text(size=20, hjust=0.5, face = 'bold'))+
    labs(title = '')+
    facet_wrap(~ biopsy_location, nrow = 2, strip.position = 'left') +
    theme(
      strip.background = element_rect(fill = "white", colour = "white"),
      strip.text = element_text(size = 15, face = 'bold'),
      strip.placement = "outside"
    )
  
  topic_plots[[i]] <- pp
}

pdf(paste0(plot.dir, lineage, '_sub', subID, '_k_',  nTopics, '_cellGroup_inf.pdf'), height = 7, width = 12)
print(topic_plots)
dev.off()


################### topic score vln distribution in different peak groups  ##########################
# use peak annotation and peak-topic. to check peak property vs topic score
df <- as.data.frame(peak_topic_mtx)

# get the peak type information in peakset
peakSet <- proj@peakSet
peak.annot <- peakSet$peakType
names(peak.annot) <- paste(seqnames(peakSet), start(peakSet), end(peakSet), sep = '_')

# match the peaktype information to topic matrix
peakType <- peak.annot[rownames(df)]
names(peakType) <- NULL
df$peakType <- peakType

##  topic by peaktype and location
topic_plots <- list()

for (i in 1:nTopics) {
  pp <- ggviolin(df, x = 'peakType', y = paste0('k', i),  
                 scale = 'width', width=0.8, 
                 trim = T, position = position_dodge(0.7), draw_quantiles = 0.5)+
    labs( x = NULL) +
    theme_pubr()+
    theme(axis.text.y = element_text(size=15),
          axis.text.x = element_text(angle=45, hjust = 1, vjust = 1),
          axis.title = element_text(size=20),
          legend.title = element_blank(),
          legend.text = element_text(size=10),
          legend.position = "bottom",
          plot.title = element_text(size=20, hjust=0.5, face = 'bold'))+
    labs(title = '')
  topic_plots[[i]] <- pp
}

pdf(paste0(plot.dir, lineage, '_sub', subID, '_k_',  nTopics, '_peakGroup_disease.pdf'), height = 4.5, width = 6)
print(topic_plots)
dev.off()



################### topic score correlation with other confounding factors  ##########################
df <- bind_cols(cell_topic_mtx, metadata1[rownames(cell_topic_mtx),])

# Number of regions selected at different lfsr cutoffs:
topic_correlation <- matrix(NA, nrow = nTopics, ncol = 5)
corr_factors <- c("nFrags", "TSSEnrichment", "PromoterRatio", "NucleosomeRatio", "DoubletEnrichment")
colnames(topic_correlation) <- corr_factors
rownames(topic_correlation) <- paste("topic", 1:nrow(topic_correlation))

for(k in 1:nTopics) {
  vector1 <- df[, k]
  topic_correlation[k,] <-
    c(cor(vector1, df$nFrags, method = 'spearman'),
      cor(vector1, df$TSSEnrichment, method = 'spearman'),
      cor(vector1, df$PromoterRatio, method = 'spearman'),
      cor(vector1, df$NucleosomeRatio, method = 'spearman'),
      cor(vector1, df$DoubletEnrichment, method = 'spearman'))
}

topic_correlation

# heatmap plot
library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c("#4393c3", "white", "#d6604d"))
col_fun(seq(-1, 1))

column_ha <- HeatmapAnnotation(foo = anno_text(corr_factors, location = 1, rot = 30, 
                                               gp = gpar(fontsize = 11, fontface='bold')))

pdf(paste0(plot.dir, lineage, '_sub', subID, '_k_',  nTopics, '_topic_confounder_corr.pdf'), height = 6, width = 5)
ht <- Heatmap(topic_correlation, 
        bottom_annotation = column_ha,
        col = col_fun,
        rect_gp = gpar(col = "black", lwd = 2),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", topic_correlation[i, j]), x, y, gp = gpar(fontsize = 10))
        },
        cluster_columns = F, cluster_rows = F,
        show_row_dend = F, show_column_dend = F, 
        show_row_names = T, show_column_names = F,
        row_names_side = "left",
        border = T,
        use_raster = F,
        show_heatmap_legend = FALSE) 
lgd <- Legend(col_fun = col_fun, title = "Spearman correlation score")
draw(ht, heatmap_legend_list = lgd)
dev.off()










                            
                                