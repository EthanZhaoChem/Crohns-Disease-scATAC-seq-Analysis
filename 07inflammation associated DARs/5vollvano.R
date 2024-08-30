

####################### given a peak, give a list of TFs ####################
celltype <- 'ECorCC'
contrast <- 'TI_remission_Control'
df <- statistics_allContrasts[[contrast]][[celltype]]
df$mlogp <- -log10(df$P.Value)

# add region information
peaks_test <- peaks_all[as.numeric(rownames(df))]
df$region <- paste0(seqnames(peaks_test), '_', start(peaks_test), '_', end(peaks_test))

# interactively select a motif to view
View(df[df$region %in% b,])

# rowID reference table to peaks
ref <- readRDS(paste0(proj@projectMetadata$outputDirectory, "/Annotations/Motif-Matches-In-Peaks.rds"))
peak_motif_table <- ref@assays@data$matches
rowranges <- ref@rowRanges
rownames(peak_motif_table) <- paste0(seqnames(rowranges ), '_', start(rowranges), '_', end(rowranges))

# for a peak, check what TFs are there
tmp.peak <- 'chr13_48121064_48121564'
colnames(peak_motif_table)[peak_motif_table[tmp.peak,]]

# for a  peak, give its corresponding information
tmp.peak <- c
peaks_vector_all <- paste0(seqnames(proj@peakSet), '_', start(proj@peakSet), '_', end(proj@peakSet))
xx <- proj@peakSet[peaks_vector_all %in% tmp.peak,]
xx[which(xx$nearestGene=='STAT4'), ]

####################### volcano plot ####################################
celltype <- 'ECorCC'
contrast <- 'TI_remission_Control'
logFC_threshold <- 1
mlogp_threshold <- 2
plot.title <- paste0('DAR_volcano_', celltype, '_', contrast, '_FC_', logFC_threshold, '_mlogP_', mlogp_threshold)

# df
df <- statistics_allContrasts[[contrast]][[celltype]]
df$mlogp <- -log10(df$P.Value)
df <- df %>%
  mutate(Significant = ifelse(abs(logFC) > logFC_threshold & mlogp > mlogp_threshold, "Significant", "Not Significant"))

p <- ggscatter(df, x = 'logFC', y = 'mlogp', color = 'Significant', size = 2) +
  theme_pubr() +
  scale_color_manual(values = c("Not Significant" = "#bdbdbd", "Significant" = "#fd8d3c")) +
  theme(axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")+
  # geom_text_repel(aes(label=label),
  #                 nudge_y = 0, max.overlaps = 1000, min.segment.length = 0.01,
  #                 size = 2, force = 20) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = mlogp_threshold, linetype = "dashed") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed") +
  guides(color=guide_legend(title=""))+
  labs(title=paste0(celltype, ':', contrast),
       x = "Log2 Fold Change",
       y = "-log10(P-value)",
       color = "Significance") +
  theme(legend.position = "bottom")

# check number
sum(df$logFC > logFC_threshold & df$mlogp > mlogp_threshold)
sum(df$logFC < -logFC_threshold & df$mlogp > mlogp_threshold)

# plot
pdf(paste0(plot.dir, plot.title, '.pdf'), width = 5, height = 5)
print(p)
dev.off()

png(paste0(plot.dir, plot.title, '.png'), width = 2000, height = 2000, res = 300)
print(p)
dev.off()


