dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(rasqualTools)
library(ArchR)
library(DESeq2)
library(beeswarm)
celltypes <- readLines('~/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt')
rasqual_dir <- '~/yuzhao1/work/atac_gca2024/19rasqual/5collect_results/5results/allPeaks_noDiseaseCo/'
counts_dir <- '~/yuzhao1/work/atac_gca2024/19rasqual/9boxplot/9.1mtx/'
gt_dir <- '~/yuzhao1/work/atac_gca2024/19rasqual/9boxplot/9.3genotype_clean/'
dir.plot <- '~/yuzhao1/work/atac_gca2024/19rasqual/9boxplot/plots/'
  
# 1. preset
rsID <- 'rs6416647'
ct <- 'TI_Goblet'
peak <- 'chr16_10871644_10872144'
zscore=TRUE
YLIM=NULL

gt_file <- paste0(gt_dir, ct, '_genotypes.tsv')
rasqual_results <- readRDS(paste0(rasqual_dir, ct, '_lead.rds'))
counts <- paste0(counts_dir,ct ,'.vst.count_matrix')


# 2. prepare
rasq <- data.frame(subset(rasqual_results, rasqual_results$Feature_ID == peak))
snp = rasq$rsID
genot = c(paste0(rasq[,5],rasq[,5]), paste0(rasq[,5],rasq[,6]), paste0(rasq[,6],rasq[,6]))
names(genot) = c(0,1,2)

df = t(rbind(read.table(pipe(paste("grep", peak, counts)))[-1],
             read.table(pipe(paste("grep", snp, gt_file)))[-1]))
df <- as.data.frame(df)
colnames(df) = c('peak_counts','genotype')
ngt <- genot[names(genot) %in% unique(df[,2])]
df$genotype <- mapvalues(df$genotype, names(ngt), ngt)
df$genotype <- factor(df$genotype, levels = ngt)

# 3. plot
col <- c('#d53e4f',  '#3288bd', '#41ab5d')
p <- ggplot(df, aes(x=genotype, y=peak_counts, fill=genotype))+
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4, fill='white')+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.125, binwidth = 0.2)+
  scale_fill_manual(values=col) +
  # theme_pubr() + 
  theme(axis.text.x = element_text(color = "black", face = "bold", size = 10,),
        axis.text.y = element_text(color = "black", size = 10, face = "bold"),
        axis.title = element_text(color = "black", size = 10, angle = 0, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),        
        legend.position = "none",
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),)+
  labs(x = rsID, y = "Normalized accessibility", title = "")

pdf(paste0(dir.plot, ct, '-', rsID, '-', snp, '.pdf'), width = 5, height = 5)
print(p)
dev.off()






