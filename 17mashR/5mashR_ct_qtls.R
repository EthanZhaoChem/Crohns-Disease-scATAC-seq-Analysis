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
source('~/yuzhao1/scripts/plot.R')
"%&%" <- function(a, b) paste0(a, b)

##################  ################ ##################  ################
## presets
celltypes <- readLines('~/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt')
dir.plots <- '~/yuzhao1/work/atac_gca2024/26mash3/5mashR_ct_qtls/'
df_lfsr <- readRDS('~/yuzhao1/work/atac_gca2024/26mash3/2mash/df_lfsr.rds')

##################  ################ ##################  ################
# save cell type specific qtls after mash
qtls_ctList <- list()
for (ct in celltypes) {
  qtls_ctList[[ct]] <- rownames(df_lfsr)[df_lfsr[, ct] <= 0.05]
}
qtls_ctList_lengths <- lapply(qtls_ctList, length) %>% unlist()

qtls_all <- unique(unlist(qtls_ctList))
qtls_counts <- setNames(rep(0, length(qtls_all)), qtls_all)
for (ct in celltypes) {
  tmp_qtls <- qtls_ctList[[ct]]
  qtls_counts[tmp_qtls] <- qtls_counts[tmp_qtls] + 1
}


##################  ################ ##################  ################
# plot proportion

df <- data.frame(matrix(0, nrow=4, ncol=length(celltypes)))
colnames(df) <- celltypes
rownames(df) <- c('1~5', '6~10', '11~15', '16~19')
for(ct in celltypes){
  tmp_qtls <- qtls_ctList[[ct]]
  tmp.counts <-  qtls_counts[tmp_qtls]
  df[1, ct] <- sum(tmp.counts>=1 & tmp.counts<=5)/length(tmp.counts)
  df[2, ct] <- sum(tmp.counts>=6 & tmp.counts<=10)/length(tmp.counts)
  df[3, ct] <- sum(tmp.counts>=11 & tmp.counts<=15)/length(tmp.counts)
  df[4, ct] <- sum(tmp.counts>=16 & tmp.counts<=19)/length(tmp.counts)
}


##################  ################ ##################  ################
# plot exact number
df <- data.frame(matrix(0, nrow=4, ncol=length(celltypes)))
colnames(df) <- celltypes
rownames(df) <- c('1~5', '6~10', '11~15', '16~19')
for(ct in celltypes){
  tmp_qtls <- qtls_ctList[[ct]]
  tmp.counts <-  qtls_counts[tmp_qtls]
  df[1, ct] <- sum(tmp.counts>=1 & tmp.counts<=5)
  df[2, ct] <- sum(tmp.counts>=6 & tmp.counts<=10)
  df[3, ct] <- sum(tmp.counts>=11 & tmp.counts<=15)
  df[4, ct] <- sum(tmp.counts>=16 & tmp.counts<=19)
}

df_long <- data.frame(matrix(nrow=nrow(df) * ncol(df), ncol=3))
colnames(df_long) <- c('celltype', 'group', 'value')
count <- 0
for (i in 1:nrow(df)) {
  for (j in 1:ncol(df)) {
    count <- count + 1
    df_long[count, 1] <- colnames(df)[[j]]
    df_long[count, 2] <- rownames(df)[[i]]
    df_long[count, 3] <- df[i, j]
  }
}

df_long$celltype <- factor(df_long$celltype, levels = names(sort(qtls_ctList_lengths, decreasing = T)))
df_long$group <- factor(df_long$group, levels = c("1~5", "6~10", "11~15", "16~19"))

p <- ggplot(df_long, aes(x = celltype, y = value, fill = group)) +
  geom_bar(position = "stack", stat = "identity", width = 0.72)+
  scale_fill_manual(values = c("#d6604d", "#4DBBD5FF", "#00A087FF", "#8073ac"))+
  theme_pubr()+
  theme(axis.text.y = element_text(size=15),
        axis.text.x = element_text(angle=45, hjust = 1, vjust = 1),
        axis.title = element_text(size=20),
        legend.title = element_text(size=13),
        legend.text = element_text(size=13),
        legend.position = "bottom",
        plot.title = element_text(size=20, hjust=0.5, face = 'bold'))+
  labs(title = "", y = "# caqtls", x = "", fill = "")+
  guides(fill=guide_legend(title="# shared celltypes"))

pdf(paste0(dir.plots, 'caqtls_shared celltype.pdf'), width = 7, height = 6, pointsize = 1)
print(p)
dev.off()


















