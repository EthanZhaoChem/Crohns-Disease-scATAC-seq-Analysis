library(mapgen)
library(tidyverse)
library(ggplot2)
library(data.table)
source('~/yuzhao1/scripts/plot.R')
out.dir <- '~/yuzhao1/work/final_GCAatac/8susie/3compare_1kg_ukbb_results/'

df_1kg <- read.csv('~/yuzhao1/work/final_GCAatac/8susie/archive_1kg/results/cd_finemapping_unifprior_95loci_L10.csv', row.names = 1)
df_ukb <- read.csv('~/yuzhao1/work/final_GCAatac/8susie/results/cd_finemapping_unifprior_L10.csv', row.names = 1)


## 1. plot all snps 
snps_all <- union(df_1kg$snp, df_ukb$snp)
df <- data.frame(snps = snps_all,
                 pip_1kg = 0,
                 pip_ukb = 0)
rownames(df) <- df$snps

df[df_1kg$snp, 'pip_1kg'] <- df_1kg$susie_pip
df[df_ukb$snp, 'pip_ukb'] <- df_ukb$susie_pip

# plot
p <- ggplot(df, aes(x=pip_1kg, y=pip_ukb)) + 
  # geom_smooth(method="lm", col="black", formula = 'y ~ x', fill = '#a6cee3', size=0.5) + 
  geom_point(size=0.5, color='#0570b0') + 
  theme_pubr() + # Use ggpubr theme
  theme(axis.title = element_text(size = 8, face = 'bold'),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "bottom") +
  labs(title = 'pip', x = "1kg", y = "ukb")


pdf(paste0(out.dir, 'ukb_1kg_pip_all.pdf'), height = 3, width = 3)
print(p)
dev.off()




## 2. plot shared pip
snps_shared <- intersect(df_1kg$snp, df_ukb$snp)
df <- data.frame(snps = snps_shared,
                 pip_1kg = 0,
                 pip_ukb = 0)
rownames(df) <- df$snps

df[snps_shared, 'pip_1kg'] <- df_1kg[snps_shared, 'susie_pip']
df[snps_shared, 'pip_ukb'] <- df_ukb[snps_shared, 'susie_pip']

# plot
p <- ggplot(df, aes(x=pip_1kg, y=pip_ukb)) + 
  # geom_smooth(method="lm", col="black", formula = 'y ~ x', fill = '#a6cee3', size=0.5) + 
  geom_point(size=0.5, color='#0570b0') + 
  theme_pubr() + # Use ggpubr theme
  theme(axis.title = element_text(size = 8, face = 'bold'),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size=8),
        legend.position = "bottom") +
  labs(title = 'pip', x = "1kg", y = "ukb")


pdf(paste0(out.dir, 'ukb_1kg_pip_shared.pdf'), height = 3, width = 3)
print(p)
dev.off()



