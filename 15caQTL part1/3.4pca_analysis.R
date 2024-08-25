dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(rasqualTools)
library(ggplot2)

# pca vector for all sample and 1kg
pcs = read.table("~/yuzhao1/work/atac_gca2024/19rasqual/3pca/plink_analysis/final_sample_1KGP_chrALL.eigenvec", header=T)

# explained variance
var = readLines("~/yuzhao1/work/atac_gca2024/19rasqual/3pca/plink_analysis/final_sample_1KGP_chrALL.eigenval")
var_expl = paste0( "PC", 1:20, " (",round((as.numeric(var)/sum(as.numeric(var)))*100,1), "%)")

# 1kg metadata
meta_1kg1 <- read.table('~/yuzhao1/work/atac_gca2024/19rasqual/3pca/1kg_metadata/1kg_phase3.tsv', header = T, sep='\t')
meta_1kg2 <- read.table('~/yuzhao1/work/atac_gca2024/19rasqual/3pca/1kg_metadata/1kg_grch38.tsv', header = T, sep='\t')
meta_1kg <- rbind(meta_1kg1, meta_1kg2)
meta_1kg <- meta_1kg[!duplicated(meta_1kg$Sample.name), ]
gca_samples <- readLines('~/yuzhao1/work/final_GCAatac/18glimpse/3.1patientIDs_all')

# build df
df <- pcs
df$sample <- pcs$FID
df$population <- meta_1kg$Superpopulation.code[match(df$sample, meta_1kg$Sample.name)]
df$population[df$sample %in% gca_samples] <- 'gca'
df$population[is.na(df$populatio)] <- 'NA'

# plot
df$population <- as.factor(df$population)
manual_colors <- rainbow(length(levels(df$population)))
df_gca = subset(df, grepl("gca", df$population))

pdf("~/yuzhao1/work/atac_gca2024/19rasqual/3pca/plots/pc1_2.pdf", height = 5, width = 5)
ggplot(df, aes(x = PC1, y = PC2, color = population)) +
  geom_point() +
  scale_color_manual(values = manual_colors) +
  labs(x = var_expl[1], y = var_expl[2], color = "Population") +
  theme_pubr()

ggplot(df_gca, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = factor(population)), size = 3) +  # Assuming 'bg' is a factor for colors
  geom_text_repel(aes(label = sample), max.overlaps = 200) +  # Adjust text position as needed
  scale_color_manual(values = 'orange') +  # Define your color values
  labs(x = var_expl[1], y = var_expl[2]) +
  theme_pubr()
dev.off()

pdf("~/yuzhao1/work/atac_gca2024/19rasqual/3pca/plots/pc3_4.pdf", height = 5, width = 5)
ggplot(df, aes(x = PC3, y = PC4, color = population)) +
  geom_point() +
  scale_color_manual(values = manual_colors) +
  labs(x = var_expl[3], y = var_expl[4], color = "Population") +
  theme_pubr()

ggplot(df_gca, aes(x = PC3, y = PC4)) +
  geom_point(aes(color = factor(population)), size = 3) +  # Assuming 'bg' is a factor for colors
  geom_text_repel(aes(label = sample), max.overlaps = 200) +  # Adjust text position as needed
  scale_color_manual(values = 'orange') +  # Define your color values
  labs(x = var_expl[3], y = var_expl[4]) +
  theme_pubr()
dev.off()


pdf("~/yuzhao1/work/atac_gca2024/19rasqual/3pca/plots/pc5_6.pdf", height = 5, width = 5)
ggplot(df, aes(x = PC5, y = PC6, color = population)) +
  geom_point() +
  scale_color_manual(values = manual_colors) +
  labs(x = var_expl[5], y = var_expl[6], color = "Population") +
  theme_pubr()

ggplot(df_gca, aes(x = PC5, y = PC6)) +
  geom_point(aes(color = factor(population)), size = 3) +  # Assuming 'bg' is a factor for colors
  geom_text_repel(aes(label = sample), max.overlaps = 200) +  # Adjust text position as needed
  scale_color_manual(values = 'orange') +  # Define your color values
  labs(x = var_expl[5], y = var_expl[6]) +
  theme_pubr()
dev.off()























