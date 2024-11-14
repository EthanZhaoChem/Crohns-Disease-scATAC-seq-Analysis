dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(rlang)
library(ggplot2)
library(cowplot)
library(Rtsne)
library(fastTopics)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/work/atac_gca2024/scripts/gca_colors.R')
source('~/yuzhao1/work/final_GCArna/scripts/gca_markers.R')

setwd('/project/gca/yuzhao1_topic/')
out.dir <- '~/yuzhao1/work/atac_gca2024/0manu/plots/4fasttopic_all_structure/'
set.seed(111)

## 1. set the order of clusters
clusters_epithelial <- c("TI_Stem", "AC_Stem", "Early_Enterocyte", "Enterocyte",
                         "Early_Colonocyte", "Colonocyte", "TI_Goblet", "AC_Goblet",
                         "BEST4", "Tuft", "EEC", "Paneth")
clusters_immune <- c("CD4T", "CD8T", "gdT", "NK", 
                     "ILCs", "GC_B", "NaiveB", "MemoryB", "Plasma", 
                     "Mast", "Macrophage", "DC", "Neutrophil")
clusters_stromal <- c("Fibroblast", "Pericyte", "Endothelium", "Glial")

cluster_names <- c("TI_Stem", "AC_Stem", "Early_Enterocyte", "Enterocyte",
                   "Early_Colonocyte", "Colonocyte", "TI_Goblet", "AC_Goblet",
                   "BEST4", "Tuft", "EEC", "Paneth",
                   "CD4T", "CD8T", "gdT", "NK", 
                   "ILCs", "GC_B", "NaiveB", "MemoryB", "Plasma", 
                   "Mast", "Macrophage", "DC", "Neutrophil",
                   "Fibroblast", "Pericyte", "Endothelium", "Glial")
topics <- c("k26", "k30", "k24", "k6", "k4", "k38", "k27", "k29", "k34", "k5", "k37", "k43", "k8", "k23", "k17", "k44", "k36", "k13", "k19", "k10", "k45", "k14", "k32", "k35", "k1", "k42", "k15", "k18", "k25", "k28", "k39", "k41", "k12", "k40", "k7", "k3", "k22", "k31", "k33", "k16", "k9", "k2", "k21", "k11", "k20")

## 2. choose what to color

# covariate <- "biopsy_location"
# colors_for_plot <- gca_colors_location


covariate <- "inflammation_status"
colors_for_plot <- gca_colors_inflammation

# Do this to remove the covariate from the plot:
# colors_for_plot <- rep("black",4)


## 3. read data and set order of topics
meta_data <- readRDS("metadata.rds")
meta_data <- transform(meta_data,
                       biopsy_location     = factor(biopsy_location,
                                                    c('TI', 'AC')),
                       inflammation_status = factor(inflammation_status,
                                                    c("Control","nonInf",
                                                      "adjInf","inf")),
                       anno1               = factor(anno1,cluster_names))

fit <- readRDS("fit_union_sub100_k45_converged.rds")
fit <- poisson2multinom(fit)
L         <- fit$L
clusters  <- meta_data$anno1
n         <- nrow(fit$L)
lineage <- rep("",n)
lineage[is.element(clusters,clusters_epithelial)] <- "epithelial"
lineage[is.element(clusters,clusters_immune)] <- "immune"
lineage[is.element(clusters,clusters_stromal)] <- "stromal"
lineage <- factor(lineage)
rows      <- order(clusters,meta_data[[covariate]])
meta_data <- meta_data[rows,]
clusters  <- clusters[rows]
lineage   <- lineage[rows]
L    <- L[rows,]
n    <- nrow(fit$L)
rows <- NULL
for (i in levels(clusters))
  rows <- c(rows,sort(sample(which(clusters == i), 150)))
meta_data <- meta_data[rows,]
clusters  <- clusters[rows]
lineage   <- lineage[rows]
L         <- L[rows,topics]
n <- nrow(L)
k <- ncol(L)
gap <- 30
x <- rep(0,n)
total <- 0
for (i in levels(clusters)) {
  N <- sum(clusters == i)
  x[clusters == i] <- total + 1:N
  total <- total + gap + N
}
pdat <- data.frame(x          = x,
                   lineage    = lineage,
                   biopsy_location = meta_data$biopsy_location,
                   inflammation_status = meta_data$inflammation_status,
                   topic      = factor(rep(topics,each = n),topics),
                   membership = as.vector(L))
pos <- tapply(x,clusters,mean)

# remove k for plot
pdat$topic <- gsub('k', '', pdat$topic)
topics_numberOnly <- gsub('k', '', topics)
pdat$topic <- factor(pdat$topic, topics_numberOnly)

p <- ggplot(pdat,aes(x = x,y = membership,
                     color = .data[[covariate]],
                     fill = .data[[covariate]])) +
  facet_grid(rows = vars(topic), scales = "free") +
  geom_col() +
  scale_x_continuous(limits = c(0,max(x)+1),breaks = pos,
                     labels = levels(clusters)) +
  scale_color_manual(values = colors_for_plot) +
  scale_fill_manual(values = colors_for_plot) +
  theme_cowplot(font_size = 7) +
  theme(strip.background = element_blank(),
        strip.text.y = element_text(size = 7,angle = 0),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "")

# pdf(paste0(out.dir, 'loc_colored.pdf'), height = 9, width = 6)
# print(p)
# dev.off()


pdf(paste0(out.dir, 'inf_colored.pdf'), height = 9, width = 6)
print(p)
dev.off()





















