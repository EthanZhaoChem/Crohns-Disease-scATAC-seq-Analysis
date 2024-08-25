# supposed to calculate all snps
dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP150.GRCh38) 
library(BSgenome.Hsapiens.UCSC.hg38)     
library(BSgenome)
library(plyr)
library(stringr)

# 1. preset
bedDir <- '~/yuzhao1/work/atac_gca2024/26motif_disruption3/1rsID_beds/'
resultsDir <- '/project/gca/yuzhao1/work/atac_gca2024/26motif_disruption3/2results_allHumanTF/'
out.dir <- '/project/gca/yuzhao1/work/atac_gca2024/26motif_disruption3/4plot/'
df_bed <- readRDS(paste0(bedDir, '00rsID_chunkID_reference.rds'))

# 2. given a rsID to plot
rsID_selected <- c('chr8:46931279:A:G', 'chr13:30697672:A:G', 'chr20:63404353:G:T', 'chr10:3095817:C:T', 'chr17:73720436:A:G',
                   'chr3:16078426:G:T', 'chr8:134557115:C:G', 'chr9:77603010:T:A',
                   'chr2:170377725:C:G', 'chr3:16078426:G:T', 'chr12:52454608:T:C',
                   'chr12:33573499:A:G', 'chr15:73390367:G:T',
                   'chr2:38148629:C:T', 'chr10:28458038:G:A', 'chr11:74756368:T:G', 'chr21:37195601:C:T',
                   'chr3:42490447:G:A', 'chr3:53218322:G:A',
                   'chr1:94543099:G:A', 'chr1:226610032:G:A', 'chr3:42490447:G:A', 
                   'chr14:54674574:C:T')

plotnames <- c('AC_Goblet_KLF4', 'AC_Goblet_KLF4', 'TI_Goblet_KLF4', 'CD8T_RUNX3', 'CD8T_RUNX3',
               'Colonocyte_CDX2', 'Colonocyte_CDX2', 'Colonocyte_CDX2',
               'Early_Colonocyte_CDX2', 'Early_Colonocyte_CDX2', 'Early_Colonocyte_CDX2',
               'Early_Enterocyte_HNF4A', 'Early_Enterocyte_HNF4A',
               'Enterocyte_HNF4A', 'Enterocyte_HNF4A', 'Enterocyte_HNF4A', 'Enterocyte_HNF4A',
               'Macrophage_SPI1', 'Macrophage_SPI1',
               'Plasma_IRF4', 'Plasma_IRF4', 'Plasma_IRF4', 
               'Plasma_TCF3')

rsID_selected <- c('chr12:52454608:T:C')

plotnames <- c('Early_Colonocyte_CDX2')

for(i in 1:length(rsID_selected)){
  rsID_plot <- rsID_selected[[i]]
  plot.name <- plotnames[[i]]
  # source back to its results chunk and read
  chunkID <- df_bed[df_bed$ID==rsID_plot, 'chunkID']
  results <- readRDS(paste0(resultsDir, chunkID, '.rds'))
  
  plot.height <- sum(results$SNP_id==rsID_plot)
  pdf(paste0(out.dir, plot.name, '-', rsID_plot, '.pdf'), width = 12, height = plot.height, pointsize = 1)
  plotMB(results = results, rsid = rsID_plot, effect = c("strong", "weak"))
  dev.off()
}














