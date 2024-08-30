library(bigsnpr)
library(mapgen) # mapgen_0.5.9
library(susieR) # susieR_0.12.35
library(data.table)
library(dplyr)
library(stringr)
gwas_folder <- '~/yuzhao1/work/final_GCAatac/0gwas/b151_hg19/'
results_folder <- '~/yuzhao1/work/final_GCAatac/8susie/archive_1kg/results/'
## 1. download and change dataformat
# download_1000G('~/yuzhao1/work/final_GCAatac/8susie/1ref_genotype', overwrite = FALSE, delete_zip = TRUE)
# snp_readBed('~/yuzhao1/work/final_GCAatac/8susie/1ref_genotype/1000G_phase3_common_norel.bed')
# snp_readBed('~/yuzhao1/resource/magma/g1000_eur.bed')


## 2. attach snp file (eur 1000g)
bigSNP <- snp_attach(rdsfile = '~/yuzhao1/resource/magma/g1000_eur.rds')
# data('Euro_LD_Chunks', package='mapgen') # version: mapgen_0.3.12
# saveRDS(LD_Blocks, '~/yuzhao1/work/final_GCAatac/8susie/LD_blocks_mapgen/LD_Blocks.rds') # save in case the package was deprecated
LD_Blocks <- readRDS('~/yuzhao1/work/final_GCAatac/8susie/LD_blocks_mapgen/LD_Blocks.rds') 

## 3. read gwas file
gwas_raw<- as.data.frame(
  fread(
    paste0(gwas_folder, "cd_build37_40266_20161107_hg19_b151_gwas_annotated_notFilterFor1KG.txt")
  )
)

gwas <- process_gwas_sumstats(sumstats = gwas_raw,
                              bigSNP=bigSNP,
                              LD_Blocks=LD_Blocks)
saveRDS(gwas, paste0(results_folder, 'gwas_processed.rds'))


  
## 4. uniform prior
n = 40266
sig.loci <- gwas %>%
  group_by(locus) %>%
  summarise(max_mlogP = max(-log10(pval))) %>%
  filter(max_mlogP > -log10(5e-8)) %>% pull(locus)

gwas.sumstats.sigloci <- gwas[gwas$locus %in% sig.loci, ]
cat(length(unique(gwas.sumstats.sigloci$locus)), "GWAS significant loci. \n")

susie_finemap_L10 <- run_finemapping(gwas.sumstats.sigloci, 
                                    bigSNP, 
                                    n = n,
                                    priortype = 'uniform',
                                    L = 10)

## save results 
saveRDS(susie_finemap_L10, paste0(results_folder, 'cd_finemapping_unifprior_95loci_L10.rds'))


## 5, plot pip
susie_finemap_L10 <- readRDS(paste0(results_folder, 'cd_finemapping_unifprior_95loci_L10.rds'))
for (i in 1:length(susie_finemap_L10)) {
  susie_locus_id <- names(susie_finemap_L10)[[i]]
  susie_locus <- susie_finemap_L10[[i]]
  png(paste0('~/yuzhao1/work/final_GCAatac/8susie/1locus_plot_L10/', susie_locus_id, '.png'), height = 2000, width = 2000, res = 300)
  susie_plot(susie_locus, y='PIP', b=NULL)
  dev.off()
}




