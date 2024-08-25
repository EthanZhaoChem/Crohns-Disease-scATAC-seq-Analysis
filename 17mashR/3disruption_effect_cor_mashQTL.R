lfsr_threshold <- 0.05

dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP150.GRCh38) 
library(BSgenome.Hsapiens.UCSC.hg38)     
library(BSgenome)
library(plyr)
library(dplyr)
library(stringr)
library(tibble)
library(Seurat)
library(ArchR)
source('~/yuzhao1/work/atac_gca2024/scripts/gca_colors.R')
source('~/yuzhao1/work/final_GCArna/scripts/gca_markers.R')
out.dir <- '~/yuzhao1/work/atac_gca2024/26mash3/3disruption_effect_cor_mashQTL/'

# helper function
check_snp_in_peak <- function(peak, snp){
  flag <- 0
  a.chr <- snp %>% strsplit(split = ':', fixed=T) %>% sapply(.,`[[`,1)
  a.pos <- snp %>% strsplit(split = ':', fixed=T) %>% sapply(.,`[[`,2) %>% as.numeric()
  b.chr <- peak %>% strsplit(split = '_', fixed=T) %>% sapply(.,`[[`,1)
  b.start <- peak %>% strsplit(split = '_', fixed=T) %>% sapply(.,`[[`,2) %>% as.numeric()
  b.end <- peak %>% strsplit(split = '_', fixed=T) %>% sapply(.,`[[`,3) %>% as.numeric()
  if(a.chr == b.chr & a.pos >= b.start & a.pos <= b.end){
    flag <- 1
  }
  return(flag)
}

# 1.read files
celltypes <- readLines('~/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt')

# in this df, each row is a tf-genesymbol pair
motif.breaks.all_unique <- readRDS('~/yuzhao1/work/atac_gca2024/26motif_disruption3/3results_summary/motif.breaks.all_unique.rds')
df_lfsr <- readRDS('~/yuzhao1/work/atac_gca2024/26mash3/2mash/df_lfsr.rds')
qtls_ctList <- list()
for (ct in celltypes) {
  qtls_ctList[[ct]] <- rownames(df_lfsr)[df_lfsr[, ct] <= lfsr_threshold]
}

# read effect sizes
prep.dir <- '~/yuzhao1/work/atac_gca2024/26mash3/1select_qtl/'
df_effect <- readRDS(paste0(prep.dir, 'df_effect_query.rds'))
df_effect$Feature_ID <- df_effect %>% rownames() %>% strsplit(split = '@', fixed=T) %>% sapply(.,`[[`,1) %>% unlist()
df_effect$rsID <- df_effect %>% rownames() %>% strsplit(split = '@', fixed=T) %>% sapply(.,`[[`,2) %>% unlist()

# 2.filter TFs
tfs_disrupted_all <- motif.breaks.all_unique$geneSymbol
tfs_disrupted_frequency <- table(tfs_disrupted_all) %>% sort(decreasing = T) %>% data.frame() %>% set_colnames(., c('TF', 'Freq'))
df_TF_regulator <- readRDS('/project/gca/yuzhao1/work/atac_gca2024/5TF/output/union/regulator_df_anno1_inflammation_status_location_cor0.4_delta0.25.rds')
TF_positive <- df_TF_regulator[df_TF_regulator$TFRegulator=='Positive', 'cisbp_matchName']
tfs_disrupted_all <- as.character(unique(tfs_disrupted_frequency$TF))
tfs_disrupted_all <- intersect(tfs_disrupted_all, TF_positive) # choose not to filter for positive regulators


# 3: calculate correlation of motif disruption score versus effect size 
min_QTLs <- 5

for (ct in celltypes) {
  dir.create(paste0(out.dir, ct), showWarnings = F)
  for (tf in tfs_disrupted_all) {
    # filter for ct
    df <- df_effect[qtls_ctList[[ct]], c(ct, 'Feature_ID', 'rsID')]
    colnames(df) <- c('Effect_size', 'Feature_ID', 'rsID')

    # filter for tf
    df_disrupt_sub <- motif.breaks.all_unique[motif.breaks.all_unique$geneSymbol == tf, ]
    df <- df[df$rsID %in% df_disrupt_sub$SNP_id, ] # QTLs for this TF
    
    # consider it not significant anyway, if fewer than min QTLs
    if(nrow(df) < min_QTLs){
      tmp.result <- 'Too few QTLs'
      saveRDS(tmp.result, paste0(out.dir, ct, '/', tf, '.rds'))
      next
    }
    
    # filter for snps in capeaks
    df$rsID_in_FeatureID <- 0
    for (i in 1:nrow(df)) {
      peak <- df$Feature_ID[[i]]
      snp <- df$rsID[[i]]
      if(check_snp_in_peak(peak, snp)){
        df$rsID_in_FeatureID[[i]] <- 1
      }
    }
    df <- df[df$rsID_in_FeatureID == 1,]
    
    # consider it not significant anyway, if fewer than min QTLs
    if(nrow(df) < min_QTLs){
      tmp.result <- 'Too few QTLs'
      saveRDS(tmp.result, paste0(out.dir, ct, '/', tf, '.rds'))
      next
    }
    
    # calculate cor score
    df$disruption_score_alleleDiff <- mapvalues(df$rsID, from=df_disrupt_sub$SNP_id, to=df_disrupt_sub$alleleDiff, warn_missing = F) %>% as.numeric()
    tmp.result <- cor.test(df$disruption_score_alleleDiff, df$Effect_size, method = 'spearman')
    saveRDS(tmp.result, paste0(out.dir, ct, '/', tf, '.rds'))
  }
}

































