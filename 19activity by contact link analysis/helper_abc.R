## 1. link a peak to a abc peak 
## prepare df_abc
# df_abc <- readRDS('~/yuzhao1/work/atac_gca2024/22abc/abc_df_ibd.rds')
# df_abc$peak <- paste0(df_abc$chr, '_', df_abc$start, '_', df_abc$end)
## name_chr_pos_hg38 could be a peak or snp
helper_abcMaxLine <- function(name_chr_pos_hg38, df_abc){
  snp_chrPos <- name_chr_pos_hg38
  gr_gwas <- GRanges(seqnames = snp_chrPos %>% strsplit(split = '_', fixed=T) %>% sapply(.,`[[`,1) %>% as.character(),
                     ranges = IRanges(start = snp_chrPos %>% strsplit(split = '_', fixed=T) %>% sapply(.,`[[`,2) %>% as.numeric(), 
                                      end = snp_chrPos %>% strsplit(split = '_', fixed=T) %>% sapply(.,`[[`,3) %>% as.numeric()))
  gr_abc_peaks <- GRanges(seqnames = df_abc$chr, ranges = IRanges(start = df_abc$start, end = df_abc$end))
  overlaps <- GenomicRanges::findOverlaps(gr_gwas, gr_abc_peaks, ignore.strand = T)
  
  if(length(overlaps@to)==0){
    return(NULL)
  }
  
  df_abc_sub <- df_abc[overlaps@to, ]
  df_abc_sub <- df_abc_sub[which.max(df_abc_sub$ABC.Score),]
  colnames(df_abc_sub) <- paste0('ABC_', colnames(df_abc_sub))
  return(df_abc_sub)
}