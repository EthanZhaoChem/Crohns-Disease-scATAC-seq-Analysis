dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(rasqualTools)
library(ArchR)
snp_dir <- '~/yuzhao1/work/atac_gca2024/19rasqual/1vcf/filtered_snpList/'
celltypes <- readLines('~/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt')
mtx.dir <- '~/yuzhao1/work/atac_gca2024/19rasqual/2.1mtx_peak_patient_raw/'
out.dir <- '~/yuzhao1/work/atac_gca2024/19rasqual/4calculation_all_peaks/4ctSpecific_rasqualInput/'
source('~/yuzhao1/work/atac_gca2024/19rasqual/4helper_covariates.R')

## 1. filter peaks, save binary files
## retained all accessible sites
list1 <- list()
for (temp.celltype in celltypes){
  a <- paste0(mtx.dir, temp.celltype, '.txt')
  tmp.mtx <- read.table(a, sep = '\t')
  nPatients <- ncol(tmp.mtx)
  # no filter
  rowSum_cutoff <- 0
  list1[[temp.celltype]] <- tmp.mtx[rowSums(tmp.mtx) >= rowSum_cutoff, ]
}
saveRasqualMatrices(list1, out.dir, file_suffix = "accessibility")



## 2. Calculate size factors
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1/")
peakset <- proj@peakSet
df_gc <- data.frame(matrix(nrow=length(peakset), ncol=2))
colnames(df_gc) <- c('gene_id', 'percentage_gc_content')
df_gc$gene_id<- paste0(seqnames(peakset), '_', start(peakset), '_', end(peakset))
df_gc$percentage_gc_content <- peakset$GC * 100

list2 <- list()
for (temp.celltype in celltypes){
  tmp.mtx <- list1[[temp.celltype]] 
  tmp.size_factors = rasqualCalculateSampleOffsets(tmp.mtx, df_gc, gc_correct = TRUE)
  list2[[temp.celltype]] <- tmp.size_factors
}
saveRasqualMatrices(list2, out.dir, file_suffix = "size_factors_gc")



## 3. facilitate rasqual input

for (ct in celltypes){
  snp_IDs <- readLines(paste0(snp_dir, ct, '.snp.list'))
  peak_IDs <- rownames(list1[[ct]])
  
  snp_coords <- data.frame(chr = snp_IDs %>% strsplit(split = ':', fixed=T) %>% sapply(.,`[[`,1),
                           pos = snp_IDs %>% strsplit(split = ':', fixed=T) %>% sapply(.,`[[`,2) %>% as.numeric(),
                           snp_id = snp_IDs)
  gene_metadata <- data.frame(gene_id = peak_IDs,
                              chr = peak_IDs %>% strsplit(split = '_', fixed=T) %>% sapply(.,`[[`,1),
                              strand = '*',
                              start = peak_IDs %>% strsplit(split = '_', fixed=T) %>% sapply(.,`[[`,2) %>% as.numeric(),
                              end = peak_IDs %>% strsplit(split = '_', fixed=T) %>% sapply(.,`[[`,3) %>% as.numeric())
  
  snp_counts = countSnpsOverlapingPeaks(gene_metadata, snp_coords, cis_window = 1e4)
  snp_counts <- snp_counts %>% unite("z", chromosome_name, range_start, sep = ":", remove = FALSE) %>% unite("region", z, range_end, sep = "-", remove = FALSE)
  snp_counts <- snp_counts %>% mutate(gene_name = gene_id)
  snp_counts_2 <- snp_counts %>% rstatix::select(gene_id, gene_name, region, cis_snp_count, feature_snp_count, exon_starts, exon_ends)
  
  # this made sure numbers can be read properly by bash
  snp_counts_2 <- lapply(snp_counts_2, function(x) if(is.numeric(x)) format(x, scientific = FALSE) else x)
  snp_counts_2 <- data.frame(snp_counts_2)
  write.table(snp_counts_2, file = paste0(out.dir, ct, '.input.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}



## 4. make covariates
pcs = read.table("~/yuzhao1/work/atac_gca2024/19rasqual/3pca/plink_analysis/final_sample_1KGP_chrALL.eigenvec", header=T)
rownames(pcs) <- pcs$FID
PCs_selected <- paste0('PC', seq(1,4))

sample_metadata <- read.table('~/yuzhao1/work/atac_gca2024/0metadata/meta_Ethan_curated_20240311.csv', sep = ',', header = T, row.names = 1)
sample_metadata <- sample_metadata[, c('patient_masked', 'disease_status', 'Age','Sex')]
patient_metadata <- sample_metadata[!duplicated(sample_metadata$patient_masked), ]
rownames(patient_metadata) <- patient_metadata$patient_masked
patient_metadata$disease_status <- mapvalues(patient_metadata$disease_status, c('CD', 'Control'), c(1, 0))
patient_metadata$Sex <- mapvalues(patient_metadata$Sex, c('Female', 'Male'), c(1, 0))

# # with disease covaraite
# Meta_selected <- c("disease_status", "Age", "Sex" )

# without disease covariate
Meta_selected <- c("Age", "Sex" )


for (ct in celltypes){
  xbin=paste0(out.dir, ct, '.covariates.bin')
  xtxt=paste0(out.dir, ct, '.covariates.txt')
  covs <- rasqualMakeCovariates(list1[[ct]], list2[[ct]])
  patients <- colnames(list1[[ct]])
  covs_pc <- pcs[patients, PCs_selected]
  covs_meta <- patient_metadata[patients, Meta_selected]
  covs_all <- bind_cols(covs, covs_pc, covs_meta)
  write.table(covs_all, xtxt, sep = ',', row.names = T, col.names = T)
  
  covs_all <- as.double(c(as.matrix(covs_all)))
  fxbin=file(xbin,"wb")
  writeBin(covs_all, fxbin)
  close(fxbin)
}

















