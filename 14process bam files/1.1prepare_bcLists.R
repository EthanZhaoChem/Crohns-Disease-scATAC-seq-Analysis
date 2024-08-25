dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)

addArchRThreads(1)
out.dir <- '~/yuzhao1/work/atac_gca2024/20bam/1bam_perBC/barcodes/'
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1/")
bcs <- rownames(proj@cellColData) %>% strsplit(split = '#', fixed=T) %>% sapply(.,`[[`,2)
anno1s <- proj$anno1


for (sample in unique(proj$Sample)) {
  bcs_sub <- bcs[proj$Sample == sample]
  anno1s_sub <- anno1s[proj$Sample == sample]
  df <- data.frame(bcs_sub, anno1s_sub)
  write.table(df, paste0(out.dir, sample, '.txt'), quote = F, col.names = F, row.names = F, sep = '\t')
}
















