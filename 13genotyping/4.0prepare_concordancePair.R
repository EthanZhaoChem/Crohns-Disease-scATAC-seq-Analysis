library(stringr)
metadata <- read.table('~/yuzhao1/work/atac_gca2024/0metadata/meta_Ethan_curated_20240211.csv', header = T, sep = ',')

############## sample sample pair #############
flag <- 0
for (patient in unique(metadata$patient_masked)) {
  samples <- unique(metadata$sample[metadata$patient_masked == patient])
  if(length(samples) < 2){
    next
  }
  if(flag == 0){
    pairs <- combn(samples, 2)
    flag <- 1
    next
  }
  pairs <- cbind(pairs, combn(samples, 2))
  
}

sink("~/yuzhao1/work/final_GCAatac/18glimpse/4.0prepare_samplePair.txt", append = F)
for (i in 1:ncol(pairs)) {
  a <- pairs[1, i]
  b <- pairs[2, i]
  cat(sep=' ', a, b, '\n')
}

sink()

############## patient sample pair #############
flag <- 0
pairs <- c()
for (patient in unique(metadata$patient_masked)) {
  samples <- unique(metadata$sample[metadata$patient_masked == patient])
  if(length(samples) < 2){
    next
  }
  pairs <- c(pairs, paste0(patient, ' ', samples))
}


sink("~/yuzhao1/work/final_GCAatac/18glimpse/4.0prepare_patient_samplePair.txt", append = F)
for (xx in pairs) {
  cat(xx, '\n')
}
sink()


############## complete sample sample pair, chr1 #############
pairs <- combn(metadata$sample, 2)

sink("~/yuzhao1/work/final_GCAatac/18glimpse/4.0prepare_samplePair_all_chr1.txt", append = F)
for (i in 1:ncol(pairs)) {
  a <- pairs[1, i]
  b <- pairs[2, i]
  cat(sep=' ', a, b, '\n')
}

sink()



