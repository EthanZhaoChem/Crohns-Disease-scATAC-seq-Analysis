library(stringr)
sample_pairs_all <- readLines('~/yuzhao1/work/final_GCAatac/18glimpse/4.0prepare_patient_samplePair.txt')
df <- data.frame(matrix(0, nrow = length(sample_pairs_all), ncol = 3))
rownames(df) <- sample_pairs_all
colnames(df) <- c('sample1', 'sample2', 'NRC')
for (sample_pair in sample_pairs_all) {
  sample1 <- strsplit(sample_pair, ' ')[[1]][[1]]
  sample2 <- strsplit(sample_pair, ' ')[[1]][[2]]
  resultPath <- paste0("~/yuzhao1/work/final_GCAatac/18glimpse/4.2patient_sample_concordance/", sample1, '_', sample2, ".rds")
  result <- readRDS(resultPath)
  df[sample_pair,] <- c(sample1, sample2, result$nonreference_concordance)
}

write.csv(df, '~/yuzhao1/work/final_GCAatac/18glimpse/4.2patient_sample_concordance/00summary_NRC.csv')











