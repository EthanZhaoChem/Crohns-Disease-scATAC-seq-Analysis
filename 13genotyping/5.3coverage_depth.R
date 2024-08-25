# consider removing 3 samples: P3(HA16-AC-CDnonInf), P7(HA26-TI-CDnonInf), P9(HA30-AC-CDnonInf)

library(stringr)
genome_length <- 3099750718
metadata <- read.table('~/yuzhao1/work/atac_gca2024/0metadata/meta_Ethan_curated_20240211.csv', header = T, sep = ',')
patients <- unique(metadata$patient_masked)
df <- data.frame(matrix(0, nrow=length(patients), ncol=2))
rownames(df) <- mixedsort(patients)
colnames(df) <- c('patient', 'depth')

for (patient in patients) {
  xx <- readLines(paste0('~/yuzhao1/work/final_GCAatac/18glimpse/5.3coverage_depth/', patient, '.bam'))
  xx %<>% as.numeric()
  df[patient, 1] <- patient
  df[patient, 2] <- xx
}

write.csv(df, '~/yuzhao1/work/final_GCAatac/18glimpse/5.3coverage_depth/00summary_patient_bam.csv')

