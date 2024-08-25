sex_inferred <- read.table('~/yuzhao1/work/final_GCAatac/18glimpse/6.3sex/plink.sexcheck', header = T)

# old metadata
metadata <- read.csv('~/yuzhao1/work/atac_gca2024/0metadata/meta_Ethan_curated_20240211.csv')
sex_metadata <- metadata[, c('patient', 'patient_masked', 'Sex')]
sex_metadata <- sex_metadata[!duplicated(sex_metadata), ]
rownames(sex_metadata) <- sex_metadata$patient_masked

sex_metadata$inferred <- sex_inferred$SNPSEX[match(sex_metadata$patient_masked, sex_inferred$FID)]
sex_metadata$F <- sex_inferred$F[match(sex_metadata$patient_masked, sex_inferred$FID)]
sex_metadata$inferred <- mapvalues(sex_metadata$inferred, c('0', '1', '2'), c('NA', 'Male', 'Female'))

sex_metadata[sex_metadata$Sex != sex_metadata$inferred,]


# curated by inferred sex and new AC sample age information
metadata <- read.table('~/yuzhao1/work/atac_gca2024/0metadata/meta_Ethan_curated_20240311.csv', sep = ',', header = T, row.names = 1)
sex_metadata <- metadata[, c('patient', 'patient_masked', 'Sex')]
sex_metadata <- sex_metadata[!duplicated(sex_metadata), ]
rownames(sex_metadata) <- sex_metadata$patient_masked

sex_metadata$inferred <- sex_inferred$SNPSEX[match(sex_metadata$patient_masked, sex_inferred$FID)]
sex_metadata$F <- sex_inferred$F[match(sex_metadata$patient_masked, sex_inferred$FID)]
sex_metadata$inferred <- mapvalues(sex_metadata$inferred, c('0', '1', '2'), c('NA', 'Male', 'Female'))

sex_metadata[sex_metadata$Sex != sex_metadata$inferred,]
