dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
setwd('~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2')

##################################################################
# create final metadata for each sample with masked ID
dir.cellRanger <- '/project/spott/yuzhao1/GEO_GCAatac/cellranger_output/'
patient_meta <- read.csv('~/yuzhao1/work/final_GCAatac/0metadata/patient_metadata_20240211.csv')
xx <- read.table('~/yuzhao1/work/final_GCAatac/0metadata/meta_Ethan_curated_20230314.csv',
                       header = T, sep = ',')
tmp3 <- c(xx$Sample_ID, 'HA98AC-ad', 'HA98AC-inf', 'HA102AC-ad', 'HA102AC-inf', 'HA104AC-ad', 'HA104AC-inf')
tmp3 %<>% gsub('-', '', .)
biopsy_location <- ifelse(grepl('TI', tmp3), 'TI', 'AC')
patient <- ifelse(grepl('TI', tmp3), str_split(tmp3, 'TI'), str_split(tmp3, 'AC')) %>% sapply(.,`[[`,1)
disease_status <- mapvalues(patient, patient_meta$Helmlsey.ID, patient_meta$Group)
inflammation_status <- ifelse(grepl('TI', tmp3), str_split(tmp3, 'TI'), str_split(tmp3, 'AC')) %>% sapply(.,`[[`,2)
inflammation_status[inflammation_status %in% c('noninf', 'non')] <- 'nonInf'
inflammation_status[inflammation_status %in% c('inf')] <- 'inf'
inflammation_status[inflammation_status %in% c('ad', 'adj', 'adjinf')] <- 'adjInf'
inflammation_status[inflammation_status == '' & disease_status == 'CD'] <- 'nonInf'

df <- data.frame('patient' = patient,
                 'biopsy_location' = biopsy_location,
                 'disease_status' = disease_status,
                 'inflammation_status' = inflammation_status)

df$patient_masked <- mapvalues(df$patient, patient_meta$Helmlsey.ID, patient_meta$MaskPatientID)
df$sample <- paste0(df$patient_masked, '-', df$biopsy_location, '-', df$disease_status, df$inflammation_status)
df$path <- paste0(dir.cellRanger, df$sample, '/outs/fragments.tsv.gz')
sum(file.exists(df$path)) == nrow(df)

##################################################################
# add other metadata
for (tmp.name in c( "Age", "Race", "Sex", "History", "Treatment")) {
  df[,tmp.name] <- mapvalues(df$patient, patient_meta$Helmlsey.ID, patient_meta[,tmp.name])
}
# write.table(df, sep = ',', file = '~/yuzhao1/work/atac_gca2024/0metadata/meta_Ethan_curated_20240211.csv')

##################################################################
metadata <- read.table('~/yuzhao1/work/atac_gca2024/0metadata/meta_Ethan_curated_20240211.csv',
                       header = T, sep = ',')

# set global parameters
filterRatio = 2.0
minFrags <- 5000
minTSS <- 6
GenomeSet <- 'hg38'
pathToMacs2 <- '/home/yuzhao1/.local/bin/macs2'
addArchRGenome(GenomeSet)
addArchRThreads(48)

#prepare framents' paths
inputFiles <- paste0(metadata$path)
names(inputFiles) <- metadata$sample

#prepare arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  outputNames = names(inputFiles),
  minTSS = minTSS,
  minFrags = minFrags,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  verbose = T,
  force = T
)


# get arrow files' names
ArrowFiles <- paste0(names(inputFiles), '.arrow')

# infer doublets for each arrow file
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1,
  nTrials = 10,
  dimsToUse = 2:50,
  force = T
)

# create, clean and reload archr project
proj <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "projdir", copyArrows = F)
proj <- filterDoublets(ArchRProj = proj, filterRatio = filterRatio)
saveArchRProject(ArchRProj = proj, load = T)




