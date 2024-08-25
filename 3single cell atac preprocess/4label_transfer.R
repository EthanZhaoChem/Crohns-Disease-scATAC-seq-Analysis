dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)
library(Seurat)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/scripts/gca_markers.R')

addArchRThreads(4)
# read files
setwd('~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2')
proj <- loadArchRProject(path = "~/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2/projdir/")
proj_ref <- loadArchRProject(path = "~/yuzhao1/work/final_GCAatac/0dataset/5kmin_6TSS_DoubletRatio2_filtered1/")

out.dir <- '/project/gca/yuzhao1/work/atac_gca2024/2preprocess/4label_transfer/'
dir.create(out.dir, showWarnings = F)


############################# Section 1: align metadata in proj_ref ################################
patient_meta <- read.csv('~/yuzhao1/work/final_GCAatac/0metadata/patient_metadata_20240211.csv')
sampleIDs_old <- unique(proj_ref$Sample)
tmp3 <- sampleIDs_old
tmp3 %<>% gsub('-', '', .)
biopsy_location <- ifelse(grepl('TI', tmp3), 'TI', 'AC')
patient <- ifelse(grepl('TI', tmp3), str_split(tmp3, 'TI'), str_split(tmp3, 'AC')) %>% sapply(.,`[[`,1)
disease_status <- mapvalues(patient, patient_meta$Helmlsey.ID, patient_meta$Group)
inflammation_status <- ifelse(grepl('TI', tmp3), str_split(tmp3, 'TI'), str_split(tmp3, 'AC')) %>% sapply(.,`[[`,2)
inflammation_status[inflammation_status %in% c('noninf', 'non')] <- 'nonInf'
inflammation_status[inflammation_status %in% c('inf')] <- 'inf'
inflammation_status[inflammation_status %in% c('ad', 'adj', 'adjinf')] <- 'adjInf'
inflammation_status[inflammation_status == '' & disease_status == 'CD'] <- 'nonInf'

df <- data.frame('sample_old' = sampleIDs_old,
                 'patient' = patient,
                 'biopsy_location' = biopsy_location,
                 'disease_status' = disease_status,
                 'inflammation_status' = inflammation_status)

df$patient_masked <- mapvalues(df$patient, patient_meta$Helmlsey.ID, patient_meta$MaskPatientID)
df$sample <- paste0(df$patient_masked, '-', df$biopsy_location, '-', df$disease_status, df$inflammation_status)

metadata <- df
proj_ref$Patient_ID <- mapvalues(proj_ref$Sample, from = metadata$sample_old, metadata$patient)
proj_ref$Sample_aligned <- mapvalues(proj_ref$Sample, from = metadata$sample_old, metadata$sample)
proj_ref$Patient_ID_masked <- mapvalues(proj_ref$Sample, from = metadata$sample_old, metadata$patient_masked)
proj_ref$biopsy_location <- mapvalues(proj_ref$Sample, from = metadata$sample_old, metadata$biopsy_location)
proj_ref$disease_status <- mapvalues(proj_ref$Sample, from = metadata$sample_old, metadata$disease_status)
proj_ref$inflammation_status <- mapvalues(proj_ref$Sample, from = metadata$sample_old, metadata$inflammation_status)
proj_ref$inflammation_status[proj_ref$inflammation_status==''] <- 'Control'

proj_ref$cellNames_aligned <- paste0(proj_ref$Sample_aligned, '#', str_split(proj_ref$cellNames , '#') %>% sapply(.,`[[`,2))

 
################# Section 2: transfer labels from previous analysis ###########
cells_overlapped <- intersect(proj_ref$cellNames_aligned, proj$cellNames)
sprintf('%.2f%% of cells from previous analysis are captured', 100*length(cells_overlapped)/length(proj_ref$cellNames))

matchedID_proj <- match(cells_overlapped, proj$cellNames)
matchedID_ref <- match(cells_overlapped, proj_ref$cellNames_aligned)
proj$transferred_anno1 <- 'NA'
proj$transferred_anno1[matchedID_proj] <- proj_ref$anno1[matchedID_ref]
proj$transferred_anno1_group <- 'NA'
proj$transferred_anno1_group[matchedID_proj] <- proj_ref$anno1_group[matchedID_ref]

saveArchRProject(ArchRProj = proj, load = T)

################# Section 3: plot ###########
png(paste0(out.dir, 'umap_transferred_anno1.png'),width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_2`
df$cluster_name <- proj$transferred_anno1
plot_df_umap_custom(df, show.label = 'name')
dev.off()


png(paste0(out.dir, 'umap_transferred_anno1_group.png'),width = 2600, height = 3000,res = 300)
df <- data.frame(proj@cellColData)
df$embedding1 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_1`
df$embedding2 <- proj@embeddings$LSI_UMAP$df$`IterativeLSI#UMAP_Dimension_2`
df$cluster_name <- proj$transferred_anno1_group
plot_df_umap_custom(df, show.label = 'name')
dev.off()








