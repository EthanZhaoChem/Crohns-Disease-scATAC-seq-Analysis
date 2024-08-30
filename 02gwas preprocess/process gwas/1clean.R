library(data.table)
library(dplyr)
library(stringr)
annotation_folder <- '~/yuzhao1/work/final_GCAatac/0gwas/b151_hg19/'
catalog_folder <- '~/yuzhao1/work/final_GCAatac/0gwas/catalog/'

####################### A, cd #####################

# 1. combine gwas with annotated rsID and INFO to combined dataframe
gwas<- as.data.frame(
  fread(
    paste0(catalog_folder, "cd_build37_40266_20161107.txt")
    )
  )

annotated_marker_table <- read.table(
  paste0(annotation_folder, "cd_build37_40266_20161107_annotated_hg19_b151.vcf"),
  sep = '\t', header=FALSE, comment.char="#")
annotated_marker_table <- annotated_marker_table[, c(1,2,3,4,5,8)]
colnames(annotated_marker_table) <- c("chr", "pos", "rsID", "A1", "A2", 'INFO')
annotated_marker_table$MarkerName <- paste0(annotated_marker_table$chr, ":", 
                                            annotated_marker_table$pos, "_", 
                                            annotated_marker_table$A1, "_",
                                            annotated_marker_table$A2)
annotation <- dplyr::select(annotated_marker_table, rsID, MarkerName, INFO)

combined <- right_join(annotation, gwas, by = c("MarkerName"))

# process INFO
gwas_info_elements <- strsplit(combined$INFO, split = ';', fixed = T) %>% unlist()
vc_vector <- grep("^VC=", gwas_info_elements, value = TRUE)
vc_vector <- gsub('VC=', '', vc_vector)

caf_vector <- grep("^CAF=", gwas_info_elements, value = TRUE)
caf_vector_numbers <- str_extract_all(caf_vector, "\\b[0-9.]+\\b")
min_values <- sapply(caf_vector_numbers, function(x) {
  if (length(x) == 0) {
    return('NA')
  } else {
    maf_temp <- min(as.numeric(x))
    maf <- min(maf_temp, 1-maf_temp)
    return(maf)
  }
})

# append the MAF and VC
combined$MAF <- min_values
combined$VC <- vc_vector

fwrite(combined, 
       paste0(annotation_folder, 'cd_build37_40266_20161107_annotated_hg19_b151_gwas_annotated.txt'),
       quote=FALSE, sep='\t', row.names=FALSE)



##################### B, uc #####################

# 1. combine gwas with annotated rsID and INFO to combined dataframe
gwas<- as.data.frame(
  fread(
    paste0(catalog_folder, "uc_build37_45975_20161107.txt")
  )
)

annotated_marker_table <- read.table(
  paste0(annotation_folder, "uc_build37_45975_20161107_annotated_hg19_b151.vcf"),
  sep = '\t', header=FALSE, comment.char="#")
annotated_marker_table <- annotated_marker_table[, c(1,2,3,4,5,8)]
colnames(annotated_marker_table) <- c("chr", "pos", "rsID", "A1", "A2", 'INFO')
annotated_marker_table$MarkerName <- paste0(annotated_marker_table$chr, ":", 
                                            annotated_marker_table$pos, "_", 
                                            annotated_marker_table$A1, "_",
                                            annotated_marker_table$A2)
annotation <- dplyr::select(annotated_marker_table, rsID, MarkerName, INFO)

combined <- right_join(annotation, gwas, by = c("MarkerName"))

# process INFO
gwas_info_elements <- strsplit(combined$INFO, split = ';', fixed = T) %>% unlist()
vc_vector <- grep("^VC=", gwas_info_elements, value = TRUE)
vc_vector <- gsub('VC=', '', vc_vector)

caf_vector <- grep("^CAF=", gwas_info_elements, value = TRUE)
caf_vector_numbers <- str_extract_all(caf_vector, "\\b[0-9.]+\\b")
min_values <- sapply(caf_vector_numbers, function(x) {
  if (length(x) == 0) {
    return('NA')
  } else {
    maf_temp <- min(as.numeric(x))
    maf <- min(maf_temp, 1-maf_temp)
    return(maf)
  }
})

# append the MAF and VC
combined$MAF <- min_values
combined$VC <- vc_vector

fwrite(combined, 
       paste0(annotation_folder, 'uc_build37_45975_20161107_annotated_hg19_b151_gwas_annotated.txt'),
       quote=FALSE, sep='\t', row.names=FALSE)







####################### C, ibd #####################

# 1. combine gwas with annotated rsID and INFO to combined dataframe
gwas<- as.data.frame(
  fread(
    paste0(catalog_folder, "ibd_build37_59957_20161107.txt")
  )
)

annotated_marker_table <- read.table(
  paste0(annotation_folder, "ibd_build37_59957_20161107_annotated_hg19_b151.vcf"),
  sep = '\t', header=FALSE, comment.char="#")
annotated_marker_table <- annotated_marker_table[, c(1,2,3,4,5,8)]
colnames(annotated_marker_table) <- c("chr", "pos", "rsID", "A1", "A2", 'INFO')
annotated_marker_table$MarkerName <- paste0(annotated_marker_table$chr, ":", 
                                            annotated_marker_table$pos, "_", 
                                            annotated_marker_table$A1, "_",
                                            annotated_marker_table$A2)
annotation <- dplyr::select(annotated_marker_table, rsID, MarkerName, INFO)

combined <- right_join(annotation, gwas, by = c("MarkerName"))

# process INFO
gwas_info_elements <- strsplit(combined$INFO, split = ';', fixed = T) %>% unlist()
vc_vector <- grep("^VC=", gwas_info_elements, value = TRUE)
vc_vector <- gsub('VC=', '', vc_vector)

caf_vector <- grep("^CAF=", gwas_info_elements, value = TRUE)
caf_vector_numbers <- str_extract_all(caf_vector, "\\b[0-9.]+\\b")
min_values <- sapply(caf_vector_numbers, function(x) {
  if (length(x) == 0) {
    return('NA')
  } else {
    maf_temp <- min(as.numeric(x))
    maf <- min(maf_temp, 1-maf_temp)
    return(maf)
  }
})

# append the MAF and VC
combined$MAF <- min_values
combined$VC <- vc_vector

fwrite(combined, 
       paste0(annotation_folder, 'ibd_build37_59957_20161107_annotated_hg19_b151_gwas_annotated.txt'),
       quote=FALSE, sep='\t', row.names=FALSE)





