library(stringr)

patients_all <- readLines('/project/gca/yuzhao1/work/final_GCAatac/18glimpse/3.5merge/samples.txt')
dir_merge <- '/project/gca/yuzhao1/work/final_GCAatac/18glimpse/3.1patient_bam_merge/'


sink("~/yuzhao1/work/final_GCAatac/18glimpse/3.5print_patient_bam_list.txt", append = F)
for (patient in patients_all) {
  merged_bam <- paste0(dir_merge, patient, '.bam')
  cat(merged_bam, "\n")
}
sink()




