dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(gtools)
dir_singleBam <- '~/yuzhao1/work/atac_gca2024/20bam/1bam_perBC/singleBams/'
dir_mergeBam <- '/project/gca/yuzhao1/work/atac_gca2024/20bam/1bam_perBC/celltype_patient_bams/'
dir_patientsPerCt <- '~/yuzhao1/work/atac_gca2024/19rasqual/0patientID_perCelltype/'
celltypes <- readLines('~/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt')
metadata <- read.table('~/yuzhao1/work/atac_gca2024/0metadata/meta_Ethan_curated_20240211.csv',
                       header = T, sep = ',')
patients_all <- unique(metadata$patient_masked)



#!!!!! we only prepared the patient files that are filtered for each cell type!!!!!

sink("~/yuzhao1/work/atac_gca2024/20bam/1.3merge_bam.sh", append = F)
cat('cd ~/yuzhao1/work/atac_gca2024/20bam\n')

for (ct in celltypes) {
  dir.create(paste0(dir_mergeBam, ct), showWarnings = F)
  patients <- readLines(paste0(dir_patientsPerCt, ct, '.txt'))
  
  ct_bam_list_filename <- paste0(dir_mergeBam, ct, '/bam.list')
  ct_bam_list <- paste0(dir_mergeBam, ct, '/', patients, '.bam')
  writeLines(ct_bam_list, ct_bam_list_filename)
  
  for (patient in mixedsort(patients)) {
    samples <- metadata$sample[metadata$patient_masked == patient]
    files_bam <- paste0(dir_singleBam, samples, '/', ct, '.bam')
    files_bam <- files_bam[file.exists(files_bam)]
    merged_bam <- paste0(dir_mergeBam, ct, '/', patient, '.bam')
    cmd1 <- paste(sep = ' ', 'samtools merge -o', merged_bam, paste(files_bam, collapse = ' '))
    cmd2 <- paste(sep = ' ', 'samtools index ', merged_bam)
    cmd <- paste(sep = '\n', cmd1, cmd2)
    jobname <- paste0(ct, '-', patient)
    cat(sep = '', 'sbatch -J ',  jobname, ' --time=36:00:00 --mem=5Gb --output=log/', jobname, '.out --error=log/', jobname, ".err --account=pi-spott  -p caslake -c 1 --wrap='", cmd, "'\n")
  }
}

sink()






