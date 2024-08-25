library(stringr)
metadata <- read.table('~/yuzhao1/work/atac_gca2024/0metadata/meta_Ethan_curated_20240211.csv',
                       header = T, sep = ',')

dir_merge <- '/project/gca/yuzhao1/work/final_GCAatac/18glimpse/3.1patient_bam_merge/'
dir_cellranger <- '/project/spott/yuzhao1/GEO_GCAatac/cellranger_output/'
patients_all <- unique(metadata$patient_masked)

sink("~/yuzhao1/work/final_GCAatac/18glimpse/3.1patient_bam_merge.sh", append = F)
cat('cd ~/yuzhao1/work/final_GCAatac/18glimpse\n')
for (patient in patients_all) {
  files_bam <- paste0(dir_cellranger, metadata$sample[metadata$patient_masked == patient], '/outs/possorted_bam.bam')
  merged_bam <- paste0(dir_merge, patient, '.bam')
  cmd1 <- paste(sep = ' ', 'samtools merge -o', merged_bam, paste(files_bam, collapse = ' '))
  cmd2 <- paste(sep = ' ', 'samtools index ', merged_bam)
  cmd <- paste(sep = '\n', cmd1, cmd2)
  cat(sep = '', 'sbatch -J ', patient, ' --time=36:00:00 --mem=5Gb --output=3.1log/', patient, '.out --error=3.1log/', patient, ".err --account=pi-spott  -p spott -c 1 --wrap='", cmd, "'\n")
}
sink()




