dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(ArchR)

# write the cell type specific patients to txt files
dir_patientsPerCt <- '~/yuzhao1/work/atac_gca2024/19rasqual/0patientID_perCelltype/'
celltypes <- readLines('~/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt')
out.dir <- '~/yuzhao1/work/atac_gca2024/20bam/2.2vcf_celltype/vcf_perCelltype_filtered_patients/'

addArchRThreads(1)

# write commands to subset vcf cell type wise to a bash file
original_vcf <- '~/yuzhao1/work/final_GCAatac/18glimpse/7.1snps_filter/allPatients_filtered.vcf.gz'
sink("~/yuzhao1/work/atac_gca2024/20bam/2.2vcf_celltype.sh", append = F)
cat('cd ~/yuzhao1/work/atac_gca2024/20bam\n')

for (ct in celltypes) {
  patients_file <- paste0(dir_patientsPerCt, ct, '.txt')
  vcf_sub_file <- paste0(out.dir, ct, '.vcf.gz')
  cmd1 <- paste(sep = ' ', 'bcftools view -Oz -S', patients_file, original_vcf, '>', vcf_sub_file)
  cmd2 <- paste(sep = ' ', 'bcftools index', vcf_sub_file)
  cmd <- paste(sep = '\n', cmd1, cmd2)
  cat(sep = '', 'sbatch -J ', ct, ' --time=36:00:00 --mem=5Gb --output=2.2vcf_celltype/log/', ct, '.out --error=2.2vcf_celltype/log/', ct, ".err --account=pi-spott  -p spott -c 1 --wrap='", cmd, "'\n")

  }
sink()












