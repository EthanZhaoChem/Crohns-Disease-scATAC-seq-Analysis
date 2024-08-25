# see the original script in rasqual package, this is a modified version for our analysis only (paired end atac)

#!/bin/bash
celltype_file=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt
VCF_INPUT_dir=/project/gca/yuzhao1/work/atac_gca2024/20bam/2.2vcf_celltype/vcf_perCelltype_filtered_patients
VCF_OUTPUT_dir=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/1vcf/maf_filtered


for ct in `cat $celltype_file | awk '{print $1}' | sort -V | uniq`;
do
  vcf1=$VCF_INPUT_dir/${ct}.vcf.gz
  vcf2=$VCF_OUTPUT_dir/${ct}.vcf.gz
	cmd="bcftools view -q 0.05:minor $vcf1 -Oz -o $vcf2
	tabix -p vcf $vcf2"
  sbatch -J $ct --time=120:00:00 --mem=10G --output=log/$ct.out --error=log/$ct.err --account=pi-spott  -p spott -c 1 --wrap="$cmd"
  # exit 1
done


