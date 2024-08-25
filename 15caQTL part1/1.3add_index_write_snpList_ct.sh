#!/bin/bash
#SBATCH --job-name=snp
#SBATCH --output=log/%x.out
#SBATCH --error=log/%x.err
#SBATCH --time=36:00:00
#SBATCH --partition=spott
#SBATCH --account=pi-spott
#SBATCH --tasks=1
#SBATCH --mem=50Gb

celltype_file=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt
VCF_dir=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/1vcf/output_as
snp_OUTPUT_dir=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/1vcf/filtered_snpList


for ct in `cat $celltype_file | awk '{print $1}' | sort -V | uniq`;
do
  vcf=$VCF_dir/${ct}.vcf.gz
  tabix -p vcf $vcf
  snpList=$snp_OUTPUT_dir/${ct}.snp.list
  bcftools query -f '%ID\n' $vcf > $snpList
done 



