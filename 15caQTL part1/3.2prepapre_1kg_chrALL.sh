#!/bin/bash
#SBATCH --job-name=co
#SBATCH --output=log/%x.out
#SBATCH --error=log/%x.err

#SBATCH --time=36:00:00
#SBATCH --partition=spott
#SBATCH --account=pi-spott
#SBATCH --tasks=8
#SBATCH --mem=50Gb

dir_output=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/3pca/1kg_vcf_single
dir_output2=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/3pca/1kg_vcf_concat

xx=$(ls -1v $dir_output/*.vcf.gz)
bcftools concat --threads 8 -Oz -o $dir_output2/1000GP.chrALL.vcf.gz  $xx 
