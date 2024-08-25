#!/bin/bash
dir_input=/project/gca/yuzhao1/software/glimpse2/ref_panel
dir_output=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/3pca/1kg_vcf_single


for chrom in {1..22} X;
do
  cmd="bcftools annotate -Oz -o $dir_output/1000GP.chr${chrom}.vcf.gz -x ID  --set-id +'%CHROM\_%POS\_%REF\_%ALT' $dir_input/1000GP.chr${chrom}.bcf
tabix -p vcf $dir_output/1000GP.chr${chrom}.vcf.gz"
  sbatch -J "$chrom" --time=12:00:00 --mem=5G --output="log/$chrom.out" --error="log/$chrom.err" --account=pi-spott  -p spott -c 1 --wrap="$cmd"
done 

