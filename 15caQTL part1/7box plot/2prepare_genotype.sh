#!/bin/bash
celltype_file=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt
VCF_INPUT_dir=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/1vcf/output_as
VCF_OUTPUT_dir=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/9boxplot/9.2genotype_raw


for ct in `cat $celltype_file | awk '{print $1}' | sort -V | uniq`;
do
  bcftools query -l $VCF_INPUT_dir/$ct.vcf.gz > $VCF_OUTPUT_dir/$ct.samples
# 	cmd="bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT[\\t%GT]\\n' $VCF_INPUT_dir/$ct.vcf.gz > $VCF_OUTPUT_dir/$ct.genotype"
#   sbatch -J $ct --time=36:00:00 --mem=5G --output=log/$ct.out --error=log/$ct.err --account=pi-spott  -p spott -c 1 --wrap="$cmd"
done 

