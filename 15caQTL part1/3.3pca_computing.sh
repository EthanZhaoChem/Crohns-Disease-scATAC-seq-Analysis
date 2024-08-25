#!/bin/bash
module load plink
cd /project/gca/yuzhao1/work/atac_gca2024/19rasqual/3pca/plink_analysis

# i used '_' in ref ID, so here I will change the snpID from out samples to ':' separated to match the ref vcf
vcf_sample_rawID=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/7.1snps_filter/allPatients_filtered.vcf.gz 
bcftools annotate -Oz -o sample.chrALL.vcf.gz -x ID  --set-id +'%CHROM\_%POS\_%REF\_%ALT' $vcf_sample_rawID

# just want to copy to the same folder (to make the folder structure look nice..)
cp /project/gca/yuzhao1/work/atac_gca2024/19rasqual/3pca/1kg_vcf_concat/1000GP.chrALL.vcf.gz ref.chrALL.vcf.gz

vcf_sample=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/3pca/plink_analysis/sample.chrALL.vcf.gz
vcf_ref=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/3pca/plink_analysis/ref.chrALL.vcf.gz

plink --vcf $vcf_sample --make-bed --out sample
plink --vcf $vcf_ref --make-bed --out ref

awk '{print $2}'  sample.bim | sort >  sample_snps.txt
awk '{print $2}' ref.bim | sort > ref_snps.txt

comm -12  sample_snps.txt ref_snps.txt > intersecting_snps_sample_1kg.txt

plink --bfile ref --extract intersecting_snps_sample_1kg.txt --indep 50 5 2 --make-bed --out ref.subset

plink --bfile ref.subset --extract intersecting_snps_sample_1kg.txt --bmerge  sample.bed  sample.bim  sample.fam --make-bed --out all_sample_1KGP

plink --bfile all_sample_1KGP \
  --extract ref.subset.prune.in \
  --out final_sample_1KGP_chrALL \
  --pca tabs header



