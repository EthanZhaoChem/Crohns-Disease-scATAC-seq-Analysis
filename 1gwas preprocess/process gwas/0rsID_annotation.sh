# download path: https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151***

module load java
java -jar /home/yuzhao1/yuzhao1/software/snpeff/snpEff/SnpSift.jar annotate -h

# cd
cut -f1 ~/yuzhao1/work/final_GCAatac/0gwas/catalog/cd_build37_40266_20161107.txt | tail -n+2 | \
  awk '{n=split($0,a,/[:_]/); print a[1]"\t"a[2]"\t.\t"a[3]"\t"a[4]"\t.\t.\t."}'  | \
  sort -k1,1V -k2,2g | \
  java -jar /home/yuzhao1/yuzhao1/software/snpeff/snpEff/SnpSift.jar \
  annotate -a /home/yuzhao1/yuzhao1/resource/vcf/human_9606_b151_GRCh37p13/00-All.vcf.gz  \
  > ~/yuzhao1/work/final_GCAatac/0gwas/b151_hg19/cd_build37_40266_20161107_annotated_hg19_b151.vcf


# uc
cut -f1 ~/yuzhao1/work/final_GCAatac/0gwas/catalog/uc_build37_45975_20161107.txt | tail -n+2 | \
  awk '{n=split($0,a,/[:_]/); print a[1]"\t"a[2]"\t.\t"a[3]"\t"a[4]"\t.\t.\t."}'  | \
  sort -k1,1V -k2,2g | \
  java -jar /home/yuzhao1/yuzhao1/software/snpeff/snpEff/SnpSift.jar \
  annotate -a /home/yuzhao1/yuzhao1/resource/vcf/human_9606_b151_GRCh37p13/00-All.vcf.gz  \
  > ~/yuzhao1/work/final_GCAatac/0gwas/b151_hg19/uc_build37_45975_20161107_annotated_hg19_b151.vcf
  
# ibd
cut -f1 ~/yuzhao1/work/final_GCAatac/0gwas/catalog/ibd_build37_59957_20161107.txt | tail -n+2 | \
  awk '{n=split($0,a,/[:_]/); print a[1]"\t"a[2]"\t.\t"a[3]"\t"a[4]"\t.\t.\t."}'  | \
  sort -k1,1V -k2,2g | \
  java -jar /home/yuzhao1/yuzhao1/software/snpeff/snpEff/SnpSift.jar \
  annotate -a /home/yuzhao1/yuzhao1/resource/vcf/human_9606_b151_GRCh37p13/00-All.vcf.gz  \
  > ~/yuzhao1/work/final_GCAatac/0gwas/b151_hg19/ibd_build37_59957_20161107_annotated_hg19_b151.vcf
  
  
  
  
  