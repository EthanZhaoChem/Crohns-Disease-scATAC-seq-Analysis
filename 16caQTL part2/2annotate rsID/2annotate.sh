# download path: https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151***

module load java

cat snps_raw.txt  | cut -f1 | \
awk '{n=split($0,a,/[:_]/); print a[1]"\t"a[2]"\t.\t"a[3]"\t"a[4]"\t.\t.\t."}'|  \
sort -k1,1V -k2,2g | \
java -jar /home/yuzhao1/yuzhao1/software/snpeff/snpEff/SnpSift.jar \
  annotate -a  /home/yuzhao1/yuzhao1/resource/vcf/human_9606_b151_GRCh38p7/00-All.vcf.gz  \
  > snps_annotated.txt