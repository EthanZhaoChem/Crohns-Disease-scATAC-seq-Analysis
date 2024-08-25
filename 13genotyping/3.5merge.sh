dir_output=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/3.5merge
sampleIDs_taskfile=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/3.1patientIDs_all
prefix="/project/gca/yuzhao1/work/final_GCAatac/18glimpse/3.4concat/"
suffix=".vcf.gz"
vcf_list="/project/gca/yuzhao1/work/final_GCAatac/18glimpse/3.5merge/tmp.txt"
samples_list="/project/gca/yuzhao1/work/final_GCAatac/18glimpse/3.5merge/samples.txt"


rm $vcf_list
rm $samples_list

for sampleID in `cat $sampleIDs_taskfile | awk '{print $1}' | sort -V | uniq`;
do
  tmp="${prefix}${sampleID}${suffix}"
  echo "$tmp" >> "$vcf_list"
  echo "$sampleID" >> "$samples_list"
done 


cmd="bcftools merge -m none -Oz -o $dir_output/allPatients.vcf.gz -l $vcf_list
bcftools index -f $dir_output/allPatients.vcf.gz"
sbatch --time=36:00:00 --mem=50Gb --output=merge.log --error=merge.err --account=pi-spott  -p caslake -c 1 --wrap="$cmd"
