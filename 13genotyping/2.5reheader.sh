sampleIDs_taskfile=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.1sampleIDs_all
input_cell=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.4concat
output_cell=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.5reheader

tmp="/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.5reheader/tmp.txt"
for sampleID in `cat $sampleIDs_taskfile | awk '{print $1}' | sort -V | uniq`;
do
    echo $sampleID >> $tmp
    bcftools reheader $input_cell/${sampleID}.vcf.gz -s $tmp -o $output_cell/${sampleID}.vcf.gz
    tabix -p vcf $output_cell/${sampleID}.vcf.gz
    rm $tmp
done

