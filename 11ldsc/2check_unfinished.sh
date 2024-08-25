module load python
source activate  /home/yuzhao1/yuzhao1/software/conda/ldsc

# fixed
bfile_path=/project/spott/yuzhao1/ldsc/minimum_required/1000G_EUR_Phase3_plink
ldsc_path=/home/yuzhao1/yuzhao1/software/ldsc
snps_path=/project/spott/yuzhao1/ldsc/minimum_required

# need to modify
annot_rootDir=/home/yuzhao1/yuzhao1/work/atac_gca2024/14ldsc/results/union_sub100_k45_daPeaks_positive_flexibleLpval30k_vsnull/annots

# script
cd ~/yuzhao1/work/atac_gca2024/14ldsc/jobs
flag=0

for singleAnnotFolder in `ls $annot_rootDir | awk '{print $1}' | sort | uniq`;
do
    annot_dir=`echo $singleAnnotFolder | awk '{print $1}'`

    for chrom in {1..22}
    do
        if [ ! -f $annot_rootDir/$annot_dir/$annot_dir.$chrom.l2.ldscore.gz ]
        then
            echo unfinished_${annot_dir}_chr${chrom}
            flag=1
        fi
    done 
done


if [[ $flag == 1 ]]; then
    echo "The tasks above are not finished."
else
    echo "All jobs are done. Congratulations!"
fi


