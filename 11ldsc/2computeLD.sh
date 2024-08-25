module load python
source activate  /home/yuzhao1/yuzhao1/software/conda/ldsc

bfile_path=/project/spott/yuzhao1/ldsc/minimum_required/1000G_EUR_Phase3_plink
ldsc_path=/home/yuzhao1/yuzhao1/software/ldsc
snps_path=/project/spott/yuzhao1/ldsc/minimum_required
cd ~/yuzhao1/work/atac_gca2024/14ldsc/jobs

# need to modify
annot_rootDir=/home/yuzhao1/yuzhao1/work/atac_gca2024/14ldsc/results/union_sub100_k45_daPeaks_positive_flexibleLpval30k_vsnull/annots

for singleAnnotFolder in `ls $annot_rootDir | awk '{print $1}' | sort | uniq`;
do
    annot_dir=`echo $singleAnnotFolder | awk '{print $1}'`
    echo $annot_dir
    
    for chrom in {1..22}
    do
        if [ ! -f $annot_rootDir/$annot_dir/$annot_dir.$chrom.l2.ldscore.gz ]
        then
            cmd="python $ldsc_path/ldsc.py --bfile $bfile_path/1000G.EUR.QC.$chrom --l2 --ld-wind-cm 1 --yes-really --annot $annot_rootDir/$annot_dir/$annot_dir.$chrom.annot.gz --print-snps $snps_path/hm3_no_MHC.list.txt --out $annot_rootDir/$annot_dir/$annot_dir.$chrom"
            sbatch -J ld_${annot_dir}_chr${chrom} --time=300:00 --mem=30G --output=ld_${annot_dir}_chr${chrom}.out --account=pi-spott --error=ld_${annot_dir}_chr${chrom}.err -p caslake -c 1 --wrap="$cmd "
        fi
    done 
    
done



