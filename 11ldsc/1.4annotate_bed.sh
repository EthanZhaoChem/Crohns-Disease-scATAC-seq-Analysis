module load python
source activate  /home/yuzhao1/yuzhao1/software/conda/ldsc
bimfile_path=/project/spott/yuzhao1/ldsc/minimum_required/1000G_EUR_Phase3_plink
script_path=~/yuzhao1/work/atac_gca2024/14ldsc
cd ~/yuzhao1/work/atac_gca2024/14ldsc/jobs

# need to modify
DirRoot=/home/yuzhao1/yuzhao1/work/atac_gca2024/14ldsc/results/union_sub100_k45_daPeaks_positive_flexibleLpval30k_vsnull
bedfilesRootDir=$DirRoot/bed_lifted_extended
annot_path=$DirRoot/annots

if [ ! -d $annot_path ]
then
    mkdir $annot_path
fi


names=`ls  $bedfilesRootDir | cut -f 1 -d '.'`
for name in $names
do
    bedname=`echo $name | awk '{print $1}'`
    echo $bedname
    if [ ! -d $annot_path/$bedname ]
    then
        mkdir $annot_path/$bedname
    fi
    cmd="python  $script_path/1customized_make_annot.py --bedname $bedname --bedfile_path $bedfilesRootDir --bimfile_path $bimfile_path --annot_path $annot_path/$bedname"
    sbatch -J $bedname --time=12:00:00 --mem=30G --output=annot_$bedname.out --error=annot_$bedname.err --account=pi-spott  -p caslake -c 1 --wrap="$cmd"
done




