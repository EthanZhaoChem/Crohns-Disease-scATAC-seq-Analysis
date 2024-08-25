module load python
source activate  /home/yuzhao1/yuzhao1/software/conda/ldsc

ldsc_path=/home/yuzhao1/yuzhao1/software/ldsc
baseline_path=/project/spott/yuzhao1/ldsc/minimum_required/1000G_Phase3_baselineLD_v2.2_ldscores
weights_path=/project/spott/yuzhao1/ldsc/minimum_required/1000G_Phase3_weights_hm3_no_MHC
freq_path=/project/spott/yuzhao1/ldsc/minimum_required/1000G_Phase3_frq
sumstats_path=/project/spott/yuzhao1/ldsc/minimum_required/sumstats

# need to modify
DirRoot=/home/yuzhao1/yuzhao1/work/atac_gca2024/14ldsc/results/union_sub100_k45_daPeaks_positive_flexibleLpval30k_vsnull

# fixed
annot_rootDir=$DirRoot/annots
output_cell=$DirRoot/output_conditioned
sumstats_taskfile=~/yuzhao1/work/atac_gca2024/14ldsc/3sumstats_tasks
allPeaks_prefix=~/yuzhao1/work/atac_gca2024/14ldsc/results/union_peakset/annots/union_peakset/union_peakset.

if [ ! -d $output_cell ]
then
    mkdir $output_cell
fi

cd ~/yuzhao1/work/atac_gca2024/14ldsc/jobs
echo $output_cell

for singleAnnotFolder in `ls $annot_rootDir | awk '{print $1}' | sort | uniq`;
do
    annot_dir=`echo $singleAnnotFolder | awk '{print $1}'`
    echo $annot_dir
    
    if [ ! -d $output_cell/$annot_dir ]
    then
            mkdir $output_cell/$annot_dir
    fi
    
    for step in `cat $sumstats_taskfile | awk '{print $1}' | sort | uniq`;
    do
        sumstats_file=`echo $step | awk '{print $1}'`
        echo $sumstats_path $sumstats_file
        
        if [ ! -f $sumstats_path/$sumstats_file ]
        then
        echo "Error: sumstats file not found" > ldsc_logfile.log
        exit 102
        fi
        
        if [ ! -f $output_cell/$annot_dir/$sumstats_file.results ]
        then
        cmd="python $ldsc_path/ldsc.py  --h2 $sumstats_path/$sumstats_file --ref-ld-chr $annot_rootDir/$annot_dir/$annot_dir.,${allPeaks_prefix},$baseline_path/baselineLD.  --frqfile-chr $freq_path/1000G.EUR.QC. --w-ld-chr $weights_path/weights.hm3_noMHC. --overlap-annot --print-coefficients --print-delete-vals  --out $output_cell/$annot_dir/$sumstats_file"
        sbatch -J h2_${annot_dir}_${sumstats_file} --time=12:00:00 --mem=30G --output=h2_${annot_dir}_${sumstats_file}.out --error=h2_${annot_dir}_${sumstats_file}.err --account=pi-spott  -p caslake -c 1 --wrap="$cmd"
        fi
        
    done
done
    
    




