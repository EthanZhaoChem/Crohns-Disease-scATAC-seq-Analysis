#!/bin/bash

# fixed
module load R/4.1.0
export OMP_NUM_THREADS=1
sumstats_path=/project/spott/yuzhao1/ldsc/minimum_required/sumstats
script_submission=~/yuzhao1/work/atac_gca2024/14ldsc/4.1post_process.R

# need to modify
root_dir=/home/yuzhao1/yuzhao1/work/atac_gca2024/14ldsc/results/union_sub100_k45_daPeaks_positive_flexibleLpval30k_vsnull

# fixed
sumstats_taskfile=~/yuzhao1/work/atac_gca2024/14ldsc/3sumstats_tasks
annot_cell=$root_dir/annots
results_cell=$root_dir/output_conditioned
post_process_cell=$root_dir/post_process_conditioned

if [ ! -d $post_process_cell ]
then
    mkdir $post_process_cell
fi

cd ~/yuzhao1/work/atac_gca2024/14ldsc/jobs
for step in `cat $sumstats_taskfile | awk '{print $1}' | sort | uniq`;
do
    trait2=`echo $step | awk '{print $1}'`
    echo $trait2
    
    if [ ! -f ${post_process_cell}/${trait2}_ldsc_postprocess.txt ]
    then
    cmd="Rscript --vanilla $script_submission $trait2 $annot_cell $results_cell $post_process_cell"
    sbatch -J $trait2 --time=12:00:00 --mem=30G --output=post_$trait2.out --error=post_$trait2.err --account=pi-spott  -p spott -c 1 --wrap="$cmd"
    fi
    
done




