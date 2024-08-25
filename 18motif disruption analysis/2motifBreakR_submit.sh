#!/bin/bash

script_submission=/project/gca/yuzhao1/work/atac_gca2024/26motif_disruption3/2motifBreakR_run.R
output_cell=/project/gca/yuzhao1/work/atac_gca2024/26motif_disruption3/2results_allHumanTF
module load R/4.1.0
module load proj
module load gsl

for bedChunkID in {1..200};
do
    if [ ! -f $output_cell/$bedChunkID.rds ]
    then
    cmd="Rscript --vanilla $script_submission $bedChunkID"
    sbatch -J "$bedChunkID" --time=36:00:00 --mem=50G --output="log/$bedChunkID.out" --error="log/$bedChunkID.err" --account=pi-spott  -p caslake -c 1 --wrap="$cmd"
    fi
done 

