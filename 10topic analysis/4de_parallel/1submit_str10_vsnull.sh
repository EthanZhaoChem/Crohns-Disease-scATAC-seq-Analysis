#!/bin/bash

# usually no need to change
Dir_root=/project/gca/yuzhao1/work/atac_gca2024/13fasttopic/rds
script_submission=/project/gca/yuzhao1/work/atac_gca2024/13fasttopic/4de_parallel/1run_vsnull.R
nThreads=20

# need to change
Dir_tmp=$Dir_root/fit_stromal_sub100_k10_converged_de_vsnull
file_counts=$Dir_root/stromal_raw_sub100.RData
file_fasttopic=$Dir_root/fit_stromal_sub100_k10_converged.rds
jobTitle=str10

# no need to change
module load R/4.1.0
module load proj
module load gsl


if [ ! -d $Dir_tmp ]
then
  mkdir $Dir_tmp
fi

if [ ! -d $Dir_tmp/log ]
then
  mkdir $Dir_tmp/log
fi


for i in `seq 1 $nThreads`;
do

  if [ -f $Dir_tmp/$i.rds ]; then
    continue
  fi

  cmd="Rscript --vanilla $script_submission $file_counts $file_fasttopic $Dir_tmp $i $nThreads"
  sbatch -J "${jobTitle}_${i}_${nThreads}" --time=36:00:00 --mem=80G --output="$Dir_tmp/log/${jobTitle}_${i}_${nThreads}.out" --error="$Dir_tmp/log/${jobTitle}_${i}_${nThreads}.err" --account=pi-spott  -p caslake -c 4 --wrap="$cmd"
done


