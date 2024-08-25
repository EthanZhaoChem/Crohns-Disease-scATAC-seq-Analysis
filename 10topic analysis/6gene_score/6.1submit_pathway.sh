#!/bin/bash
module load R/4.1.0
script_submission=/project/gca/yuzhao1/work/atac_gca2024/13fasttopic/6gene_score/6.1union_calculate.R


for jobid in {1..45}
# for jobid in 33
do
  cmd="Rscript --vanilla $script_submission $jobid "
  sbatch -J "job${jobid}"  -o "log/job${jobid}.out" \
  -e "log/job${jobid}.err" \
  --time=12:00:00 --mem=80G --account=pi-spott -p caslake -c 1 --wrap="$cmd"
done

