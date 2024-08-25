#!/bin/bash
module load proj
module load gsl
module load hdf5
module load R/4.2.0

scriptPath=/project/gca/yuzhao1/work/atac_gca2024/7dreamlet/2pairwise_dreamlet_union.R
lineage=union

for i in 1 2 3;
do
  cmd="Rscript --vanilla $scriptPath $i"
  sbatch -J "${lineage}_${i}" -o log/${lineage}_${i}.out -e log/${lineage}_${i}.err --time=300:00:00 --mem=60G --account=pi-spott -p spott -c 1 --wrap="$cmd" 
done










