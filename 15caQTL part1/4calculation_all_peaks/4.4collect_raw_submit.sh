#!/bin/bash
script_submission=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/4calculation_all_peaks/4.4collect_raw_run.sh
celltype_file=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt

for ct in `cat $celltype_file | awk '{print $1}' | sort -V | uniq`;
do

  for chr in {1..22} X
  do
    cmd="bash $script_submission $ct $chr"
    sbatch -J "perm2${ct}_chr${chr}"  -o "log/${ct}_chr${chr}.out" -e "log/${ct}_chr${chr}.err" \
    --time=12:00:00 --mem=5G --account=pi-onibasu -p caslake -c 1 --wrap="$cmd"
  done

done

