#!/bin/bash

module load R/4.1.0
module load proj
module load gsl

# manual
input_file=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/4.0prepare_samplePair_all_chr1.txt
dir_log=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/4.3sample_concordance_all_chr1/log
script_submission=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/4.3sample_concordance_all_chr1.R
dir_out=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/4.3sample_concordance_all_chr1



# Read each line from the file
while read -r line; do
  # Split the line into two parts and assign them to variables 'a' and 'b'
  IFS=' ' read -ra ADDR <<< "$line"
  a="${ADDR[0]}"
  b="${ADDR[1]}"
  
  if [ ! -f $dir_out/${a}_${b}.rds ]
  then
  cmd="Rscript --vanilla $script_submission $a $b"
  sbatch -J ${a}_${b} --time=36:00:00 --mem=20G --output=$dir_log/${a}_${b}.out --error=$dir_log/${a}_${b}.err --account=pi-spott  -p caslake -c 1 --wrap="$cmd"
  fi
done < "$input_file"












