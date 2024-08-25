#!/bin/bash
permID=perm3

vcf_dir=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/1vcf/output_as
celltype_file=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt
rasqual_input_dir=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/4calculation_all_peaks/4ctSpecific_rasqualInput
script_submission=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/4calculation_all_peaks/4.3perm_wrapper.batch
rasqual_output_dir_root=/scratch/midway3/yuzhao1/rasqual/allPeaks_noDiseaseCo/$permID
module load python
module load gsl
module load gcc/7.4.0
cd /project/gca/yuzhao1/work/atac_gca2024/19rasqual

# can customize it
n_blocks=10
tmp_dir=log/$permID


if [ ! -d $tmp_dir ]
then
    mkdir $tmp_dir
fi

for ct in `cat $celltype_file | awk '{print $1}' | sort -V | uniq`;
do
  output_cell=$rasqual_output_dir_root/$ct
  if [ ! -d $output_cell ]
  then
      mkdir $output_cell
  fi
  
  input_txt=$rasqual_input_dir/$ct.input.txt
  n_lines=$(wc -l $input_txt | cut -d" " -f1)
  blocksize=$(( $n_lines/($n_blocks - 1) )) # default is floor rounding
  blocksize_end=$(( $n_lines-($n_blocks-1)*$blocksize ))

  for block_index in $(seq $n_blocks);
  do 
      if [ $block_index -lt $n_blocks ]
      then
        lines=$(seq $((1+($block_index-1)*$blocksize)) $(($block_index*$blocksize)) )
        if [ -z "$lines" ]; then
            echo "lines is empty"
            continue
        fi
        
        tmp=$tmp_dir/${ct}_${block_index}.tmpLines.txt
        echo $lines > $tmp
        
        my_message=$(echo block index $block_index , range $((1+($block_index-1)*$blocksize)) -  $(($block_index*$blocksize)) )
        cmd="echo $my_message
        bash $script_submission $ct $block_index"

      else
        # last block
        lines=$(seq $((1+($block_index-1)*$blocksize)) $n_lines)
        if [ -z "$lines" ]; then
            echo "lines is empty"
            continue
        fi
        
        tmp=$tmp_dir/${ct}_${block_index}.tmpLines.txt
        echo $lines > $tmp
        
        my_message=$(echo block index $block_index , range, $((1+($block_index-1)*$blocksize)) -  $n_lines)
        cmd="echo $my_message
        bash $script_submission $ct $block_index"
      fi
      
      sbatch -J ${permID}_${ct}_${block_index} --time=36:00:00 --mem=5G --output=$tmp_dir/${ct}_${block_index}.out --error=$tmp_dir/${ct}_${block_index}.err --account=pi-onibasu  -p caslake -c 10 --wrap="$cmd"
  done
done 



