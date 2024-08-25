#!/bin/bash
ct=$1
chr=$2
rasqual_output_dir_root=/scratch/midway3/yuzhao1/rasqual/allPeaks_noDiseaseCo/perm3
Dir_save=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/4calculation_all_peaks/4ctSpecific_compliedOutput/perm3/$ct

if [ ! -d $Dir_save ]
then
    mkdir -p $Dir_save
fi

[ -e file ] && rm $Dir_save/chr${chr}

for file in $rasqual_output_dir_root/$ct/chr$chr*; do
    cat "$file" >> $Dir_save/chr${chr}
done



