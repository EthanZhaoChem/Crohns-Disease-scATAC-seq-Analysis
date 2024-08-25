#!/bin/bash
# nCovariates is automatically calculated by size of X file and nSamples in C.
rasqual_input_dir=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/4calculation_all_peaks/4ctSpecific_rasqualInput
nPatients_perCt_dir=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/0patientID_perCelltype
param_file=$1
line_num=$2
Y=$3
K=$4
X=$5
vcf_file=$6
out_dir=$7
ct=$8

param=($(cat $param_file | sed "${line_num}q;d"))
gene_id=${param[0]}
gene_name=${param[1]}
region=${param[2]}
n_rsnp=${param[3]}
n_fsnp=${param[4]}
exon_start_positions=${param[5]}
exon_end_positions=${param[6]}

feat_id=$(grep $gene_id -n $rasqual_input_dir/$ct.accessibility.txt | cut -d":" -f1,1)
n_sample=$(cat $nPatients_perCt_dir/$ct.txt | wc -l)
window_size=20502

if [[ -s $out_dir/${gene_id}_${gene_name}.txt ]]; then
	# echo $out_dir/${gene_id}_${gene_name}.txt exists! skipping...
	exit 1
else
  echo id: $gene_id
  echo name: $gene_name
  echo region: $region
  echo reference snps: $n_rsnp
  echo feature snps: $n_fsnp
  echo feature id: $feat_id
  
	tabix $vcf_file $region | \
	rasqual \
	-y $Y -k $K -x $X \
	-n $n_sample -j $feat_id -l $n_rsnp -m $n_fsnp \
	-s $exon_start_positions -e $exon_end_positions \
	--cis-window-size $window_size --feature-name $gene_name \
	--n_threads 10 -v --force --random-permutation > $out_dir/${gene_id}_${gene_name}.txt
fi

echo --- --- --- --- --- --- 






