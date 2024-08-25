# see the original script in rasqual package, this is a modified version for our analysis only (paired end atac)

#!/bin/bash
celltype_file=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt

PAIRED_OR_SINGLE_END=paired_end
RASQUALDIR=/project/gca/yuzhao1/software/rasqual/rasqual-master
VCF_INPUT_dir=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/1vcf/maf_filtered
VCF_OUTPUT_dir=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/1vcf/output_as
ASSAY_TYPE="atac"
BAM_LIST_file=/project/gca/yuzhao1/work/atac_gca2024/20bam/1bam_perBC/celltype_patient_bams


export RASQUALDIR=/project/gca/yuzhao1/software/rasqual/rasqual-master

for ct in `cat $celltype_file | awk '{print $1}' | sort -V | uniq`;
do
	cmd="bash $RASQUALDIR/src/ASVCF/createASVCF.sh $PAIRED_OR_SINGLE_END $BAM_LIST_file/$ct/bam.list $VCF_INPUT_dir/${ct}.vcf.gz $VCF_OUTPUT_dir/${ct}.vcf.gz $ASSAY_TYPE"
  sbatch -J $ct --time=36:00:00 --mem=50G --output=log/$ct.out --error=log/$ct.err --account=pi-spott  -p caslake -c 1 --wrap="$cmd"
  # exit 1
done 



