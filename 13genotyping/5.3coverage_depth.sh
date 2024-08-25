#!/bin/bash

cd /project/gca/yuzhao1/work/final_GCAatac/18glimpse

# Specify the directory containing your BAM files
BAM_DIR="/project/gca/yuzhao1/work/final_GCAatac/18glimpse/3.1patient_bam_merge"

# Output file to store the results
OUTPUT_dir="/project/gca/yuzhao1/work/final_GCAatac/18glimpse/5.3coverage_depth"

if [ ! -d $OUTPUT_dir/log ]
then
        mkdir $OUTPUT_dir/log
fi
    
# Iterate over each BAM file in the directory
for BAM_FILE in "$BAM_DIR"/*.bam; do
    # Extract the filename without the directory path
    BAM_FILENAME=$(basename "$BAM_FILE")
    OUTPUT_FILE=$OUTPUT_dir/$BAM_FILENAME
    > "$OUTPUT_FILE"
    cmd="samtools depth $BAM_FILE | awk '{sum+=\$3} END { print \"Average = \",sum/NR}' > $OUTPUT_FILE"
    sbatch -J $BAM_FILENAME --time=36:00:00 --mem=5Gb --output=$OUTPUT_dir/log/$BAM_FILENAME.out --error=$OUTPUT_dir/log/$BAM_FILENAME.err --account=pi-spott  -p spott -c 1 --wrap="$cmd"

    done
    
