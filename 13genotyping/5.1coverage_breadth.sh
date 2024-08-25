#!/bin/bash

cd /project/gca/yuzhao1/work/final_GCAatac/18glimpse

# Specify the directory containing your BAM files
BAM_DIR="/project/gca/yuzhao1/work/final_GCAatac/18glimpse/3.1patient_bam_merge"

# Output file to store the results
OUTPUT_dir="/project/gca/yuzhao1/work/final_GCAatac/18glimpse/5.1coverage_breadth/"

# Minimum coverage depth
MIN_COVERAGE_DEPTH=1

# Iterate over each BAM file in the directory
for BAM_FILE in "$BAM_DIR"/*.bam; do
    # Extract the filename without the directory path
    BAM_FILENAME=$(basename "$BAM_FILE")
    OUTPUT_FILE=$OUTPUT_dir/$BAM_FILENAME
    > "$OUTPUT_FILE"
    cmd="samtools mpileup $BAM_FILE | awk -v X=$MIN_COVERAGE_DEPTH '\$4>=X' | wc -l > $OUTPUT_FILE"
    sbatch --time=36:00:00 --mem=20Gb --output=$OUTPUT_dir/log/$BAM_FILENAME.out --error=$OUTPUT_dir/log/$BAM_FILENAME.err --account=pi-spott  -p caslake -c 1 --wrap="$cmd"

    done

