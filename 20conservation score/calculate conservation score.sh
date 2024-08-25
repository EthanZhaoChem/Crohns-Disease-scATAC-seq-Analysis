#!/bin/bash
bigWigFile="/project/gca/yuzhao1/resource/UCSC/hg38.phastCons30way.bw"
bedFile="/project/gca/yuzhao1/work/atac_gca2024/23conservation/union_peakset.bed"
outputFile="/project/gca/yuzhao1/work/atac_gca2024/23conservation/phastCons_scores.txt"
ucsc_tool_folder="/project/gca/yuzhao1/resource/UCSC/"

while read -r line; do
  chrom=$(echo "$line" | cut -f1)
  start=$(echo "$line" | cut -f2)
  end=$(echo "$line" | cut -f3)
  score=$($ucsc_tool_folder/bigWigSummary $bigWigFile $chrom $start $end 1)
  echo -e "$chrom\t$start\t$end\t$score" >> $outputFile
done < $bedFile