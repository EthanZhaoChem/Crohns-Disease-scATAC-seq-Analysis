#!/bin/bash

bigWigFileDir="/project/gca/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1/GroupBigWigs/anno1_individual"
bedFile="/project/gca/yuzhao1/work/atac_gca2024/4peaks/6pseudo_bulk_accessibility/query.bed"
outDir="/project/gca/yuzhao1/work/atac_gca2024/4peaks/6pseudo_bulk_accessibility/results"
outputFile="$outDir/all.bed"
ucsc_tool_folder="/project/gca/yuzhao1/resource/UCSC/"
suffix="-TileSize-25-normMethod-ReadsInTSS-ArchR.bw"


while read -r line; do
  for tmpBW in `ls $bigWigFileDir | awk '{print $1}' | sort | uniq`; # folder name
  do
    bigWigFile=$bigWigFileDir/$tmpBW
    pseudo_bulk="${tmpBW%$suffix}"
    chrom=$(echo "$line" | cut -f1)
    start=$(echo "$line" | cut -f2)
    end=$(echo "$line" | cut -f3)
    score=$($ucsc_tool_folder/bigWigSummary $bigWigFile $chrom $start $end 1)
    echo -e "$chrom\t$start\t$end\t$pseudo_bulk\t$score" >> $outputFile
  done
done < $bedFile

