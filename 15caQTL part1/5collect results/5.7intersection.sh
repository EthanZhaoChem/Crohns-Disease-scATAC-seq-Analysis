bedfilesDir=~/yuzhao1/work/atac_gca2024/19rasqual/5results/snps/hg19/bed_lifted_extended
outDir=~/yuzhao1/work/atac_gca2024/19rasqual/5results/snps/hg19/intersection

bedtools intersect -a $bedfilesDir/all_30385.bed -b $bedfilesDir/snps_nt2017_4312.bed > $outDir/caQTL_nt2017SNP_overlap.bed


