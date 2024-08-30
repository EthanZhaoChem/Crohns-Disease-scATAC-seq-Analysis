software_path=/project/gca/yuzhao1/software/liftover

# need to modify
bedfilesRootDir=/project/gca/yuzhao1/work/final_GCAatac/0gwas/tmp

${software_path}/liftOver $bedfilesRootDir/hg19.bed ${software_path}/hg19ToHg38.over.chain $bedfilesRootDir/hg38.bed $bedfilesRootDir/trash.bed
bedtools sort -i $bedfilesRootDir/hg38.bed





