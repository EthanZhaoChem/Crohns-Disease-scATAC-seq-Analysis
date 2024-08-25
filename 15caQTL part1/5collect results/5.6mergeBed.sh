# make sure the bed files are not overlapping

# need to modify
bedfilesRootDir=~/yuzhao1/work/atac_gca2024/19rasqual/5results/snps/hg19/bed_lifted_extended


names=`ls  $bedfilesRootDir | cut -f 1 -d '.'`
for name in $names
do
    bedname=`echo $name | awk '{print $1}'`
    echo $bedname
    bedtools merge -i $bedfilesRootDir/$bedname.bed > $bedfilesRootDir/temp.bed
    mv $bedfilesRootDir/temp.bed $bedfilesRootDir/$bedname.bed
done




