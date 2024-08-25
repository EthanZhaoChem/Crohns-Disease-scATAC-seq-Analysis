software_path=/home/yuzhao1/yuzhao1/software/liftover
hg19_chr_size_path=/home/yuzhao1/yuzhao1/resource/UCSC/hg19.chrom.sizes

# need to modify
bedfilesRootDir=~/yuzhao1/work/atac_gca2024/19rasqual/5results/snps/hg38
DirRoot=~/yuzhao1/work/atac_gca2024/19rasqual/5results/snps/hg19

# fixed
bedfilesNewDir=$DirRoot/bed_lifted
bedfilesNewDir2=$DirRoot/bed_trash
bedfilesExtendDir=$DirRoot/bed_lifted_extended

if [ ! -d $DirRoot ]
then
    mkdir $DirRoot
fi


if [ ! -d $bedfilesNewDir ]
then
    mkdir $bedfilesNewDir
fi

if [ ! -d $bedfilesNewDir2 ]
then
    mkdir $bedfilesNewDir2
fi

if [ ! -d $bedfilesExtendDir ]
then
    mkdir $bedfilesExtendDir
fi

names=`ls  $bedfilesRootDir | cut -f 1 -d '.'`
for name in $names
do
    bedname=`echo $name | awk '{print $1}'`
    echo $bedname
    ${software_path}/liftOver $bedfilesRootDir/$bedname.bed ${software_path}/hg38ToHg19.over.chain $bedfilesNewDir/$bedname.bed $bedfilesNewDir2/${bedname}_unlifted.bed
    bedtools sort -i $bedfilesNewDir/$bedname.bed
    
    bedtools slop -i $bedfilesNewDir/$bedname.bed -g $hg19_chr_size_path -b 1000 > $bedfilesExtendDir/$bedname.bed

done




