software_path=/home/yuzhao1/yuzhao1/software/liftover
hg19_chr_size_path=/home/yuzhao1/yuzhao1/resource/UCSC/hg19.chrom.sizes

# need to modify
bedfilesDir1=~/yuzhao1/work/atac_gca2024/19rasqual/5results/snps/hg19/bed_lifted
bedfilesDir2=~/yuzhao1/work/atac_gca2024/19rasqual/5results/snps/hg19/bed_lifted_extended


names=`ls  $bedfilesDir1 | cut -f 1 -d '.'`
for name in $names
do
    bedname=`echo $name | awk '{print $1}'`
    echo $bedname
    
    bedtools slop -i $bedfilesDir1/$bedname.bed -g $hg19_chr_size_path -b 1000 > $bedfilesDir2/$bedname.bed

done




