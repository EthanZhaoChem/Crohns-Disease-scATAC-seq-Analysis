library(stringr)
sink("~/yuzhao1/work/final_GCAatac/18glimpse/1.3ref_extract_sites.sh", append = F)

cat('cd /home/yuzhao1/yuzhao1/software/glimpse2\n')


for(chr in c(1:22, 'X')){
  cmd = paste0(
    'bcftools view -G -Oz -o ref_panel/1000GP.chr',chr,
    '.sites.vcf.gz ref_panel/1000GP.chr',chr,'.bcf\n',
    'bcftools index -f ref_panel/1000GP.chr',chr,'.sites.vcf.gz'
  )
  
  cat(sep = '', 'sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_', chr, '.out --error=log/1.3_', chr, ".err --account=pi-spott  -p caslake -c 1 --wrap='", cmd, "'\n")
}


sink()





