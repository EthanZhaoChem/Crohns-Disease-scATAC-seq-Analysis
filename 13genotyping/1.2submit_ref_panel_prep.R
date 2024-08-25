library(stringr)
sink("~/yuzhao1/work/final_GCAatac/18glimpse/1.2submit_ref_panel_prep.sh", append = F)

cat('cd /home/yuzhao1/yuzhao1/software/glimpse2\n')


for(chr in 1:22){
  cmd = paste0(
    'bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr',
    chr,
    '.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr', chr,'.bcf\n',
    'bcftools index -f ref_panel/1000GP.chr', chr, '.bcf --threads 4'
  )
  
  cat(sep = '', 'sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_', chr, '.out --error=log/1.2_', chr, ".err --account=pi-spott  -p caslake -c 4 --wrap='", cmd, "'\n")
}

chr <- 'X'
cmd = paste0(
  'bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr',
  chr,
  '.filtered.eagle2-phased.v2.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr', chr,'.bcf\n',
  'bcftools index -f ref_panel/1000GP.chr', chr, '.bcf --threads 4'
)

cat(sep = '', 'sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_', chr, '.out --error=log/1.2_', chr, ".err --account=pi-spott  -p caslake -c 4 --wrap='", cmd, "'\n")


sink()





