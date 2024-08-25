library(stringr)
sink("~/yuzhao1/work/final_GCAatac/18glimpse/1.4ref_chunks.sh", append = F)

cat('cd /home/yuzhao1/yuzhao1/software/glimpse2\n')


for(chr in c(1:22, 'X')){
  cmd = paste0(
    'GLIMPSE2_chunk_static --input ref_panel/1000GP.chr',chr,
    '.sites.vcf.gz --region chr',chr,
    ' --sequential --output chunks/chunks.chr',chr,
    '.txt --map github/maps/genetic_maps.b38/chr',chr,
    '.b38.gmap.gz'
  )
  
  cat(sep = '', 'sbatch --time=36:00:00 --mem=40Gb --output=log/1.4_', chr, '.out --error=log/1.4_', chr, ".err --account=pi-spott  -p caslake -c 1 --wrap='", cmd, "'\n")
}


sink()





