library(stringr)
dir_fastq <- '/project/spott/yuzhao1/GEO_GCAatac/fastq_all_20240211'
df <- read.csv('~/spott/yuzhao1/GEO_GCAatac/fastq_metadata_20240211.csv', row.names = 1)

sink("~/yuzhao1/work/final_GCAatac/1upstream/1submit_atac.sh", append = F)
cat('cd ~/spott/yuzhao1/GEO_GCAatac/cellranger_output\n')


for(sample in unique(df$sample)){
  tmp <- unique(df[df$sample == sample, 'seq_library'])
  if(length(tmp) > 1){
    tmp2 <- paste0(tmp[[1]], ',', tmp[[2]])
  }
  else{
    tmp2 <- tmp
  }
  
  cmd=paste0('/project2/gca/software/cellranger-atac-2.0.0/cellranger-atac count --id=', sample,
             ' --reference=/project2/gca/software/cellranger-atac-2.0.0/reference-data/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/',
             ' --fastqs=', dir_fastq,
             ' --sample=', tmp2,
             ' --localcores=4 --localmem=56')
  
  cat(sep = '', 'sbatch --time=36:00:00 --mem=64Gb --output=00log/', sample, '.out --error=00log/', sample, ".err --account=pi-spott  -p caslake -c 4 --wrap='", cmd, "'\n")
}

sink()





