dir_output=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.3eagle
dir_ligate=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.2ligate
dir_log=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.3eagle_log
sampleIDs_taskfile=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.1sampleIDs_run
script_submission=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.3eagle.batch

    
for sampleID in `cat $sampleIDs_taskfile | awk '{print $1}' | sort | uniq`;
do

    ligate_output_cell=$dir_ligate/${sampleID}
    output_cell=$dir_output/${sampleID}
    outputLog_cell=$dir_log/${sampleID}
    
    if [ ! -d $output_cell ]
    then
        mkdir $output_cell
    fi
    
    if [ ! -d $outputLog_cell ]
    then
        mkdir $outputLog_cell
    fi
    
    for chrom in {1..22} X;
    do
        sbatch -J "${sampleID}/${chrom}"  $script_submission $sampleID $chrom $output_cell $ligate_output_cell
    done 
    
done


