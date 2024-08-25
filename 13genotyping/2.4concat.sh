dir_concat=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.4concat
dir_eagle=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.3eagle
dir_log=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.4concat_log
sampleIDs_taskfile=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.1sampleIDs_run
script_submission=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.4concat.batch


for sampleID in `cat $sampleIDs_taskfile | awk '{print $1}' | sort | uniq`;
do

    eagle_output_cell=$dir_eagle/${sampleID}
    output_cell=$dir_concat
    outputLog_cell=$dir_log
    
    if [ ! -d $output_cell ]
    then
        mkdir $output_cell
    fi
    
    if [ ! -d $outputLog_cell ]
    then
        mkdir $outputLog_cell
    fi
    

    sbatch -J "${sampleID}"  $script_submission $sampleID $output_cell $eagle_output_cell
done




