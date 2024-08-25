dir_output=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.2ligate
dir_phase=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.1phase
dir_log=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.2ligate_log
sampleIDs_taskfile=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.1sampleIDs_run
ligate_submission=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.2ligate.batch

for sampleID in `cat $sampleIDs_taskfile | awk '{print $1}' | sort | uniq`;
do

output_cell=$dir_output/${sampleID}
if [ ! -d $output_cell ]
then
    mkdir $output_cell
fi

outputLog_cell=$dir_log/${sampleID}
if [ ! -d $outputLog_cell ]
then
    mkdir $outputLog_cell
fi

phase_output_cell=$dir_phase/${sampleID}

    for chrom in {1..22} X;
    do
        sbatch -J "${sampleID}/${chrom}"  $ligate_submission $sampleID $chrom  $output_cell $phase_output_cell 
    done 
    
done
