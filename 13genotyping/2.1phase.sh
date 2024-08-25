dir_cellranger=/project/spott/yuzhao1/GEO_GCAatac/cellranger_output
dir_glimpse=/project/gca/yuzhao1/software/glimpse2
dir_output=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.1phase
dir_log=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.1phase_log
sampleIDs_taskfile=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.1sampleIDs_run
phase_submission=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.1phase.batch

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

    for chrom in {1..22} X;
    do
        sbatch -J "${sampleID}/${chrom}"  $phase_submission $dir_glimpse $dir_cellranger $output_cell $outputLog_cell $sampleID $chrom 
    done 
    
done
