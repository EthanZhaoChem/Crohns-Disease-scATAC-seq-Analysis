dir_cellranger=/project/spott/yuzhao1/GEO_GCAatac/cellranger_output
dir_output=/project/gca/yuzhao1/work/atac_gca2024/20bam/1bam_perBC/singleBams
dir_log=/project/gca/yuzhao1/work/atac_gca2024/20bam/1bam_perBC/log
sampleIDs_taskfile=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/2.1sampleIDs_all
script_submission=/project/gca/yuzhao1/work/atac_gca2024/20bam/1.2split_barcode.batch

if [ ! -d $dir_log ]
then
    mkdir $dir_log
fi

if [ ! -d $dir_output ]
then
    mkdir $dir_output
fi

for sampleID in `cat $sampleIDs_taskfile | awk '{print $1}' | sort | uniq`;
# for sampleID in "C1-AC-Control";
do

    output_cell=$dir_output/${sampleID}
    if [ ! -d $output_cell ]
    then
        mkdir $output_cell
    fi

    # sbatch -J "${sampleID}" -o "$output_cell" -e "$outputLog_cell" $script_submission  $dir_cellranger $output_cell  $sampleID
    sbatch -J "${sampleID}" -o "$dir_log/$sampleID.out" -e "$dir_log/$sampleID.err" $script_submission  $dir_cellranger $output_cell  $sampleID

done


