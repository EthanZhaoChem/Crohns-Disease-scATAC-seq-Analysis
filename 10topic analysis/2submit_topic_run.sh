module load R/4.1.0
module load proj
module load gsl

# manual1
flag_continue=1

all_lineages="union"
# all_lineages="epithelial immune stromal"
# all_lineages="epithelial immune"
# all_lineages="immune"
# all_lineages="stromal"

# nTopic_sequence=$(echo {5..30..5})
nTopic_sequence=$(echo {35..50..5})
# nTopic_sequence=$(echo {15..25..5})
# nTopic_sequence="15 20 25 30"
# nTopic_sequence="30"

# manual2
nIteration1=80  
nIteration2=20
iterations_before=650
iterations_add=100

# pretty fixed
dir_rds=/project/gca/yuzhao1/work/atac_gca2024/13fasttopic/rds
dir_log=/project/gca/yuzhao1/work/atac_gca2024/13fasttopic/log
script_submission=/project/gca/yuzhao1/work/atac_gca2024/13fasttopic/2topic_run.R


if [ "$flag_continue" -eq 0 ]; then
  nIteration_aiming=$((nIteration1 + nIteration2))
fi

if [ "$flag_continue" -eq 1 ]; then
  nIteration_aiming=$((iterations_before + iterations_add))
fi

# it will run as long as there is not a file suffixed with 'converged.rds'
for lineage in $all_lineages; do
    for nTopic in $nTopic_sequence; do
        
        oldFile="fit_${lineage}_sub100_k${nTopic}_${iterations_before}iterations.rds" 

        # previous iteration not finished, skip
        if [[ $flag_continue == 1 && ! -f "$dir_rds/$oldFile" ]]; then
            continue
        fi
        
        # converged or new file already existed, skip
        convergedFile="fit_${lineage}_sub100_k${nTopic}_converged.rds"
        newFile="fit_${lineage}_sub100_k${nTopic}_${nIteration_aiming}iterations.rds" 
        if [[ ! -f "$dir_rds/$convergedFile" && ! -f "$dir_rds/$newFile" ]]
        then
        cmd="Rscript --vanilla $script_submission $lineage $nTopic $nIteration1 $nIteration2 $flag_continue $iterations_before $iterations_add"
        sbatch -J ${lineage}_${nTopic}_${nIteration_aiming} --time=36:00:00 --mem=180G --output=$dir_log/${lineage}_${nTopic}_$nIteration_aiming.out --error=$dir_log/${lineage}_${nTopic}_$nIteration_aiming.err --account=pi-spott  -p caslake -c 24 --wrap="$cmd"
        fi
    done
done












