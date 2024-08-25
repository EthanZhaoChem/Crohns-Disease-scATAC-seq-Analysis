module load R/4.1.0
script_submission=/project/gca/yuzhao1/work/atac_gca2024/19rasqual//5collect_results/5.1collect_raw.R
celltype_file=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt

for ct in `cat $celltype_file | awk '{print $1}' | sort -V | uniq`;
do

  # for jobid in 1 2 3
  for jobid in 4
  do
    cmd="Rscript --vanilla $script_submission $jobid $ct"
    sbatch -J "job${jobid}_${ct}"  -o "log/job${jobid}_${ct}.out" \
    -e "log/job${jobid}_${ct}.err" \
    --time=12:00:00 --mem=30G --account=pi-onibasu -p caslake -c 1 --wrap="$cmd"
  done
  
done

