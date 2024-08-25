module load R/4.1.0
script_submission=/project/gca/yuzhao1/work/atac_gca2024/19rasqual//5collect_results/5.2filter_FDR.R
celltype_file=/project/gca/yuzhao1/work/atac_gca2024/19rasqual/00celltypes_filtered.txt

for ct in `cat $celltype_file | awk '{print $1}' | sort -V | uniq`;
do
  cmd="Rscript --vanilla $script_submission $ct"
  sbatch -J "${ct}"  -o "log/${ct}.out" \
  -e "log/${ct}.err" \
  --time=12:00:00 --mem=100G --account=pi-spott -p caslake -c 1 --wrap="$cmd"
done

