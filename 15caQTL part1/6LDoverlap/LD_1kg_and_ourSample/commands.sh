sinteractive --partition=caslake --account=pi-onibasu --ntasks-per-node=1 --cpus-per-task=8 --mem=100G --time=36:00:00

# copy genotype data: 1kg and our samples
cd /home/yuzhao1/yuzhao1/work/atac_gca2024/19rasqual/8LDoverlap/ld_1kg_and_ourSample
cp ../../3pca/plink_analysis/all_sample_1KGP.* .

# compute ld 
module load plink
plink --bfile all_sample_1KGP --r2 --ld-window-r2 0.6 --ld-window 99999 --ld-window-kb 1000 --out ld_0.6/all_sample_1KGP_ld --threads 8
plink --bfile all_sample_1KGP --r2 --ld-window-r2 0.8 --ld-window 99999 --ld-window-kb 1000 --out ld_0.8/all_sample_1KGP_ld --threads 8

