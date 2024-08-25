# this uses an interactive job
split_submission=/project/gca/yuzhao1/work/final_GCAatac/18glimpse/1.5split.batch

mkdir /project/gca/yuzhao1/software/glimpse2/ref_panel/split


for chrom in {1..22} X;
do
    sbatch -J "1.5_${chrom}" $split_submission ${chrom}

done

