cd /home/yuzhao1/yuzhao1/software/glimpse2
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_1.out --error=log/1.2_1.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr1.bcf
bcftools index -f ref_panel/1000GP.chr1.bcf --threads 4'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_2.out --error=log/1.2_2.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr2.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr2.bcf
bcftools index -f ref_panel/1000GP.chr2.bcf --threads 4'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_3.out --error=log/1.2_3.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr3.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr3.bcf
bcftools index -f ref_panel/1000GP.chr3.bcf --threads 4'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_4.out --error=log/1.2_4.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr4.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr4.bcf
bcftools index -f ref_panel/1000GP.chr4.bcf --threads 4'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_5.out --error=log/1.2_5.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr5.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr5.bcf
bcftools index -f ref_panel/1000GP.chr5.bcf --threads 4'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_6.out --error=log/1.2_6.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr6.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr6.bcf
bcftools index -f ref_panel/1000GP.chr6.bcf --threads 4'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_7.out --error=log/1.2_7.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr7.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr7.bcf
bcftools index -f ref_panel/1000GP.chr7.bcf --threads 4'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_8.out --error=log/1.2_8.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr8.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr8.bcf
bcftools index -f ref_panel/1000GP.chr8.bcf --threads 4'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_9.out --error=log/1.2_9.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr9.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr9.bcf
bcftools index -f ref_panel/1000GP.chr9.bcf --threads 4'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_10.out --error=log/1.2_10.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr10.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr10.bcf
bcftools index -f ref_panel/1000GP.chr10.bcf --threads 4'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_11.out --error=log/1.2_11.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr11.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr11.bcf
bcftools index -f ref_panel/1000GP.chr11.bcf --threads 4'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_12.out --error=log/1.2_12.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr12.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr12.bcf
bcftools index -f ref_panel/1000GP.chr12.bcf --threads 4'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_13.out --error=log/1.2_13.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr13.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr13.bcf
bcftools index -f ref_panel/1000GP.chr13.bcf --threads 4'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_14.out --error=log/1.2_14.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr14.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr14.bcf
bcftools index -f ref_panel/1000GP.chr14.bcf --threads 4'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_15.out --error=log/1.2_15.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr15.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr15.bcf
bcftools index -f ref_panel/1000GP.chr15.bcf --threads 4'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_16.out --error=log/1.2_16.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr16.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr16.bcf
bcftools index -f ref_panel/1000GP.chr16.bcf --threads 4'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_17.out --error=log/1.2_17.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr17.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr17.bcf
bcftools index -f ref_panel/1000GP.chr17.bcf --threads 4'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_18.out --error=log/1.2_18.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr18.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr18.bcf
bcftools index -f ref_panel/1000GP.chr18.bcf --threads 4'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_19.out --error=log/1.2_19.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr19.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr19.bcf
bcftools index -f ref_panel/1000GP.chr19.bcf --threads 4'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_20.out --error=log/1.2_20.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr20.bcf
bcftools index -f ref_panel/1000GP.chr20.bcf --threads 4'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_21.out --error=log/1.2_21.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr21.bcf
bcftools index -f ref_panel/1000GP.chr21.bcf --threads 4'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_22.out --error=log/1.2_22.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chr22.bcf
bcftools index -f ref_panel/1000GP.chr22.bcf --threads 4'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.2_X.out --error=log/1.2_X.err --account=pi-spott  -p caslake -c 4 --wrap='bcftools norm -m -any ref_downloaded/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps  --threads 4 -Ob -o ref_panel/1000GP.chrX.bcf
bcftools index -f ref_panel/1000GP.chrX.bcf --threads 4'
