cd /home/yuzhao1/yuzhao1/software/glimpse2
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_1.out --error=log/1.3_1.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chr1.sites.vcf.gz ref_panel/1000GP.chr1.bcf
bcftools index -f ref_panel/1000GP.chr1.sites.vcf.gz'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_2.out --error=log/1.3_2.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chr2.sites.vcf.gz ref_panel/1000GP.chr2.bcf
bcftools index -f ref_panel/1000GP.chr2.sites.vcf.gz'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_3.out --error=log/1.3_3.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chr3.sites.vcf.gz ref_panel/1000GP.chr3.bcf
bcftools index -f ref_panel/1000GP.chr3.sites.vcf.gz'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_4.out --error=log/1.3_4.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chr4.sites.vcf.gz ref_panel/1000GP.chr4.bcf
bcftools index -f ref_panel/1000GP.chr4.sites.vcf.gz'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_5.out --error=log/1.3_5.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chr5.sites.vcf.gz ref_panel/1000GP.chr5.bcf
bcftools index -f ref_panel/1000GP.chr5.sites.vcf.gz'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_6.out --error=log/1.3_6.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chr6.sites.vcf.gz ref_panel/1000GP.chr6.bcf
bcftools index -f ref_panel/1000GP.chr6.sites.vcf.gz'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_7.out --error=log/1.3_7.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chr7.sites.vcf.gz ref_panel/1000GP.chr7.bcf
bcftools index -f ref_panel/1000GP.chr7.sites.vcf.gz'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_8.out --error=log/1.3_8.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chr8.sites.vcf.gz ref_panel/1000GP.chr8.bcf
bcftools index -f ref_panel/1000GP.chr8.sites.vcf.gz'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_9.out --error=log/1.3_9.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chr9.sites.vcf.gz ref_panel/1000GP.chr9.bcf
bcftools index -f ref_panel/1000GP.chr9.sites.vcf.gz'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_10.out --error=log/1.3_10.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chr10.sites.vcf.gz ref_panel/1000GP.chr10.bcf
bcftools index -f ref_panel/1000GP.chr10.sites.vcf.gz'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_11.out --error=log/1.3_11.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chr11.sites.vcf.gz ref_panel/1000GP.chr11.bcf
bcftools index -f ref_panel/1000GP.chr11.sites.vcf.gz'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_12.out --error=log/1.3_12.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chr12.sites.vcf.gz ref_panel/1000GP.chr12.bcf
bcftools index -f ref_panel/1000GP.chr12.sites.vcf.gz'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_13.out --error=log/1.3_13.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chr13.sites.vcf.gz ref_panel/1000GP.chr13.bcf
bcftools index -f ref_panel/1000GP.chr13.sites.vcf.gz'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_14.out --error=log/1.3_14.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chr14.sites.vcf.gz ref_panel/1000GP.chr14.bcf
bcftools index -f ref_panel/1000GP.chr14.sites.vcf.gz'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_15.out --error=log/1.3_15.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chr15.sites.vcf.gz ref_panel/1000GP.chr15.bcf
bcftools index -f ref_panel/1000GP.chr15.sites.vcf.gz'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_16.out --error=log/1.3_16.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chr16.sites.vcf.gz ref_panel/1000GP.chr16.bcf
bcftools index -f ref_panel/1000GP.chr16.sites.vcf.gz'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_17.out --error=log/1.3_17.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chr17.sites.vcf.gz ref_panel/1000GP.chr17.bcf
bcftools index -f ref_panel/1000GP.chr17.sites.vcf.gz'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_18.out --error=log/1.3_18.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chr18.sites.vcf.gz ref_panel/1000GP.chr18.bcf
bcftools index -f ref_panel/1000GP.chr18.sites.vcf.gz'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_19.out --error=log/1.3_19.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chr19.sites.vcf.gz ref_panel/1000GP.chr19.bcf
bcftools index -f ref_panel/1000GP.chr19.sites.vcf.gz'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_20.out --error=log/1.3_20.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chr20.sites.vcf.gz ref_panel/1000GP.chr20.bcf
bcftools index -f ref_panel/1000GP.chr20.sites.vcf.gz'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_21.out --error=log/1.3_21.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chr21.sites.vcf.gz ref_panel/1000GP.chr21.bcf
bcftools index -f ref_panel/1000GP.chr21.sites.vcf.gz'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_22.out --error=log/1.3_22.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chr22.sites.vcf.gz ref_panel/1000GP.chr22.bcf
bcftools index -f ref_panel/1000GP.chr22.sites.vcf.gz'
sbatch --time=36:00:00 --mem=40Gb --output=log/1.3_X.out --error=log/1.3_X.err --account=pi-spott  -p caslake -c 1 --wrap='bcftools view -G -Oz -o ref_panel/1000GP.chrX.sites.vcf.gz ref_panel/1000GP.chrX.bcf
bcftools index -f ref_panel/1000GP.chrX.sites.vcf.gz'
