dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
args <- commandArgs(trailingOnly = T)
jobid <- as.numeric(args[1])
ct <- as.character(args[2])

library(stringr)
library(dplyr)
library(tidyverse)
library(data.table)
"%&%" <- function(a, b) paste0(a, b)

new_colnames <- c("Feature_ID", "rsID", "Chromosome", "SNP_position", "Ref_allele",
                  "Alt_allele", "Freq", "HWE_Chisquare", "Imp_quality", "Log10_BH_Q",
                  "Chisquare", "Effect_size", "Delta", "Phi", "Overdispersion",
                  "SNP_id_region", "Num_feature_SNPs", "Num_tested_SNPs", "Num_iterations_null", "Num_iterations_alt",
                  "Random_ties", "Log_likelihood_null", "Convergence_status", "Sq_corr_fSNPs", "Sq_corr_rSNP")


## helper utilities
result.rootDir <- '~/yuzhao1/work/atac_gca2024/19rasqual/4calculation_all_peaks/4ctSpecific_compliedOutput/'
rasqual.dir <- result.rootDir %&% 'result/'
perm1.dir <- result.rootDir %&% 'perm1/'
perm2.dir <- result.rootDir %&% 'perm2/'
perm3.dir <- result.rootDir %&% 'perm3/'
perm4.dir <- result.rootDir %&% 'perm4/'
perm5.dir <- result.rootDir %&% 'perm5/'
result.dirs <- c(rasqual.dir, perm1.dir, perm2.dir, perm3.dir, perm4.dir, perm5.dir)

## select which iteration to work on
result.dir <- result.dirs[[jobid]]
cat(paste0('reading...job ', jobid, ' ', ct, '\n'))

df_list <- list()
for (chr in paste0('chr', c(1:22, 'X'))) {
  df_list[[chr]] <- read.table(result.dir %&% ct %&% '/' %&% chr, sep = '\t')
}
df <- as_tibble(do.call(rbind, df_list))
colnames(df) <- new_colnames

## save results
cat('saving...\n')
saveRDS(df, paste0(result.dir, ct, '.rds'))
cat('done!\n')






