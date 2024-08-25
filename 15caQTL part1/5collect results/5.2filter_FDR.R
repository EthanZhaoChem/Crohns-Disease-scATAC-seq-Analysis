dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
args <- commandArgs(trailingOnly = T)
ct <- as.character(args[1])
total_subDirs <- 4 # number of all rasqual runs (1 + #perm)

library(stringr)
library(dplyr)
library(tidyverse)
library(data.table)
"%&%" <- function(a, b) paste0(a, b)

result.rootDir <- '~/yuzhao1/work/atac_gca2024/19rasqual/4calculation_all_peaks/4ctSpecific_compliedOutput/'
rasqual.dir <- result.rootDir %&% 'result/'
perm1.dir <- result.rootDir %&% 'perm1/'
perm2.dir <- result.rootDir %&% 'perm2/'
perm3.dir <- result.rootDir %&% 'perm3/'
perm4.dir <- result.rootDir %&% 'perm4/'
perm5.dir <- result.rootDir %&% 'perm5/'
result.dirs <- c(rasqual.dir, perm1.dir, perm2.dir, perm3.dir, perm4.dir, perm5.dir)

## presets
out.calculated.rootDir <- '~/yuzhao1/work/atac_gca2024/19rasqual/5collect_results/5results/allPeaks_noDiseaseCo/'
dir.plots <- '~/yuzhao1/work/atac_gca2024/19rasqual/5collect_results/5results/allPeaks_noDiseaseCo/plots/'


## helper functions
# q1 : real lead Q-value vector for all peaks from RASQUAL
# q0 : permutated Q-value vector
# alpha : FDR threshold
# This function returns the P-value threshold corresponding to FDR=alpha.
# see: https://github.com/natsuhiko/rasqual/issues/21
getFDR <- function(q1,
                   q0,
                   alpha = 0.1,
                   z = NULL,
                   subset = NULL) {
  if (is.null(z)) {
    a = 0
    for (itr in 1:10) {
      a = getFDR(q1, q0, alpha, rev(a + 0:100 / 100 ^ itr), subset)
    }
    a
  } else{
    if (!is.null(subset)) {
      q1 = q1[subset]
      q0 = q0[subset]
    }
    q1 = q1[!is.na(q1)]
    q0 = q0[!is.na(q0)]
    x = NULL
    
    for (i in z) {
      x = c(x, sum(q0 < i) / length(q0) / (sum(q1 < i) / length(q1)))
    }
    
    max(c(0, z[x < alpha]), na.rm = T)
  }
}


ct.tmpList <- list()
ct.tmpList.sorted <- list()
for (i in 1:total_subDirs) {
  # only converted the cols that are used here
  cols_to_convert <- c('Log10_BH_Q', 'Chisquare', 'Effect_size')
  tmp_raw <- readRDS(paste0(result.dirs[[i]], ct, '.rds'))
  tmp_raw <- tmp_raw[tmp_raw$Convergence_status == 0, ] # only use converged results
  xx <- tmp_raw %>%
    mutate(across(all_of(cols_to_convert), as.numeric))
  ct.tmpList[[i]] <- xx %>% filter(rsID != "SKIPPED") %>% mutate(P = pchisq(Chisquare, 1, lower = F)) # https://github.com/natsuhiko/rasqual/issues/22
  ct.tmpList.sorted[[i]] <- ct.tmpList[[i]] %>% arrange(P) %>% mutate(q = 10^(Log10_BH_Q))
}

# just to make it clear that the first df is results without perturbation
rm(ct.tmpList, tmp_raw, xx)
gc()
names(ct.tmpList.sorted) <- c('results', 'perm1', 'perm2', 'perm3')

# use lead variants to calculate FDR (already sorted, just remove duplicated regions)
results_lead  = ct.tmpList.sorted[['results']][!duplicated(ct.tmpList.sorted[['results']]$Feature_ID),]
perm1_lead  = ct.tmpList.sorted[['perm1']][!duplicated(ct.tmpList.sorted[['perm1']]$Feature_ID),]
perm2_lead  = ct.tmpList.sorted[['perm2']][!duplicated(ct.tmpList.sorted[['perm2']]$Feature_ID),]
perm3_lead  = ct.tmpList.sorted[['perm3']][!duplicated(ct.tmpList.sorted[['perm3']]$Feature_ID),]


# test hypothesis
q0 <- c(perm1_lead$q, perm2_lead$q, perm3_lead$q)
q1 <- results_lead$q

# #Plot qq plot, Cap q-values at min.q for drawing purposes
png(paste0(dir.plots, ct, '.png'), height = 1500, width = 2000, res=300)
null.q <- perm1_lead$q
obs.q <- q1
min.q <- 1e-16
obs.q[obs.q < min.q]  <- min.q
null.q[null.q < min.q] <- min.q
n.test <- nrow(results_lead)
qqplot(-log10(null.q), -log10(obs.q), xlab = "Permuted Log10(q-value)", ylab = "Log10(q-value) caQTL")
abline(a=0, b=1)
dev.off()


#Mark significant caQTLs (10% FDR), True = significant QTLs
fdr_value <- getFDR(q1, q0, 0.1)

flag_fdr01_all <- ct.tmpList.sorted$results$q < fdr_value
flag_fdr01_lead <- results_lead$q < fdr_value

table(flag_fdr01_all)
table(flag_fdr01_lead)

# save all peaks all snps with FDR label
ct.tmpList.sorted[['results']]$keep = flag_fdr01_all
saveRDS(ct.tmpList.sorted[['results']], paste0(out.calculated.rootDir, ct, '.rds'))

# save lead snp results
xx <- ct.tmpList.sorted[['results']][!duplicated(ct.tmpList.sorted[['results']]$Feature_ID),]
saveRDS(xx, paste0(out.calculated.rootDir, ct, '_lead.rds'))

# subset for rows passing FDR cutoff
FDR_results <- ct.tmpList.sorted[['results']] %>% filter(keep == TRUE)

filename_csv <- paste0(out.calculated.rootDir, ct, '_Perm_three_times_FDR_0.1.csv')
filename_rds <- paste0(out.calculated.rootDir, ct, '_Perm_three_times_FDR_0.1.rds')
saveRDS(FDR_results, filename_rds)
write.csv(FDR_results, filename_csv)




