args <- commandArgs(trailingOnly = TRUE) 
sample1 <- args[1]
sample2 <- args[2]

library(VariantAnnotation)
resultPath <- paste0("~/yuzhao1/work/final_GCAatac/18glimpse/4.1sample_concordance/", sample1, '_', sample2, ".rds")
resultPath_txt <- paste0("~/yuzhao1/work/final_GCAatac/18glimpse/4.1sample_concordance/", sample1, '_', sample2, ".txt")


# Load the VCF file
vcfFilePath1 <- paste0("~/yuzhao1/work/final_GCAatac/18glimpse/2.5reheader/", sample1, ".vcf.gz")
vcfFilePath2 <- paste0("~/yuzhao1/work/final_GCAatac/18glimpse/2.5reheader/", sample2, ".vcf.gz")

vcf1 <- readVcf(vcfFilePath1, genome="hg38") 
vcf2 <- readVcf(vcfFilePath2, genome="hg38") 

# Extract genotypes for "sample1" and "sample2"
genotypes1 <- geno(vcf1)$GT
genotypes_sample1 <- genotypes1
genotypes2 <- geno(vcf2)$GT
genotypes_sample2 <- genotypes2

# Function to convert genotype strings to numeric codes (0, 1, 2)
convert_genotypes <- function(genotype_str) {
  sapply(genotype_str, function(g) {
    if (g == "0/0" || g == "0|0") {
      return(0)
    } else if (g == "0/1" || g == "1/0" || g == "0|1" || g == "1|0") {
      return(1)
    } else if (g == "1/1" || g == "1|1") {
      return(2)
    } else {
      return(NA) # for cases with missing or complex genotypes
    }
  })
}

# Convert genotype strings to numeric codes
codes_sample1 <- convert_genotypes(genotypes_sample1)
codes_sample2 <- convert_genotypes(genotypes_sample2)

# Calculate concordance
concordance_table <- table(codes_sample1, codes_sample2)

# Calculate metrics
a <- concordance_table[1,1]
e <- concordance_table[2,2]
i <- concordance_table[3,3]
total <- sum(concordance_table)

overall_concordance <- (a + e + i) / total
nrc <- (e + i) / (total - concordance_table[1,1])



# Perform calculation
concordance_results <- list(nonreference_concordance = nrc, 
                            overall_concordance = overall_concordance, 
                            concordance_table = concordance_table)

# save results
saveRDS(concordance_results, resultPath)

# print
sink(resultPath_txt, append = F)
print(concordance_results)
sink()














