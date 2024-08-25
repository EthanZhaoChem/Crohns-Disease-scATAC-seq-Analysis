library(stringr)
library(ArchR)
library(Seurat)
library(TFBSTools)
library(JASPAR2020)
library(seqLogo)
library(ggseqlogo)

#############################################################################
# plot a logo from jaspar database
pfm2ppm <- function(pfm) {
  ppm <- apply(pfm, 2, function(col) {
    col / sum(col)} ) 
  return(ppm)
}

ppm2pwm <- function(ppm){
  pwm <- log2(ppm/0.25)
  pwm[is.infinite(pwm)] <- 0
  return(pwm)
}

# # example usage
# ppm_set1 <- getMatrixSet(x = JASPAR2020, opts = list(all_versions = FALSE, species = 'Homo sapiens', collection="CORE", matrixtype="PFM")) 
# ppm_set2 <- getMatrixSet(x = JASPAR2020, opts = list(all_versions = FALSE, species = 'Homo sapiens', collection="UNVALIDATED", matrixtype="PFM")) 
# ppm_set <- c(ppm_set1, ppm_set2)
# 
# 
# p <- pfm2ppm(ppm_set2$UN0110.1@profileMatrix)
# seqLogo(p, ic.scale=F)
# seqLogo(p, ic.scale=T)




















