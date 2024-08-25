dyn.load('/software/geos-3.9.1-el8-x86_64/lib64/libgeos_c.so')
library(stringr)
library(tidyverse)
library(data.table)
library(rasqualTools)

randomize <- function(x,g=NULL){
  if(is.null(g)){
    n=ncol(x);
    t(apply(x,1,function(xx){xx[order(runif(n))]}))
  }else{
    for(i in unique(g)){
      x[,g==i]=randomize(x[,g==i,drop=F])
    }
    x
  }
}

rasqualMakeCovariates <- function(counts, size_factors) {
  
  #Map parameters to Natsuhiko's variables
  Y = counts
  K = size_factors
  n=ncol(Y)
  
  # fpm calculation
  fpkm=t(t(Y/K+1)/apply(Y/K,2,sum))*1e6 #  /len*1e9
  
  # Singular value decomposition
  fpkm.svd   = svd((log(fpkm)-apply(log(fpkm),1,mean))/apply(log(fpkm),1,sd))
  fpkm.svd.r = svd(randomize((log(fpkm)-apply(log(fpkm),1,mean))/apply(log(fpkm),1,sd)))
  
  # Covariate selection
  sf=log(apply(Y,2,sum))
  covs=fpkm.svd$v[,1:sum(fpkm.svd$d[-n]>fpkm.svd.r$d[-n])]
  if(cor(sf,covs[,1])^2<0.9){covs=cbind(sf, covs)}
  
  # Write covariates
  return(covs)
}


























