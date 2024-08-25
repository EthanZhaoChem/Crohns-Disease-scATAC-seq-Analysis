args <- commandArgs(trailingOnly = T)
trait2 <- args[1]
annot_cell <- args[2]
results_cell <- args[3]
post_process_cell <- args[4]

library(ggplot2)
library(cowplot)
library(tidyverse)
library(plyr)
library(dplyr)
library(stringr)
library(rmeta)
library(data.table)
source('~/yuzhao1/scripts/plot.R')
source('~/yuzhao1/work/final_GCArna/scripts/gca_colors.R')

########################### part 1 ##############################
# check parameter M: all snps used for heritability analysis
Mref <- 0
for (chr in 1:22){
  snps_tmp <- readLines(paste0('~/spott/yuzhao1/ldsc/minimum_required/1000G_Phase3_ldscores/LDscore.', chr, '.l2.M_5_50'))
  snps_tmp <- as.numeric(snps_tmp)
  Mref <- Mref + snps_tmp
}
# after calculation, it is 5961159


########################### part 2 ##############################
annot_names <- list.dirs(annot_cell, recursive = F, full.names = F)


########################### part 3 ##############################
get_sd_annot = function(cell_path, annot_index = 1){
  num = 0
  den = 0
  ll <- list.files(cell_path, pattern = ".annot.gz")
  for(m in 1:length(ll)){
    dat <- data.frame(fread(cmd=paste0("zcat ", cell_path, "/", ll[m])))
    num = num  + (nrow(dat)-1) * var(dat[,4+annot_index])
    den = den + (nrow(dat)-1)
    rm(dat)
  }
  estd_sd_annot = sqrt(num/den)
  return(estd_sd_annot)
}

########################### part 4 ##############################
index_in_results=1 
base_index = NULL
if(is.null(base_index)){base_index = index_in_results}

final_df = c()
for(aa in 1:length(annot_names)){
  cat(aa, '\n')
  cell_path = paste0(annot_cell, "/", annot_names[aa])
  sd_annot1=get_sd_annot(cell_path, annot_index=base_index)
  result.file=paste0(results_cell, "/", annot_names[aa], "/", trait2, ".part_delete")
  new_table=read.table(result.file,header=F)
  logfile = paste(results_cell, "/", annot_names[aa], "/", trait2,".log", sep="")
  log = read.table(logfile,h=F,fill=T)
  h2g = as.numeric(as.character(log[which(log$V4=="h2:"),5]))
  coef1=sd_annot1*Mref/h2g
  sc=c()
  for(i in 1:dim(new_table)[1]){
    tau1=as.numeric(new_table[i,base_index])
    taus1=tau1*coef1
    sc=c(sc,taus1)
  }
  mean_sc=mean(sc)
  se_sc=sqrt(199**2/200*var(sc))
  p_sc=pnorm(abs(mean_sc/se_sc), 0, 1, lower.tail=F)*2
  
  result2.file=paste0(results_cell, "/", annot_names[aa], "/", trait2, ".results")
  res2=read.table(result2.file,header=T)
  myenr    =  res2$Enrichment[base_index] #step1
  myenr_sd = res2$Enrichment_std_error[base_index] #step3
  myenr_p = res2$Enrichment_p[base_index]
  final_df = rbind(final_df, c(mean_sc, se_sc, p_sc, myenr, myenr_sd, myenr_p))
}

rownames(final_df) = annot_names
colnames(final_df) = c("tau-star", "se-tau-star", "p-tau-star", "E", "se(E)", "p(E)")
write.table(final_df, file = paste0(post_process_cell, '/', trait2, "_ldsc_postprocess.txt"), quote=F, sep = "\t")



                            
                                