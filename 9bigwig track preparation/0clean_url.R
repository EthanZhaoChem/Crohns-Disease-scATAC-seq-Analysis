xx.1 <- 'http://194.113.195.165/data/gca20240417/anno1/'
folder.bigwig <- '/home/yuzhao1/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1/GroupBigWigs/anno1/'
celltypes <- list.files(folder.bigwig)
for (celltype in celltypes){
  cat(paste0(xx.1,celltype,'\n'))
}

xx.1 <- 'http://194.113.195.165/data/gca20240417/anno1_inf/'
folder.bigwig <- '/home/yuzhao1/yuzhao1/work/atac_gca2024/0dataset/5kmin_6TSS_DoubletRatio2_filtered1/GroupBigWigs/anno1_inf/'
celltypes <- list.files(folder.bigwig)
for (celltype in celltypes){
  cat(paste0(xx.1,celltype,'\n'))
}

