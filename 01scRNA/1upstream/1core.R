library(sys)
library(dplyr)
library(stringr)
patients <- list.dirs('/home/yuzhao1/gca/GCA_scRNA/fastq/', recursive = F, full.names = F)
patients.paths <- paste0('/home/yuzhao1/gca/GCA_scRNA/fastq/', patients)
problematic.samples <- c()

# cellranger commands
for(i in 1:length(patients)){
  i.patient <- patients[[i]]
  i.patient.path <- paste0('/home/yuzhao1/gca/GCA_scRNA/fastq/', i.patient)
  i.patient.files <- list.files(path = i.patient.path, recursive = F)
  i.patient.samples <- c()
  
  # get unique sample names from a patient folder
  for (j.file in i.patient.files) {
    lane.start <- j.file %>% gregexpr('L00', .) %>% unlist()
    
    # identify samples that can't be grouped by lanes
    if(lane.start < 1){
      problematic.samples <- unique(c(problematic.samples, j.file))
      next
    }
    
    # identify the existed sample name
    sample.name <- substr(j.file, 1, lane.start-2)
    i.patient.samples <- unique(c(i.patient.samples, sample.name))
  }
  
  # print commands to submit job for each sample
  for (k.sample in i.patient.samples) {
    # clean up sampleID
    SampleID.start <- max(k.sample %>% gregexpr('HA', .) %>% unlist(),
                          k.sample %>% gregexpr('OR', .) %>% unlist())
    temp1 <- substr(k.sample, SampleID.start, nchar(k.sample))
    SampleID <- strsplit(temp1, split = '_S')[[1]][[1]]
    
    # remove S string in k.sample
    samplename.end <- as.numeric(str_locate(k.sample, SampleID))[2]
    k.sample.updated <- substr(k.sample, 1, samplename.end)
    cat('sbatch -J ', SampleID, ' cellranger.batch ',
        SampleID, " /home/yuzhao1/gca/GCA_scRNA/fastq/", i.patient, ' ', k.sample.updated,
        '\n',sep = '')
  }
}

# print samples that can't be grouped by lanes
for (i.problem in problematic.samples) {
  cat(i.problem, '\n')
}


# write SampleIDs to metedata folder
SampleIDs <- list.dirs('~/yuzhao1/work/final_GCArna/upstream', recursive = F, full.names = F)
write.table(SampleIDs, '~/yuzhao1/work/final_GCArna/metadata/SampleIDs.csv',
            col.names = 'SampleID', row.names = F)
SampleIDs <- read.csv("yuzhao1/work/final_GCArna/metadata/SampleIDs.csv")


# check cellranger out slurm file, successful or not
SampleIDs <- read.csv("yuzhao1/work/final_GCArna/metadata/SampleIDs.csv")
SampleIDs <- SampleIDs$SampleID
out.files <- paste0('~/yuzhao1/work/final_GCArna/upstream/', SampleIDs, '.out')
successes <- 0
for (i.outfile in out.files) {
  temp <- readLines(i.outfile)
  flag <- grepl('Pipestance completed successfully!', temp) %>% sum()
  if(flag >0){
    successes <- successes+1
  }
  if(flag < 1){
    cat(i.outfile, '\n')
  }
}

# check the number of successfull samples is same with #samples: all successes, mv output logs to log folder
library(filesstrings)
cat(successes)
if(successes == length(SampleIDs)){
  SampleIDs <- read.csv("yuzhao1/work/final_GCArna/metadata/SampleIDs.csv")
  SampleIDs <- SampleIDs$SampleID
  out.files <- paste0('~/yuzhao1/work/final_GCArna/upstream/', SampleIDs, '.out')
  err.files <- paste0('~/yuzhao1/work/final_GCArna/upstream/', SampleIDs, '.err')
  for (log in c(out.files, err.files)){
    file.move(log, '~/yuzhao1/work/final_GCArna/upstream/log_cellranger')
  }
}
length(list.files('~/yuzhao1/work/final_GCArna/upstream/log_cellranger', recursive = F))


# velocyto commands
SampleIDs <- read.csv("yuzhao1/work/final_GCArna/metadata/SampleIDs.csv")
SampleIDs <- SampleIDs$SampleID
for (SampleID in SampleIDs){
  sample.dir <- paste0('~/yuzhao1/work/final_GCArna/upstream/', SampleID)
  
  cat('sbatch -J ', SampleID,
  ' velocyto.batch ',sample.dir,'\n')
}

# check velocyto out slurm file, successful or not
SampleIDs <- read.csv("~/yuzhao1/work/final_GCArna/metadata/SampleIDs.csv")
SampleIDs <- SampleIDs$SampleID
out.files <- paste0('~/yuzhao1/work/final_GCArna/upstream/', SampleIDs, '.out')
successes <- 0
for (i.outfile in out.files) {
  temp <- readLines(i.outfile)
  flag <- grepl('Terminated Succesfully!', temp) %>% sum()
  if(flag >0){
    successes <- successes+1
  }
  if(flag < 1){
    cat(i.outfile, '\n')
  }
}

# check the number of successfull samples is same with #samples: all successes, mv output logs to log folder
library(filesstrings)
cat(successes)
if(successes == length(SampleIDs)){
  SampleIDs <- read.csv("yuzhao1/work/final_GCArna/metadata/SampleIDs.csv")
  SampleIDs <- SampleIDs$SampleID
  out.files <- paste0('~/yuzhao1/work/final_GCArna/upstream/', SampleIDs, '.out')
  err.files <- paste0('~/yuzhao1/work/final_GCArna/upstream/', SampleIDs, '.err')
  for (log in c(out.files, err.files)){
    file.move(log, '~/yuzhao1/work/final_GCArna/upstream/log_velo/')
  }
}
length(list.files('~/yuzhao1/work/final_GCArna/upstream/log_velo/', recursive = F))















