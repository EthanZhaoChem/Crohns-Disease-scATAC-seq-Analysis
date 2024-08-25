library(plyr)
library(dplyr)

out.dir <- '~/yuzhao1/work/atac_gca2024/0manu/plots/2UCSC_celltype_multiRegion/'
marker_regions_union <- c(
  'chr12:71,438,802-71,440,655 LGR5
chr13:53,027,655-53,029,876 OLFM4
chr6:52,762,926-52,764,137 GSTA2
chr11:116,822,858-116,823,961 APOA4
chr8:85,464,813-85,466,594 CA2
chr19:41,686,846-41,689,650 CEACAM7
chr11:1,073,859-1,075,702 MUC2
chr1:44,787,367-44,789,085 BEST4
chr11:120,238,222-120,241,878 POU2F3
chr14:92,922,258-92,924,073 CHGA
chr8:7,056,277-7,057,211 DEFA5
chr11:118,342,389-118,343,014 CD3D
chr12:6,789,190-6,789,800 CD4
chr2:86,785,040-86,787,691 CD8A
chr19:54,905,370-54,906,658 NCR1
chr12:9,826,129-9,828,685 KLRF1
chr13:67,047,722-67,049,322 PCDH9
chr12:122,023,957-122,025,239 BCL7A
chr16:28,931,504-28,932,366 CD19
chr2:20,206,456-20,208,010 SDC1
chr3:148,877,609-148,878,972 CPA3
chr1:22,636,123-22,636,840 C1QA
chr12:10,051,223-10,054,125 CLEC9A
chr1:153,390,533-153,391,574 S100A8
chr17:34,284,927-34,286,129 CCL11
chr1:163,202,125-163,203,834 RGS5
chr17:64,327,010-64,328,714 PECAM1
chr21:46,604,136-46,606,327 S100B'
)


temp_lineage <- 'union'

marker_regions <- paste0('marker_regions_', temp_lineage) %>% as.name(.) %>% eval(.)
marker_regions <- marker_regions %>% gsub('-',' ', .) %>% gsub(':',' ', .) %>% gsub(',','', .) 
marker_region_seqs <- marker_regions %>% strsplit(., '\n') %>% unlist() %>% strsplit(., ' ') %>% sapply(.,`[[`,1)
marker_region_lef <- marker_regions %>% strsplit(., '\n') %>% unlist() %>% strsplit(., ' ') %>% sapply(.,`[[`,2) %>% as.numeric()
marker_region_right <- marker_regions %>% strsplit(., '\n') %>% unlist() %>% strsplit(., ' ') %>% sapply(.,`[[`,3) %>% as.numeric()
marker_region_centers <- floor((marker_region_lef + marker_region_right)/2)

out_data <- cbind(marker_region_seqs,
                  marker_region_centers - 300,
                  marker_region_centers + 300)
write.table(out_data,file=paste0(out.dir, "markerRegions_", temp_lineage, ".bed"),quote=F,row.names=F,col.names=F,sep="\t")













