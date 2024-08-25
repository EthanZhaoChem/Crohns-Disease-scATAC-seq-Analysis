library(biomaRt)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)

##-- 1.1 convert differential expression dataframe into a named and ordered logFC vector and perform GSEA
# # example usage
# ddf <- data.frame(gene=rownames(df_gene_scores),
#                   z_logFC=gvec)
# gvec <- GetGenelist(ddf, by = 'z_logFC')
# reslist <- myEnrichGSEA_msigdb(gvec)
# resDf <- do.call("rbind", lapply(reslist, function(sobj) sobj@result))

GetGenelist <- function(ddf, by=NULL){
  if (is.data.frame(ddf)){
    if (is.null(by)){
      xx <- colnames(ddf)[grep("log",colnames(ddf))] 
    } else {
      xx <- by
    }
    gvec <- ddf[[xx]]
    if (!("gene" %in% colnames(ddf))) ddf$gene <- row.names(ddf)
  } else {
    gvec <- ddf
    ddf <- data.frame(gene=ddf, stringsAsFactors = F)
  }
  gdf <- biomaRt::select(org.Hs.eg.db, ddf$gene, "ENTREZID", "SYMBOL")
  gdf <- gdf[!base::duplicated(gdf$SYMBOL),]
  names(gvec) <- as.character(gdf$ENTREZID)
  gvec <- sort(gvec, decreasing = T)
  gvec <- gvec[!is.na(names(gvec))]
  return(gvec)
}


myEnrichGSEA_msigdb <- function(gvec, pthresh=0.05){
  require(msigdbr)
  reslist <- list()
  # catlbl <- c("H","C2","C5", "C7", "C8")
  catlbl <- c("H","C2","C5")
  for (cat in catlbl){
    m_df <- msigdbr(species = "Homo sapiens", category = c(cat)) %>% 
      dplyr::select(gs_name, entrez_gene)
    pathway <- GSEA(gvec, TERM2GENE = m_df, pvalueCutoff = pthresh, pAdjustMethod = "fdr", minGSSize = 10, maxGSSize = 5000)
    pathway <- setReadable(pathway, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
    reslist[[cat]] <- pathway
  }
  return(reslist)
}

##-- 1.2 perform OverRepresentation
# # genes are gene symbols
# gdf <- biomaRt::select(org.Hs.eg.db, genes, "ENTREZID", "SYMBOL")
# gdf <- gdf[!base::duplicated(gdf$SYMBOL),]
# genes_id <- as.character(gdf$ENTREZID)
# resDf <- do.call("rbind", lapply(reslist, function(sobj) sobj@result))

myOverRepresentation_msigdb <- function(genes_id, pvalueCutoff=0.05){
  reslist <- list()
  # catlbl <- c("H","C2","C5", "C7", "C8")
  catlbl <- c("H","C2","C5")
  for (cat in catlbl){
    m_t2g <- msigdbr(species = "Homo sapiens", category = c(cat)) %>% 
      dplyr::select(gs_name, entrez_gene)
    em <- enricher(genes_id, TERM2GENE=m_t2g, pAdjustMethod = "fdr", minGSSize = 10, maxGSSize = 5000, pvalueCutoff = pvalueCutoff)
    pathway <- setReadable(em, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
    reslist[[cat]] <- pathway
  }
  return(reslist)
}
