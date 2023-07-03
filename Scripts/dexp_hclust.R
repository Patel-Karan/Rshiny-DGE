options(java.parameters = "-Xmx30g")
library("DESeq2")


d_exp <- function(dframe, nC, nT, hclust_method, dist_method){
  countData  <- dframe
  nControl <- nC
  nCase <- nT
  condition <- c(rep("Control",nControl),rep("Case",nCase))
  colData <- data.frame(row.names=colnames(countData), condition=factor(condition, levels=c('Control','Case')))
  dataset <- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~condition)
  dds <- DESeq(dataset)
  result <- results(dds)
  colnames(countData) <- paste(colnames(countData),"Raw.Read.Count", sep = "_")
  Normalized_count=(counts(dds, normalized=TRUE))
  colnames(Normalized_count) <- paste(colnames(Normalized_count),"Normalized.Read.Count", sep = "_")
  res <- data.frame(countData,Normalized_count,result, check.names = FALSE)
  res <- res[ order(res$log2FoldChange ), ]
  res <- res[, grep('Normalized.Read.Count', names(res))]
  names(res)=gsub(x=names(res),pattern="_Normalized.Read.Count",replacement="")
  res = res[rowSums(res[])>0,]
  
  
  # CLUSTERING
  res_t=t(res)
  clust_dist = dist(res_t,dist_method)
  clust_dist_df = melt(as.matrix(clust_dist), varnames = c("Sample1", "Sample2"),value.name = "Distance") %>% 
      filter(Distance != 0) %>% mutate(across(where(is.numeric), round, 2))
  set.seed(240)
  modelname<-hclust(clust_dist, method = hclust_method)
  
  # PCA
  pca <- prcomp(res_t)
  scores <- data.frame(pca$x[,1:2], check.names = FALSE)
  scores$NAMES <- rownames(scores)
  percent <- round((((pca$sdev)^2 / sum(pca$sdev^2))*100)[1:2])
  
  # OBJECT LIST'S
  obj <- list(model_df=list(m=modelname,df=res),distance=clust_dist_df,pca=list(pca_score=scores, p_percent = percent))
  
  return(obj)
}