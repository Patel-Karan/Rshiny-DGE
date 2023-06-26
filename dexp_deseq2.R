options(java.parameters = "-Xmx30g")
#library(xlsx)
library("DESeq2")


d_exp <- function(dframe, nC, nT){
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
  res <- data.frame(countData,Normalized_count,result)
  Final <- res[ order(res$log2FoldChange ), ]
  return(Final)
}