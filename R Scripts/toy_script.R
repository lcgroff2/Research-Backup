library(readxl)
library(plyr)
library(dplyr)
library(ggplot2)
library(readr)
library(xlsx)
library(stringr)
library(cluster)
library(vegan)
library(factoextra)
library(dendextend)
library(pheatmap)
library(gplots)

dataNumbers <- as.data.frame(read_csv("C:/Users/CLowe/Desktop/dataNumbers.csv", col_types = cols(NORMAL = col_double())))
row.names(dataNumbers) <- as.character(dataNumbers[,1])
dataNumbers[,1] <- NULL

dataSimple <- dataNumbers

for (i in 1:nrow(dataSimple)) {
  for (j in 1:ncol(dataSimple)) {
    if(dataSimple[i,j]  > 1) {
      dataSimple[i,j] <- 1
    }
    if(dataSimple[i,j] < 1) {
      dataSimple[i,j] <- 0
    }
  }
}

data2pairs <- dataNumbers

for (i in 1:nrow(dataNumbers)) {
  if (i+1 < nrow(dataNumbers)) {
    for (j in (i+1):nrow(dataNumbers)) {
      data2pairs[paste(colnames(data2pairs[i]),colnames(data2pairs[j]),sep = " ")] <- data2pairs[i]*data2pairs[j]
    }
  }  
}


pheatmap(as.matrix(dataSimple),cluster_cols = TRUE,cutree_rows = 4, cutree_cols = 6)

pheatmap(as.matrix(data2pairs),cluster_cols = TRUE, cluster_rows = TRUE,cutree_rows = 4, cutree_cols = 6,clustering_method = "ward.D",main = "ward.D")
pheatmap(as.matrix(data2pairs),cluster_cols = TRUE, cluster_rows = TRUE,cutree_rows = 4, cutree_cols = 6,clustering_method = "ward.D2",main = "ward.D2")
pheatmap(as.matrix(data2pairs),cluster_cols = TRUE, cluster_rows = TRUE,cutree_rows = 4, cutree_cols = 6,clustering_method = "single",main = "single")
pheatmap(as.matrix(data2pairs),cluster_cols = TRUE, cluster_rows = TRUE,cutree_rows = 4, cutree_cols = 6,clustering_method = "complete",main = "complete")
pheatmap(as.matrix(data2pairs),cluster_cols = TRUE, cluster_rows = TRUE,cutree_rows = 4, cutree_cols = 6,clustering_method = "average",main = "average")
pheatmap(as.matrix(data2pairs),cluster_cols = TRUE, cluster_rows = TRUE,cutree_rows = 4, cutree_cols = 6,clustering_method = "median",main = "median")
pheatmap(as.matrix(data2pairs),cluster_cols = TRUE, cluster_rows = TRUE,cutree_rows = 4, cutree_cols = 6,clustering_method = "centroid",main = "centroid")

test <- data2pairs[c("POISON DRAGON","STEEL DRAGON","GHOST DRAGON","BUG ELECTRIC",'FIRE DRAGON',"PSYCHIC FAIRY")]
pheatmap(as.matrix(test),cluster_cols = TRUE)