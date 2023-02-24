library("clusterProfiler")
library("topGO")
library("Rgraphviz")
library(org.EcK12.eg.db) 

data<-read.csv('CIP-degs.csv')

ego <- enrichGO(
  gene = data$CIP_2h_n,
  OrgDb = org.EcK12.eg.db,
  keyType = 'SYMBOL',
  ont = "ALL", 
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05
)
dotplot(ego,font.size=12)+theme_classic()


