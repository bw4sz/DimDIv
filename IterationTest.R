#load data from cluster
load("C:/Users/Jorge/Dropbox/Shared Ben and Catherine/DimDivRevision/ClusterResults/DimDivRevisionCluster.RData")

head(Null_dataframe)

null_rows<-Null_dataframe[Null_dataframe$To==7 & Null_dataframe$From==13,]


hist(null_rows$Phylosor.Phylo)

hist(sample(null_rows$Phylosor.Phylo,100),add=TRUE,col="red")


qI<-lapply(seq(100,500,50),function(x){
  replicate(1000,quantile(sample(null_rows$Phylosor.Phylo,x),c(.95,.05)))
  })
names(qI)<-seq(100,500,50)

qI<-melt(qI)

ggplot(qI,aes(L1,value,col=X1)) + geom_point() + xlab("Iterations") + ylab("Phylosor.Phylo") + labs(col="Quantile")
ggsave("C:/Users/Jorge/Dropbox/Shared Ben and Catherine/DimDivRevision/Results/Iterationtest.svg")
