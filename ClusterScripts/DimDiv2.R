#Taxonomic, Phylogenetic and Functional Diversity in South American Hummingbirds
#Manuscript submitted to the American Naturalist

#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW
#Boris A. Tinoco, Juan L. Parra, PhD, Leone M. Brown, PhD, Gary Stiles, PhD, Jim A. McGuire, PhD, Catherine H. Graham, PhD

require(Rmpi)
require(doSNOW)
require(vegan)
require(reshape)
require(ape)
require(picante)
require(ggplot2)
library(snow)
library(foreach)
require(stringr)
require(scales)
library(foreach)
require(raster)
require(boot)


##########################################################
#Read in data
##########################################################

#Load in the data
#droppath<-"C:/Users/Jorge/Dropbox/"

#Cluster
load("/home1/02443/bw4sz/DimDiv/DimDivRevision.RData")

#Locally
#load(paste(droppath,"Shared Ben and Catherine/DimDivRevision/Results/DimDivRevision.RData",sep=""))

#Create Cluster

#cluster<-makeCluster(4,"SOCK")
cluster <- getMPIcluster()

# Print the hostname for each cluster member
sayhello <- function()
{
  info <- Sys.info()[c("nodename", "machine")]
  paste("Hello from", info[1], "with CPU type", info[2])
}

names <- clusterCall(cluster, sayhello)
print(unlist(names))

#registerDoSNOW(cluster)

#list loaded packages
(.packages())

#Subset assemblage matrix for testing. 
comm<-siteXspp[,]

dim(comm)

ls()

#Source Functions
source("/home1/02443/bw4sz/DimDiv/ClusterSourceFunctions.R")

load("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")


#################################################################################
#Null Models of Betadiversity - Step 2: Null distribution of Phylogentic
#################################################################################

phylo_HL<-data.merge[!data.merge$observedP %in% "Random",]

#Richness of the assemblages
richness_sites<-apply(siteXspp,1,sum)

#Composition of each site
sp.list<-apply(siteXspp,1,function(x){
  names(x[which(x==1)])
})

splist<-colnames(siteXspp)

clusterExport(cluster, list("phylo_HL", "splist","sp.list","siteXspp","richness_sites","beta_all","tree","mon"))
clusterEvalQ(cluster,source("/home1/02443/bw4sz/DimDiv/ClusterSourceFunctions.R"))

Null_dataframe0<-clusterApply(cluster,1:500,Phylo_N)

#If this has failed, read and aggregate from file 
head(Null_dataframe0[[1]])

Null_dataframe<-rbind.fill(Null_dataframe0)

#Keep desired columns, ignoring nestedness components
Null_dataframe<-Null_dataframe[,colnames(Null_dataframe) %in% c("To","From","MNTD","Phylosor.Phylo","Sorenson","Iteration") ]

#fil<-list.files("/home1/02443/bw4sz/DimDiv/Iterations",full.names=TRUE)
#fil.l<-lapply(fil,read.csv,row.names=1)
#Null_dataframe<-rbind.fill(fil.l)

save.image("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")
