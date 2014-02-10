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

#Calculate Observed Betadiverisity
system.time(beta_metrics<-beta_all(comm=comm,tree=tree,traits=mon))

#Visualizations of the beta metrics
head(beta_metrics)

#Get rid of the nestedness components
beta_metrics<-beta_metrics[,colnames(beta_metrics) %in% c("To","From","MNTD","Phylosor.Phylo","Sorenson") ]

#####################################################
#Merge Betadiversity and Environmnetal Dissimilairity
#####################################################

data.merge<-merge(compare.env,beta_metrics,by=c("To","From"))

#save.image("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")

print("Observed Betadiversity Found")

#Begin Part 2

#load("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")

#################################################################################
#Null Models of Betadiversity - Step 1: Observed Deliniations of High and Low
#################################################################################

#Deliniate observed betadiversity based on the .3 and .7 quantiles 
highlowdelimF<-quantile(data.merge$MNTD,c(.3,.7))
highlowdelimP<-quantile(data.merge$Phylosor.Phylo,c(.3,.7))
highlowdelimS<-quantile(data.merge$Sorenson,c(.3))

data.merge$observedT<-cut(data.merge$Sorenson,breaks=c(0,highlowdelimS,.9999,1),include.lowest=TRUE,labels=c("Low","Random","High"))
data.merge$observedP<-cut(data.merge$Phylosor.Phylo,breaks=c(0,highlowdelimP,1),include.lowest=TRUE,labels=c("Low","Random","High"))

#For trait betadiversity get the .3 and .7 quantile
highlowdelim<-quantile(data.merge$MNTD,c(.3,.7))
data.merge$observedF<-cut(data.merge$MNTD,breaks=c(0,highlowdelim,max(data.merge$MNTD)),include.lowest=TRUE,labels=c("Low","Random","High"))

print("Deliniated high and low betadiversity")

#################################################################################
#Null Models of Betadiversity - Step 2: Taxonomic
#################################################################################

#Null model of taxonomic diversity with respect to richness
#Remove the random assemblages from the null model computation
sorenson_HL<-data.merge[!data.merge$observedT %in% "Random",]

paste(nrow(sorenson_HL),"taxonomic assemblage comparisons to evaluate")

#Export objects to cluster
clusterExport(cluster, list("sorenson_HL","comm"))

#source functions on each cluster
clusterEvalQ(cluster,source("/home1/02443/bw4sz/DimDiv/ClusterSourceFunctions.R"))

#Compute taxonomic betadiversity nulls
null_taxlists<-rbind.fill(clusterApply(cluster,1:nrow(sorenson_HL),TaxComp))

save.image("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")