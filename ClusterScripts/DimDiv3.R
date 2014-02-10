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
load("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")

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

#################################################################################
#Null Models of Betadiversity - Step 3: Comparing Null phylogenetic to observed
#################################################################################

#How many rows?
paste(nrow(phylo_HL),"phylogenetic assemblage comparisons to evaluate")

#Test inner function
system.time(Null_P(1))

Null_dataframe[1:100,]
dim(Null_dataframe)

Null_dataframe<-Null_dataframe[!is.na(Null_dataframe$Phylosor.Phylo),]

head(Null_dataframe)
dim(Null_dataframe)

cluster <- getMPIcluster()

#Export objects to cluster
clusterExport(cluster, list("phylo_HL", "Null_dataframe"))
invisible(clusterEvalQ(cluster,source("/home1/02443/bw4sz/DimDiv/ClusterSourceFunctions.R")))

null_Phylo<-rbind.fill(clusterApply(cluster,1:nrow(phylo_HL),Null_P))

#Bind together the null model outputs
colnames(null_Phylo)<-c("To","From",paste(colnames(Null_dataframe)[!colnames(Null_dataframe) %in% c("To","From","Iteration","Sorenson")],"Null",sep="_"))
save.image("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")
