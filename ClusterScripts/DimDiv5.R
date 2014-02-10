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
#Null Models of Betadiversity - Step 5: Comparing Null trait to observed
#################################################################################

#How many rows?
paste(nrow(Null_dataframeTraitF),"assemblage comparisons to evaluate")

#Test inner function
Null_T(1)

#Export objects to cluster
clusterExport(cluster, list("trait_HL", "Null_dataframeTraitF"))
clusterEvalQ(cluster,source("/home1/02443/bw4sz/DimDiv/ClusterSourceFunctions.R"))

null_Trait<-rbind.fill(clusterApply(cluster,1:nrow(trait_HL),Null_T))

#Bind together the null model outputs
colnames(null_Trait)<-c("To","From",paste(colnames(Null_dataframeTraitF)[!colnames(Null_dataframeTraitF) %in% c("To","From","Iteration","Sorenson")],"Null",sep="_"))

print(nrow(null_Trait))

save.image("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")

#################################################################################
#Null Models of Betadiversity - Step 6: Combine outputs
#################################################################################

#Combine the environmental, observed and null metrics into a huge dataframe
data.df.null<-merge(data.merge,null_taxlists,by=c("To","From"),all=TRUE)
data.df.null<-merge(data.df.null,null_Phylo,by=c("To","From"),all=TRUE)
data.df.null<-merge(data.df.null,null_Trait,by=c("To","From"),all=TRUE)

dim(data.df.null)

#Legacy correction, data.merge is data.df, sorry
data.df<-data.merge

setwd("/home1/02443/bw4sz/")
#Write to file
write.csv(data.df,"/home1/02443/bw4sz/DimDiv/FinalData.csv")
write.csv(data.df.null,"/home1/02443/bw4sz/DimDiv/FinalDataNull.csv")

#Or save data
save.image("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")

#Data Generation Complete
##########################################################################################
##########################################################################################

stopCluster(cluster)
#Quit
q()
#Don't save
n
#exit terminal
exit
