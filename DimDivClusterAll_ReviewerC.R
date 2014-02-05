#Taxonomic, Phylogenetic and Functional Diversity in South American Hummingbirds
#Manuscript submitted to the American Naturalist

#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW
#Boris A. Tinoco, Juan L. Parra, PhD, Leone M. Brown, PhD, Gary Stiles, PhD, Jim A. McGuire, PhD, Catherine H. Graham, PhD

#require(Rmpi)
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

#Cluster Path
#load("/home1/02443/bw4sz/DimDiv/DimDivRevision.RData")

#Locally
load(paste(droppath,"Shared Ben and Catherine/DimDivRevision/Results/DimDivRevision.RData",sep=""))

#Create Cluster

cluster<-makeCluster(2,"SOCK")
#cluster <- getMPIcluster()

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
comm<-siteXspp[1:30,]

dim(comm)

ls()

#Source Functions Cluster
source("/home1/02443/bw4sz/DimDiv/ClusterSourceFunctions.R")

#Source Functions Locally
source("ClusterSourceFunctions.R")

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

#Deliniate observed betadiversity
data.merge$observedT<-cut(data.merge$Sorenson,breaks=c(0,.3,.7,1),include.lowest=TRUE,labels=c("Low","Random","High"))
data.merge$observedP<-cut(data.merge$Phylosor.Phylo,breaks=c(0,.3,.7,1),include.lowest=TRUE,labels=c("Low","Random","High"))

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

#Source functions locally
invisible(clusterEvalQ(cluster,source("ClusterSourceFunctions.R")))

#Compute taxonomic betadiversity nulls

system.time(null_taxlists<-rbind.fill(clusterApply(cluster,1:nrow(sorenson_HL),TaxComp)))

save.image("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")

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

j <- getMPIcluster()

clusterExport(j, list("phylo_HL", "splist","sp.list","siteXspp","richness_sites","beta_all","tree","mon"))
clusterEvalQ(j,source("/home1/02443/bw4sz/DimDiv/ClusterSourceFunctions.R"))

Null_dataframe0<-clusterApply(j,1:100,Phylo_N)

#If this has failed, read and aggregate from file 
head(Null_dataframe0[[1]])

Null_dataframe<-rbind.fill(Null_dataframe0)

#Keep desired columns, ignoring nestedness components
Null_dataframe<-Null_dataframe[,colnames(Null_dataframe) %in% c("To","From","MNTD","Phylosor.Phylo","Sorenson","Iteration") ]

#fil<-list.files("/home1/02443/bw4sz/DimDiv/Iterations",full.names=TRUE)
#fil.l<-lapply(fil,read.csv,row.names=1)
#Null_dataframe<-rbind.fill(fil.l)
save.image("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")

#################################################################################
#Null Models of Betadiversity - Step 3: Comparing Null phylogenetic to observed
#################################################################################

#How many rows?
paste(nrow(phylo_HL),"phylogenetic assemblage comparisons to evaluate")

#Test inner function
Null_P(1)

cluster <- getMPIcluster()

#Export objects to cluster
clusterExport(cluster, list("phylo_HL", "Null_dataframe"))
clusterEvalQ(cluster,source("/home1/02443/bw4sz/DimDiv/ClusterSourceFunctions.R"))

null_Phylo<-rbind.fill(clusterApply(cluster,1:nrow(phylo_HL),Null_P))

#Bind together the null model outputs
colnames(null_Phylo)<-c("To","From",paste(colnames(Null_dataframe)[!colnames(Null_dataframe) %in% c("To","From","Iteration","Sorenson")],"Null",sep="_"))
save.image("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")

#################################################################################
#Null Models of Betadiversity - Step 4: Null Distribution of Trait betadiversity
#################################################################################

trait_HL<-data.merge[!data.merge$observedF %in% "Random",]
clusterExport(cluster, list("trait_HL", "splist","sp.list","siteXspp","richness_sites","beta_all","tree","mon"))
clusterEvalQ(cluster,source("/home1/02443/bw4sz/DimDiv/ClusterSourceFunctions.R"))

Null_dataframeTrait<-clusterApply(cluster,1:2,Trait_N)

#If this has failed, read and aggregate from file 
head(Null_dataframeTrait[[1]])

Null_dataframeTraitF<-rbind.fill(Null_dataframeTrait)

#Keep desired columns, ignoring nestedness components
Null_dataframeTraitF<-Null_dataframeTraitF[,colnames(Null_dataframeTraitF) %in% c("To","From","MNTD","Phylosor.Phylo","Sorenson","Iteration") ]

#fil<-list.files("/home1/02443/bw4sz/DimDiv/Iterations",full.names=TRUE)

#fil.l<-lapply(fil,read.csv,row.names=1)

#Null_dataframeTraitF<-rbind.fill(fil.l)

save.image("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")

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
data.df.null<-merge(data.merge,null_taxlists,by=c("To","From"))
data.df.null<-merge(data.df.null,null_Phylo,by=c("To","From"))
data.df.null<-merge(data.df.null,null_Trait,by=c("To","From"))

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
