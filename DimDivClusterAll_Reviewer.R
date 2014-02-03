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


#Source Functions
source(ClusterSourceFunctions.R)

#Create Cluster

#cluster<-makeCluster(4,"SOCK")
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

##########################################################
#Read in data
##########################################################

#Load in the data
#load("/home1/02443/bw4sz/DimDiv/DimDivRevision.RData")
load(paste(droppath,"Shared Ben and Catherine/DimDivRevision/500Iterations/DimDivRevisionCluster.RData",sep=""))

comm<-siteXspp[1:50,]

dim(comm)

ls()

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
###
#End Part 1
###

#Begin Part 2

#load("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")

#################################################################################
#Null Models of Betadiversity - Step 1: Observed Deliniations of High and Low
#################################################################################

#Deliniate observed betadiversity
data.merge$observedT<-cut(data.merge$Sorenson,breaks=c(0,.3,.7,1),include.lowest=TRUE,labels=c("Low","Random","High"))
data.merge$observedP<-cut(data.merge$Phylosor.Phylo,breaks=c(0,.3,.7,1),include.lowest=TRUE,labels=c("Low","Random","High"))
data.merge$observedF<-cut(data.merge$MNTD,breaks=c(0,.3,.7,1),include.lowest=TRUE,labels=c("Low","Random","High"))


#Null model of taxonomic diversity with respect to richness
#Remove the random assemblages from the null model computation
sorenson_HL<-data.merge[!data.merge$observedT %in% "Random",]

#Export objects to cluster
#clusterExport(cluster, list("sorenson", "tax_nullM"))

#null_taxlists<-rbind.fill(clusterApply(cluster,1:nrow(sorenson),Sorenson_N))

#For each pair of assemblages create 1000 randomizations of assemblages of the same size
Sorenson_N<-lapply(1:nrow(sorenson_HL),function(x){
  
  #Get row in matrix
  r<-sorenson_HL[x,]
  
  null_rows<-t(replicate(10,TaxN(r),simplify="matrix"))
  colnames(null_rows)<-c("To","From","Sorenson")
  
  #Create a distribution of null values
  test_stat<-ecdf(null_rows[,"Sorenson"]) (r[,"Sorenson"])
  
  if(test_stat >= .95) {answer<-"High"}
  if(test_stat <= .05) {answer<-"Low"}
  if( test_stat <= .95 & test_stat >= .05 ) {answer<-"Random"}
  
  return(data.frame(To=r[,"To"],From=r[,"From"],Sorenson_Null=answer))
})

#Create a dataframe
Sorensond_df<-rbind.fill(Sorenson_N)

#merge with values
data.merge<-merge(data.merge,Sorensond_df,all=TRUE)

####Set up information for phylogenetic and trait betadiversity

#Richness of the assemblages
richness_sites<-apply(siteXspp,1,sum)

#Composition of each site
sp.list<-apply(siteXspp,1,function(x){
  names(x[which(x==1)])
})

#Get list of asssemblages to perform function on, only want non random tax, no need to waste the calculations
data.merge<-merge(data.merge,Sorensond_df,all=TRUE)

#Only compute null model on assemblage combinations that have observed high or low betadiversity, and null model for tax is still high or low
null_taxlists<-data.merge[data.merge$Sorenson_Null %in% c("High","Low") & data.merge$observedP %in% c("High","Low") & data.merge$observedF %in% c("High","Low"),]

head(null_taxlists)

##############################################################################
#Calculate Phylogenetic and Trait Diversity with respect to Taxonomic Diversity
##############################################################################

splist<-colnames(siteXspp)

############################
#Null model for taxonomic phylogenetic and trait
############################

#clusterExport(cluster, list("null_taxlists", "splist","sp.list","siteXspp","richness_sites","beta_all","tree","mon"))
#Null_dataframe0<-clusterApply(cluster,1:500,Phylo_TraitN)

#If this has failed, read and aggregate from file 
#Null_dataframe0[[1]]
#Null_dataframe<-rbind.fill(Null_dataframe0)

#Keep desired columns, ignoring nestedness components
#Null_dataframe<-Null_dataframe[,colnames(Null_dataframe) %in% c("To","From","MNTD","Phylosor.Phylo","Sorenson","Iteration") ]

#fil<-list.files("/home1/02443/bw4sz/DimDiv/Iterations",full.names=TRUE)

#fil.l<-lapply(fil,read.csv,row.names=1)

#Null_dataframe<-rbind.fill(fil.l)

save.image("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")

#Begin Part 3
nrow(Null_dataframe)
######################################################################
#######################Null Model Analysis, defining "High and Low"
######################################################################

##################################################################################
#Define Function
#For each pair of assemblage compare the observed metrics to the null distribution. 
###################################################################################

print(nrow(data.merge))

#Test inner function
Null_PT(1)

#Export objects to cluster
clusterExport(cluster, list("data.merge", "Null_dataframe","Null_PT"))

null_PhyloTrait<-rbind.fill(clusterApply(cluster,1:nrow(data.merge),Null_PT))

#Bind together the null model outputs
colnames(null_PhyloTrait)<-c("To","From",paste(colnames(Null_dataframe)[!colnames(Null_dataframe) %in% c("To","From","Iteration","Sorenson")],"Null",sep="_"))

#Combine the environmental, observed and null metrics into a huge dataframe
data.df.null<-merge(data.merge,null_PhyloTrait,by=c("To","From"))
data.df.null<-merge(data.df.null,null_taxlists,by=c("To","From"))

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
