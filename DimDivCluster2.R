#Taxonomic, Phylogenetic and Functional Diversity in South American Hummingbirds
#Manuscript submitted to the American Naturalist

#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW
#Boris A. Tinoco, Juan L. Parra, PhD, Leone M. Brown, PhD, Gary Stiles, PhD, Jim A. McGuire, PhD, Catherine H. Graham, PhD

library(Rmpi)
require(vegan)
require(reshape)
require(ape)
require(picante)
require(ggplot2)
library(doSNOW)
library(foreach)
require(GGally)
require(stringr)
require(scales)
library(foreach)
require(raster)
require(boot)


cluster <- getMPIcluster()
 
# Print the hostname for each cluster member
 sayhello <- function()
 {
 	info <- Sys.info()[c("nodename", "machine")]
 	paste("Hello from", info[1], "with CPU type", info[2])
 }
 
names <- clusterCall(cluster, sayhello)
print(unlist(names))


registerDoSNOW(cluster)

###########################
###############Read in data
###########################

###Define Source Functions

#list loaded packages
(.packages())

load("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")
###
#Begin Part 2

#################################################################################
table(null_taxlists$Sorenson_Null)

#Richness of the assemblages
richness_sites<-apply(siteXspp,1,sum)

#Composition of each site
sp.list<-apply(siteXspp,1,function(x){
  names(x[which(x==1)])
})

nrow(null_taxlists)

head(null_taxlists)
tail(null_taxlists)
head(splist)


####################

splist<-colnames(siteXspp)
x=2

  #Richness of the two assemblages
  richness_To<-richness_sites[names(richness_sites) %in% null_taxlists[x,]$To]
  richness_From<-richness_sites[names(richness_sites) %in% null_taxlists[x,]$From]
  
  #Create slots for each community
  nullTo<-matrix(nrow=richness_To,ncol=1)
  nullFrom<-matrix(nrow=richness_From,ncol=1)
  
  #number of species overlapping
  species.To<-sp.list[names(sp.list) %in% null_taxlists[x,]$To][[1]]
  species.From<-sp.list[names(sp.list) %in% null_taxlists[x,]$From][[1]]
  overlapping<-sum(species.To %in% species.From)
  
overlap.species<-NULL
  #If there are species overlaps
  if (overlapping > 0){
    #Fill slots in both assemblages with the same randomly drawn species
    overlap.species<-sample(splist,overlapping,replace=FALSE)
    nullTo[1:overlapping,]<-overlap.species
    nullFrom[1:overlapping,]<-overlap.species
  }
  
  #Fill the rest of the assemblage with randomly drawn species
  #If communities are entirely nested, there is no more species draw.
  #Take out the species in overlap species, cant be picked twice
  if(!overlapping==richness_To){
    nullTo[(overlapping+1):richness_To,]<-sample(splist[!splist %in% overlap.species],richness_To-overlapping,replace=FALSE)}
  

speciestopick<-splist[!splist %in% overlap.species & !splist %in% nullTo]

  if(!overlapping==richness_From){
    nullFrom[(overlapping+1):richness_From,]<-sample(speciestopick,richness_From-overlapping,replace=FALSE)}
  
  #Create siteXspp table for new assemblages
  null_a<-melt(list(nullTo,nullFrom))
  
  #create binary matrix
  null_siteXspp<-(!is.na(cast(null_a,L1~value)[,-1]))*1
  rownames(null_siteXspp)<-c(as.character(null_taxlists[x,]$To),as.character(null_taxlists[x,]$From))
  
  #Compute phylogenetic and trait metrics on null assemblages
 print(beta_all(comm=null_siteXspp,tree=tree,traits=mon))

############################
#Null model for taxonomic phylogenetic and trait
############################

Null_dataframe<-foreach(p=1:500,.combine=rbind) %dopar%{
sink("/home1/02443/bw4sz/DimDiv/DimSink.txt")
require(reshape)
require(picante)
require(foreach)

print(paste("Iteration",p))
print(Sys.time())

nullPTframe<-foreach(x=1:nrow(null_taxlists),.combine=rbind) %do% {

splist<-colnames(siteXspp)

  #Richness of the two assemblages
  richness_To<-richness_sites[names(richness_sites) %in% null_taxlists[x,]$To]
  richness_From<-richness_sites[names(richness_sites) %in% null_taxlists[x,]$From]
  
  #Create slots for each community
  nullTo<-matrix(nrow=richness_To,ncol=1)
  nullFrom<-matrix(nrow=richness_From,ncol=1)
  
  #number of species overlapping
  species.To<-sp.list[names(sp.list) %in% null_taxlists[x,]$To][[1]]
  species.From<-sp.list[names(sp.list) %in% null_taxlists[x,]$From][[1]]
  overlapping<-sum(species.To %in% species.From)
  
overlap.species<-NULL
  #If there are species overlaps
  if (overlapping > 0){
    #Fill slots in both assemblages with the same randomly drawn species
    overlap.species<-sample(splist,overlapping,replace=FALSE)
    nullTo[1:overlapping,]<-overlap.species
    nullFrom[1:overlapping,]<-overlap.species
  }
  
 
  #Fill the rest of the assemblage with randomly drawn species
  #If communities are entirely nested, there is no more species draw.
  #Take out the species in overlap species, cant be picked twice
  if(!overlapping==richness_To){
    nullTo[(overlapping+1):richness_To,]<-sample(splist[!splist %in% overlap.species],richness_To-overlapping,replace=FALSE)}
  

speciestopick<-splist[!splist %in% overlap.species & !splist %in% nullTo]

  if(!overlapping==richness_From){
    nullFrom[(overlapping+1):richness_From,]<-sample(speciestopick,richness_From-overlapping,replace=FALSE)}
  
  #Create siteXspp table for new assemblages
  null_a<-melt(list(nullTo,nullFrom))
  
  #create binary matrix
  null_siteXspp<-(!is.na(cast(null_a,L1~value)[,-1]))*1
  rownames(null_siteXspp)<-c(as.character(null_taxlists[x,]$To),as.character(null_taxlists[x,]$From))
  
  #Compute phylogenetic and trait metrics on null assemblages
  null_beta_metrics<-beta_all(comm=null_siteXspp,tree=tree,traits=mon)
}

print(head(nullPTframe))
sink()
out<-data.frame(nullPTframe,Iteration=p)
return(out)
}

#Keep desired columns, ignoring nestedness components
Null_dataframe<-Null_dataframe[,colnames(Null_dataframe) %in% c("To","From","MNTD","Phylosor.Phylo","Sorenson","Iteration") ]

save.image("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")
