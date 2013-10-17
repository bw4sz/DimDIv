#Taxonomic, Phylogenetic and Functional Diversity in South American Hummingbirds
#Manuscript submitted to the American Naturalist

#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW
#Boris A. Tinoco, Juan L. Parra, PhD, Leone M. Brown, PhD, Gary Stiles, PhD, Jim A. McGuire, PhD, Catherine H. Graham, PhD

#require(RMPI)
#require(doSNOW)
require(vegan)
require(reshape)
require(ape)
require(picante)
require(ggplot2)
library(snow)
library(foreach)
require(GGally)
require(stringr)
require(scales)
library(foreach)
require(raster)
require(boot)

cluster<-makeCluster(4,"SOCK")
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

###########################
###############Read in data
###########################

###Define Source Functions

#list loaded packages
(.packages())

##########################################################
#Read in data
##########################################################

#Load in the data
#load("/home1/02443/bw4sz/DimDiv/DimDivRevision.RData")

comm<-siteXspp[1:20,]

dim(comm)

ls()


######################################################
#Create a function for computing betadiversity metrics
#######################################################

#Find taxonomic betadiversity

#####################################
##Taxonomic Betadiversity
d<-as.matrix(vegdist(comm,binary=TRUE,upper=FALSE,diag=FALSE))
sorenson<-melt(d)[melt(upper.tri(d))$value,]
colnames(sorenson)<-c("To","From","Sorenson")
######################################

#Null model off taxonomic diversity with respect to richness
tax_nullM<-rbind.fill(lapply(1:5,function(x){
  null_comm<-commsimulator(comm,method="r0")
  d<-as.matrix(vegdist(null_comm,binary=TRUE,upper=FALSE,diag=FALSE))
  sorenson<-data.frame(melt(d)[melt(upper.tri(d))$value,],Iteration=x)
  colnames(sorenson)<-c("To","From","Sorenson","Iteration")
  return(sorenson)}))

#For each pair of assemblage compare the observed metrics to the null distribution. 

nrow(tax_nullM)
nrow(sorenson)
head(sorenson)

#Define function that runs through sorenson index.
Sorenson_N<-function(x){ 
  #Select Row
  rowS<-sorenson[x,]
  
  #Grab all the iteration rows that match these richness 
  null_rows<-tax_nullM[tax_nullM$To==rowS$To & tax_nullM$From==rowS$From,]
  
  #Create a distribution of null values
  test_stat<-ecdf(null_rows[,"Sorenson"]) (rowS[,"Sorenson"])
  
  if(test_stat >= .95) {answer<-"High"}
  if(test_stat <= .05) {answer<-"Low"}
  if( test_stat <= .95 & test_stat >= .05 ) {answer<-"Random"}
  
  return(data.frame(t(c(To=rowS$To,From=rowS$From,Sorenson_Null=answer))))
}

#Export objects to cluster
clusterExport(cluster, list("sorenson", "tax_nullM"))

null_taxlists<-rbind.fill(clusterApply(cluster,1:10,Sorenson_N))

stopCluster(cluster)

#####################################################
#Perform Phylogenetic and Trait Betadiversity
#####################################################

cluster<-makeCluster(4,"SOCK")
#cluster <- getMPIcluster()

beta_all<-function(comm,tree,traits){
  
  #####################################
  ##Taxonomic Betadiversity
  d<-as.matrix(vegdist(comm,binary=TRUE,upper=FALSE,diag=FALSE))
  sorenson<-melt(d)[melt(upper.tri(d))$value,]
  colnames(sorenson)<-c("To","From","Sorenson")
  ######################################
  
  #Phylosor Calculation see Bryant 2008
  phylo.matrix<-as.matrix(phylosor(comm,tree))
  diag(phylo.matrix)<-NA
  Phylosor.phylo<-melt(phylo.matrix)
  colnames(Phylosor.phylo)<-c("To","From","Phylosor.Phylo")
  Phylosor.phylo$Phylosor.Phylo<-1-Phylosor.phylo$Phylosor.Phylo
  
  #############################################
  #Non-dendrogram approach, functional approach employed by villeger 2013
  #############################################
  
  #Trait frame needs to match siteXSpp table
  mon_cut<-traits[rownames(traits) %in% colnames(comm),]
  
  #There are eight species without traits, take them out for just this portion of the analysis, keep the assemblage lsit
  siteXspp_traits<-comm[,colnames(comm) %in% rownames(mon_cut)]
  
  
  prc_traits<-prcomp(mon_cut)
  newSGdist <- dist(prc_traits$x)
  #source("/home1/02443/bw4sz/DimDiv/BenHolttraitDiversity.R")
  source("BenHolttraitDiversity.R")
  
  #create sp.list
  sp.list<-lapply(rownames(siteXspp_traits),function(k){
    x<-siteXspp_traits[k,]
    names(x[which(x==1)])
  })
  
  names(sp.list)<-rownames(siteXspp_traits)
  
  dists <- as.matrix(newSGdist)
  
  rownames(dists) <- rownames(mon_cut)
  colnames(dists) <- rownames(mon_cut)
  
  sgtraitMNTD <- sapply(rownames(siteXspp_traits),function(i){
    
    #Iterator count
    #print(round(which(rownames(siteXspp_traits)==i)/nrow(siteXspp_traits),3))
    
    #set iterator
    A<-i
    
    #
    out<-lapply(rownames(siteXspp_traits)[1:(which(rownames(siteXspp_traits) == i))], function(B) {MNND(A,B,sp.list=sp.list,dists=dists)})
    names(out)<-rownames(siteXspp_traits)[1:(which(rownames(siteXspp_traits) == i))]
    return(out)
  })
  
  names(sgtraitMNTD) <- rownames(siteXspp_traits)
  melt.MNTD<-melt(sgtraitMNTD)
  
  colnames(melt.MNTD)<-c("MNTD","To","From")
  
  #Combine with other metrics
  Allmetrics0<-merge(Phylosor.phylo,melt.MNTD,by=c("To","From"))
  Allmetrics<-merge(Allmetrics0,sorenson,by=c("To","From"))
  
  return(Allmetrics)}

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

##############################################################################
#Calculate Phylogenetic and Trait Diversity with respect to Taxonomic Diversity
##############################################################################

####################
# Define function
#####################

Phylo_TraitN<-function(p){
  #sink("/home1/02443/bw4sz/DimDiv/DimSink.txt")
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

splist<-colnames(siteXspp)

#test function
print(Phylo_TraitN(2))

############################
#Null model for taxonomic phylogenetic and trait
############################

clusterExport(cluster, list("null_taxlists", "splist","sp.list","siteXspp","richness_sites","beta_all","tree","mon"))
Null_dataframe<-rbind.fill(clusterApply(cluster,1:10,Phylo_TraitN))

stopCluster(cluster)
#Keep desired columns, ignoring nestedness components
Null_dataframe<-Null_dataframe[,colnames(Null_dataframe) %in% c("To","From","MNTD","Phylosor.Phylo","Sorenson","Iteration") ]

#save.image("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")

#Begin Part 3
nrow(Null_dataframe)
######################################################################
#######################Null Model Analysis, defining "High and Low"
######################################################################

#################################################
#Perform Null Model on Null Iterations
#################################################

cluster<-makeCluster(4,"SOCK")
#cluster <- getMPIcluster()

##################################################################################
#Define Function
#For each pair of assemblage compare the observed metrics to the null distribution. 
###################################################################################

print(nrow(data.merge))

Null_PT<-function(x){
  
  #Select Row
  rowS<-data.merge[x,]
  
  #Grab all the iteration rows that match these richness 
  null_rows<-Null_dataframe[Null_dataframe$To==rowS$To & Null_dataframe$From==rowS$From,]
  
  null_stats<-sapply(colnames(null_rows)[!colnames(null_rows) %in% c("To","From","Iteration","Sorenson")],function(y){
    
    if(!is.finite(rowS[,y])) {return(NA)}
    #Create a distribution of null values
    test_stat<-ecdf(null_rows[,y]) (rowS[,y])
    
    if(test_stat >= .95) {return(answer<-"High")}
    if(test_stat <= .05) {return(answer<-"Low")}
    if( test_stat <= .95 & test_stat >= .05 ) {return(answer<-"Random")}
    return(answer)
  })
  out<-data.frame(t(c(To=rowS$To,From=rowS$From,null_stats)))
  return(out)
}

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

#setwd("/home1/02443/bw4sz/")
#Write to file
#write.csv(data.df,"/home1/02443/bw4sz/DimDiv/FinalData.csv")
#write.csv(data.df.null,"/home1/02443/bw4sz/DimDiv/FinalDataNull.csv")

#Or save data
#save.image("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")

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
