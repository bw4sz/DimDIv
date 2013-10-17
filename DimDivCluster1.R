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

##########################################################
#Read in data
##########################################################

#Load in the data
load("/home1/02443/bw4sz/DimDiv/DimDivRevision.RData")

comm<-siteXspp[,]

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
tax_nullM<-rbind.fill(lapply(1:500,function(x){
null_comm<-commsimulator(comm,method="r0")
d<-as.matrix(vegdist(null_comm,binary=TRUE,upper=FALSE,diag=FALSE))
sorenson<-data.frame(melt(d)[melt(upper.tri(d))$value,],Iteration=x)
colnames(sorenson)<-c("To","From","Sorenson","Iteration")
return(sorenson)}))

#For each pair of assemblage compare the observed metrics to the null distribution. 

nrow(tax_nullM)
nrow(sorenson)
head(sorenson)


null_taxlists<-foreach(x=1:nrow(sorenson),.combine=rbind) %dopar% { 
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

#####################################################
#Perform Phylogenetic and Trait Betadiversity
#####################################################

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
  source("/home1/02443/bw4sz/DimDiv/BenHolttraitDiversity.R")
  
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

save.image("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")
###
#End Part 1