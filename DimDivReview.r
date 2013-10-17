#Taxonomic, Phylogenetic and Functional Diversity in South American Hummingbirds
#Manuscript submitted to the American Naturalist

#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW
#Boris A. Tinoco, Juan L. Parra, PhD, Leone M. Brown, PhD, Gary Stiles, PhD, Jim A. McGuire, PhD, Catherine H. Graham, PhD

######################Below are all the packages needed for this work, if any of them flag an error, just say install.package( "name of package")
require(vegan)
require(Matrix)
require(reshape)
require(maptools)
require(raster)
require(rgdal)
require(ape)
require(picante)
require(ggplot2)
require(gdistance)
library(doSNOW)
library(foreach)
require(fields)
require(GGally)
require(stringr)
require(scales)

#load data from cluster
load("C:/Users/Ben//Dropbox/Shared Ben and Catherine/DimDivRevision/1000Iterations//DimDivRevisionCluster.RData")


#Goal of this code is to test the assumptions of the taxonomic null model for phylogenetic and trait betadiversity

#Observed Data is data.df
head(data.df)

#Null data is Null_dataframe
nrow(Null_dataframe)

#pick a row to inspect
rowS<-data.df[sample(1:nrow(data.df),1),]

#Grab all the iteration rows that match these richness 
null_rows<-Null_dataframe[Null_dataframe$To==rowS$To & Null_dataframe$From==rowS$From,]

#There should be 1000 rows
nrow(null_rows)

head(null_rows)

#Sorenson should be identical for all iterations!
table(null_rows$Sorenson)

#That's not right!



###Inspect the randomization matrix
#Pick a row
SorTest<-function(x){

#Richness of the two assemblages
richness_To<-richness_sites[names(richness_sites) %in% data.df[x,]$To]
richness_From<-richness_sites[names(richness_sites) %in% data.df[x,]$From]

#Create slots for each community
nullTo<-matrix(nrow=richness_To,ncol=1)
nullFrom<-matrix(nrow=richness_From,ncol=1)

#number of species overlapping
species.To<-sp.list[names(sp.list) %in% data.df[x,]$To][[1]]
species.From<-sp.list[names(sp.list) %in% data.df[x,]$From][[1]]
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
  nullTo[(overlapping+1):richness_To,]<-sample(splist[!splist %in% overlap.species],richness_To-overlapping,replace=FALSE)
  }

#Fill the From community, but cannot pick from any species in the filled To community

speciestopick<-splist[!splist %in% overlap.species & !splist %in% nullTo]

if(!overlapping==richness_From){
  nullFrom[(overlapping+1):richness_From,]<-sample(splist[,richness_From-overlapping,replace=FALSE)}

#Create siteXspp table for new assemblages
null_a<-melt(list(nullTo,nullFrom))

#create binary matrix
null_siteXspp<-(!is.na(cast(null_a,L1~value)[,-1]))*1
rownames(null_siteXspp)<-c(as.character(data.df[x,]$To),as.character(data.df[x,]$From))

#compute sorensons
S<-vegdist(null_siteXspp,binary=TRUE)
return(S)}


sapply(1:10,function(x){
  SorTest(1)})

#Visualize outputs
ggplot(null_rows,aes(x=Sorenson)) + geom_histogram() + geom_vline(aes(xintercept=rowS$Sorenson),col="Red")

null_stats<-sapply(colnames(null_rows)[!colnames(null_rows) %in% c("To","From","Iteration","Sorenson")],function(y){
  
  if(!is.finite(rowS[,y])) {return(NA)}
  #Create a distribution of null values
  test_stat<-ecdf(null_rows[,y]) (rowS[,y])
  
  if(test_stat >= .95) {return(answer<-"High")}
  if(test_stat <= .05) {return(answer<-"Low")}
  if( test_stat <= .95 & test_stat >= .05 ) {return(answer<-"Random")}
  ggplot(null_rows,aes(x=y)) + geom_histogram() + geom_vline(aes(xintercept=rowS[,y]),col="Red")
  
  return(answer)
})
out<-data.frame(t(c(To=rowS$To,From=rowS$From,null_stats)))
print(out)