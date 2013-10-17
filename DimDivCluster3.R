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
load("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")

nrow(Null_dataframe)
######################################################################
#######################Null Model Analysis, defining "High and Low"
######################################################################
#################################################
#Perform Null Model on Clusters
#################################################

#For each pair of assemblage compare the observed metrics to the null distribution. 

print(nrow(data.merge))

null_PhyloTrait<-foreach(x=1:nrow(data.merge),.combine=rbind) %dopar% {
   
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
q()
n
exit
