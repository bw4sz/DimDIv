#Source Cluster Functions for DimDiv Cluster



######################################################
#Create a function for computing betadiversity metrics
#######################################################

#Define Betadiversity function

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
  #Trait Betadiversity
  #############################################
  
  #Trait frame needs to match siteXSpp table
  mon_cut<-traits[rownames(traits) %in% colnames(comm),]
  
  #There are eight species without traits, take them out for just this portion of the analysis, keep the assemblage lsit
  siteXspp_traits<-comm[,colnames(comm) %in% rownames(mon_cut)]
  
  #Zscores, standardized by sd and subtracted means
  means<-apply(mon_cut,2,mean)
  
  Bill<-mon_cut$Bill - means["Bill"]/sd(mon_cut$Bill)
  Mass<-mon_cut$Mass - means["Mass"]/sd(mon_cut$Mass)
  WingChord<-(mon_cut$WingChord - means["WingChord"])/sd(mon_cut$WingChord)
  
  z.scores<-data.frame(Bill,Mass,WingChord)
  rownames(z.scores)<-rownames(mon_cut)
  
  newSGdist<-dist(z.scores,method="euclidean")
  
  #If you wanted a PCA 
  #prc_traits<-prcomp(mon_cut)
  #newSGdist <- dist(prc_traits$x)
  
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


####################
# Define function to compare observed phylogenetic and trait betadiversity to a null model maintaining tax betadiveristy
#####################

Phylo_TraitN<-function(p){
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
  out<-data.frame(nullPTframe,Iteration=p)
  #write.csv(data.frame(nullPTframe,Iteration=p),paste("/home1/02443/bw4sz/DimDiv/Iterations/",paste(p,"iteration.csv",sep="_"),sep=""))
  return(out)
}

###########
#Function to compare taxonomic betaiversity  to null maintaining richness
###########

TaxN<-function(r){
  null_comm<-commsimulator(comm,method="r0")
  #Get new null assemblages and find sorenson
  null_t<-null_comm[rownames(null_comm) %in% as.character(r[,c("To","From")]),]
  d<-as.matrix(vegdist(null_t,upper=FALSE,diag=FALSE))
  sorenson<-melt(d)[melt(upper.tri(d))$value,]
  colnames(sorenson)<-c("To","From","Sorenson")
  return(as.matrix(sorenson))
}

print("Functions Sourced")


#Deliniate high or low betadiversity compared to null

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
