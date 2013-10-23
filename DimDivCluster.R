#Taxonomic, Phylogenetic and Functional Diversity in South American Hummingbirds
#Manuscript submitted to the American Naturalist

#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW
#Boris A. Tinoco, Juan L. Parra, PhD, Leone M. Brown, PhD, Gary Stiles, PhD, Jim A. McGuire, PhD, Catherine H. Graham, PhD

######################Below are all the packages needed for this work, if any of them flag an error, just say install.package( "name of package")
require(vegan)
require(parallel)
require(reshape)
require(ape)
require(picante)
require(ggplot2)
library(doSNOW)
library(foreach)
require(GGally)
require(stringr)
require(scales)
library(multicore)
library(foreach)
require(raster)
require(boot)
#library(Rmpi)

#cluster <- getMPIcluster()
 
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
Sys.setlocale('LC_ALL','C') 


##########################################################
#Read in data
##########################################################

#Load in the data
load("/home1/02443/bw4sz/DimDiv/DimDivRevision.RData")

comm<-siteXspp[,]

dim(comm)
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
null_comm<-commsimulator(comm,method="r1")
d<-as.matrix(vegdist(null_comm,binary=TRUE,upper=FALSE,diag=FALSE))
sorenson<-data.frame(melt(d)[melt(upper.tri(d))$value,],Iteration=x)
colnames(sorenson)<-c("To","From","Sorenson","Iteration")
return(sorenson)}))

#For each pair of assemblage compare the observed metrics to the null distribution. 

system.time(null_taxlists<-foreach(x=1:nrow(sorenson),.combine=rbind) %dopar% {
 
  #Select Row
  rowS<-sorenson[x,]
      
  #Grab all the iteration rows that match these richness 
  null_rows<-tax_nullM[tax_nullM$To==rowS$To & tax_nullM$From==rowS$From,]
  
  #Create a distribution of null values
  test_stat<-ecdf(null_rows[,"Sorenson"]) (rowS[,"Sorenson"])
    
    if(test_stat >= .95) answer<-"High"
    if(test_stat <= .05) answer<-"Low"
    if( test_stat <= .95 & test_stat >= .05 ) answer<-"Random"
  
return(data.frame(t(c(To=rowS$To,From=rowS$From,Sorenson_Null=answer))))
  })

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

#################################################################################
#The Null model on 219*219 assemblages is 23871 comparisons, this is too extreme to do in a reasonable time
#I propose, we just randomize those species comparisons which have greater or less than taxonomic betadiversity given the null model of richness
table(null_taxlists$Sorenson_Null)

#Get just the high and low comparisons, the number of null comparisons is equal to the number of rows in this matrix
taxHL<-null_taxlists[null_taxlists$Sorenson_Null=="High"|null_taxlists$Sorenson_Null=="Low",]
nrow(taxHL)

#Richness of the assemblages
richness_sites<-apply(siteXspp,1,sum)

#Composition of each site
sp.list<-apply(siteXspp,1,function(x){
  names(x[which(x==1)])
})
############################
#Null model for taxonomic phylogenetic and trait
############################
system.time(
Null_dataframe<-foreach(j=1:500,.combine=rbind,.packages=c("foreach","reshape","picante")) %dopar% {
  
nullPT<-lapply(1:nrow(taxHL),function(x){
print(x)
  #Richness of the two assemblages
  richness_To<-richness_sites[names(richness_sites) %in% taxHL[x,]$To]
  richness_From<-richness_sites[names(richness_sites) %in% taxHL[x,]$From]
  
  #Create slots for each community
  nullTo<-matrix(nrow=richness_To,ncol=1)
  nullFrom<-matrix(nrow=richness_From,ncol=1)
  
  #number of species overlapping
  species.To<-sp.list[names(sp.list) %in% taxHL[x,]$To][[1]]
  species.From<-sp.list[names(sp.list) %in% taxHL[x,]$From][[1]]
  overlapping<-sum(species.To %in% species.From)
  
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
    nullTo[(overlapping+1):richness_To,]<-sample(splist[!splist %in% overlap.species],richness_To-overlapping,replace=FALSE,prob=sp.prob[!names(sp.prob) %in% overlap.species])}
  
  if(!overlapping==richness_From){
    nullFrom[(overlapping+1):richness_From,]<-sample(splist[!splist %in% overlap.species],richness_From-overlapping,replace=FALSE,prob=sp.prob[!names(sp.prob) %in% overlap.species])}
  
  #Create siteXspp table for new assemblages
  null_a<-melt(list(nullTo,nullFrom))
  
  #create binary matrix
  null_siteXspp<-(!is.na(cast(null_a,L1~value)[,-1]))*1
  rownames(null_siteXspp)<-c(as.character(taxHL[x,]$To),as.character(taxHL[x,]$From))
  
  #Compute phylogenetic and trait metrics on null assemblages
  null_beta_metrics<-beta_all(comm=null_siteXspp,tree=tree,traits=mon)
})

nullPT<-rbind.fill(nullPT)
return(data.frame(nullPT,Iteration=j))
})


save.image("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")

#Keep desired columns, ignoring nestedness components
Null_dataframe<-Null_dataframe[,colnames(Null_dataframe) %in% c("To","From","MNTD","Phylosor.Phylo","Sorenson","Iteration") ]

nrow(Null_dataframe)
######################################################################
#######################Null Model Analysis, defining "High and Low"
######################################################################

#Which rows of observed data are in the null model data
data.cNull<-data.merge[,]

index_datamerge<-paste(data.merge$To,data.merge$From,sep=".")
index_taxHL<-paste(taxHL$To,taxHL$From,sep=".")

data.cNull<-data.merge[index_datamerge %in% index_taxHL,]

#################################################
#Perform Null Model on Clusters
#################################################
#See ParallelRandomization.R for this script, contact the author
#Not included since it is specific the particular supercomputing cluster used, and not transferable.

#For each pair of assemblage compare the observed metrics to the null distribution. 

print(nrow(data.cNull))

null_PhyloTrait<-foreach(x=1:nrow(data.cNull),.combine=rbind) %dopar% {
   
  #Select Row
  rowS<-data.cNull[x,]
  
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
  print(out)
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

##########################################################################################
#Tables and Statistics
#########################################################################################

setwd("/home1/02443/bw4sz/DimDiv")
dir.create("/home1/02443/bw4sz/DimDivResults")
setwd("DimDivResults")
#Get the bounds of each 
range_metrics<-list()

for(x in 0:2){
  print(x)
  range_min<-aggregate(data.df.null[,12+x],list(data.df.null[,15+x]),min,na.rm=TRUE)
  range_max<-aggregate(data.df.null[,12+x],list(data.df.null[,15+x]),max,na.rm=TRUE)
  
  range_val<-data.frame(Index=range_min[,1],Min=range_min[,2],Max=range_max[,2])
  range_metrics[[x+1]]<-range_val
}

names(range_metrics)<-colnames(data.df)[12:14]
range_metrics<-melt(range_metrics)

write.csv(range_metrics,"Range_Metrics.csv")

#Find Prevalence of each combination
data_prev<-lapply(colnames(data.df.null)[15:17],function(x){
  range_prev<-table(data.df.null[,x])/(nrow(comm)*(nrow(comm)-1))/2})

names(data_prev)<-colnames(data.df.null)[15:17]
data_prev<-melt(data_prev)
data_prev<-cast(data_prev,L1~Var.1)
rownames(data_prev)<-data_prev[,1]
data_prev<-data_prev[,-1]

write.csv(round(data_prev,3)*100,"NullPrevalence.csv")

##################################################################
#######################ScatterPlots###############################
##################################################################

#Phylo phylosor and Hulls
p<-ggplot(data.df,aes(y=MNTD,x=Phylosor.Phylo,col=Elev)) + geom_point() 
p<-p+ theme_bw() + scale_color_gradient("Elevation",low="gray90",high="black")
p<-p+ ylab("Trait Convex Hull") + xlab("Phylogenetic Phylosor") + coord_equal()
p
ggsave("PhylosorPhylovConvexHull_Elevation.svg",height=7,width=7.5,dpi=300)


##########################################
#Plot each of the metrics with their ranges
##########################################
range_plots<-lapply(12:14,function(x){
  print(colnames(data.df.null)[x])
  print(colnames(data.df.null)[x+3])
  p<-ggplot(data.df.null,aes(y=data.df.null[,colnames(data.df.null)[x]],x=data.df.null[,colnames(data.df.null)[x+3]])) + geom_boxplot()
  p<-p+labs(y=colnames(data.df.null)[x],x=colnames(data.df.null)[x+3])
  filnam<-paste(colnames(data.df.null)[x],"_range.jpeg")
  ggsave(filnam,plot=p,height=7,width=5,dpi=300)
  return(p)})

###################################################
#Combinations of the dimensions of betadiversity
###################################################

#Create a function that takes the input of which null metrics you want to use to create output lists
#Create multiple options for the hyplist, hold them in a list and spit them to file seperately

Hyplist.func<-function(Tax,Phylo,Func){
  setwd("/home1/02443/bw4sz/DimDiv/DimDivResults")
  #Create directory
  dir.store<-dir.create(paste(Tax,Phylo,Func,sep="_"))
  setwd(paste(Tax,Phylo,Func,sep="_"))
  
  #Delinate the combinations of betadiversity
  Hyplist<-split(data.df.null,list(data.df.null[,Tax],data.df.null[,Phylo],data.df.null[,Func]))
  
  allL<-names(Hyplist)
  
  Hyplist<-Hyplist[sapply(Hyplist,nrow) > 0]
  remove.level<-names(Hyplist)[str_detect(names(Hyplist),"Random")]
  Hyplist<-Hyplist[!names(Hyplist) %in% remove.level]
  HypBox<-melt(Hyplist, id.vars=c("To","From"),measure.vars=c("Euclid","CostPathCost","AnnualPrecip","Elev"))
  
  #Proportion of combinations for each null model, later make a recursive object that grabs this output for plotting
  prop.Hyp<-melt(sapply(Hyplist,nrow)/nrow(data.df))*100
  prop.Hyp$Hyp<-as.character(rownames( prop.Hyp))
  write.csv(prop.Hyp,"ProportionHypotheses.csv")
  
  ###########################As i see it we need two randomization test, to say the difference in elevation is "high" we should compare that to the dataset as a whole.
  #To say that env differences are between hypothesis, we need to swap those labels.
  #1) Are the environmental differences "high" compared to the dataset?
  #Step 1, get the actual mean values for the env variables for each hypothesis
  
  true.median<-sapply(Hyplist,function(x){
    apply(x[,colnames(x) %in% c("Euclid","CostPathCost","AnnualPrecip","Elev")],2,median,na.rm=TRUE)})
  true.median<-as.data.frame(melt(true.median))
  colnames(true.median)<-c("Env","Hyp","value")
  
  #Step 2, get the number of cases in each hypothesis
  todraw<-sapply(Hyplist,nrow)
  
  #The goal is to draw this many rows from the entire dataset at random
  #create parallel cluster
  
  #run 1000 iterations
  boot.run<-foreach(x=1:1000,.export="data.df") %dopar% {
    
    require(reshape)
    require(boot)
    
    #drawn the desired number of rows
    boot.m<-lapply(as.numeric(todraw),function(x){ 
      rowstodraw<-sample(rownames(data.df),x)
      boot.rows<-data.df[rowstodraw,]
      
      #find median value for each env variables
      cbind(apply(boot.rows[,colnames(boot.rows) %in% c("Euclid","CostPathCost","AnnualPrecip","Elev")],2,median,na.rm=T))
    })
    names(boot.m)<-names(Hyplist)
    boot.mean<-melt(boot.m)[,-2]
    colnames(boot.mean)<-c("Env","Value","Hyp")
    list(boot.mean)}
  
  boot.run<-melt(boot.run,id.vars=c("Env","Hyp"))
  
  #Step 3, split boots and true into each combination of Env, Hyp
  split.boot<-split(boot.run,list(boot.run$Env,boot.run$Hyp))
  
  #For each position in the list, i'd like to know the distribution, and whether the true value falls outside the 95% interval. 
  p.values<-lapply(split.boot,function(x){ 
    
    CI<-boot.ci(boot(x$value,function(j,i) median(j[i],na.rm=T), R=1000)) 
    upper<-1-ecdf(x$value) (true.median[true.median$Env==levels(factor(x$Env)) & true.median$Hyp==levels(factor(x$Hyp)),"value"])
    lower<-ecdf(x$value) (true.median[true.median$Env==levels(factor(x$Env)) & true.median$Hyp==levels(factor(x$Hyp)),"value"])
    
    out<-data.frame(levels(factor(x$Env)),levels(factor(x$Hyp)),upper,lower,CI$basic[[4]],CI$basic[[5]])
    colnames(out)<-c("Env","Hyp","Upper","Lower","CIupper","CIlower")
    rownames(out)<-NULL
    return(out)
  })
  
  #now merge the data with the true median
  mean_Env_CI<-merge(rbind.fill(p.values),true.median,by=c("Env","Hyp"))
  mean_Env_CI<-mean_Env_CI[!mean_Env_CI$Env %in% c("Biome","H_mean","Tree","AnnualTemp"),]
  
  #just get the data we use for the paper precip, elev, H-mean, Biome, Euclid, CostPath
  ggplot(data=mean_Env_CI) + geom_point(aes(x=Hyp,y=value,col=Upper < .05 | Lower < .05),size=5) + facet_wrap(~ Env,scales="free_y",nrow=2) + theme_bw() + xlab("Combination of Betadiversity") + ylab("Median") + scale_color_discrete("Significant (.05)") + theme(axis.text.x = element_text(angle = 90))
  ggsave("Appendix2.jpeg",dpi=300,height=6,width=8)
  
  #Now you have a data frame with the difference in means for all reps and the true stat at the end
  #Our goal is to figure out if its within the 95 CI interval and then cast it as a matrix
  
  ######################################
  #Figure Creation and Export
  ######################################
  
  for(x in 1:length(Hyplist)){
    n<-names(Hyplist[x])
    write.csv(Hyplist[[x]],paste(n,"rowSwap.csv"),row.names=FALSE)}
  
  #####################################
  #Boxplots
  #####################################
  
  #Create Boxplots for all variables across all hypothesis and entire dataset
  dir.create("2wayplots")
  setwd("2wayplots")
  
  plots.hold<-list()
  for (x in 1:length(levels(HypBox$variable))){
   
    p<-qplot(data=HypBox[HypBox$variable==levels(HypBox$variable)[[x]],],x=" ",y=value, geom="boxplot") + facet_grid(.~L1,scale="free_x") + theme_bw() + ylab(paste("Dissim:",names(Evar[x]))) + xlab("") + geom_hline(aes(yintercept=median(data.df[,levels(HypBox$variable)[x]],na.rm=TRUE)),col='red')
    ggsave(plot=p,width=12,height=7,paste(levels(HypBox$variable)[[x]],"jpeg",sep="."))
    plots.hold[[x]]<-p
  }
  
  #Show three way combinations only.
  #Remove all random combinations
  remove.level<-levels(as.factor(HypBox$L1))[str_detect(levels(as.factor(HypBox$L1)),"Random")]
  HypBox<-HypBox[!HypBox$L1 %in% remove.level,]
  
  setwd("/home1/02443/bw4sz/DimDiv/DimDivResults")
  setwd(paste(Tax,Phylo,Func,sep="_"))
  dir.create("3WayBoxplots")
  setwd("3WayBoxplots")
  
  ###############################################
  #To do, remove outlier values from CostPath!!
  ###############################################
  
  head(HypBox)
  HypBox[HypBox$variable=="CostPathCost" & HypBox$value > 1.5e9,"value"]<-1.5e9
  
 #Get the medians
  colnames(HypBox)<-c("To","From","Diss","value","Hyp")
  intercepts<-aggregate(HypBox$value,list(HypBox$Diss),median,na.rm=T)
  colnames(intercepts)<-c("Diss","value")
  
  #Plot all simultanously, need to get intercepts on the plot?
  #add in empty levels?
  #Which levels are we missing
  levels_add<-allL[!allL %in% levels(as.factor(HypBox$Hyp)) & !allL %in% allL[str_detect(allL,"Random")]]
  
  HypBox$Hyp<-factor(HypBox$Hyp,levels=c("Low.Low.Low","High.Low.Low","High.High.Low","Low.High.Low","High.High.High","Low.High.High","Low.Low.High","High.Low.High"))
  
  p<-ggplot(HypBox,aes(x=Hyp,y=value)) + geom_boxplot()+ facet_wrap(~Diss,scales="free",drop=FALSE) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(data=intercepts,aes(yintercept=value,group=Diss), linetype="dashed",col='grey40') + theme_bw() + scale_x_discrete(drop=FALSE)
  p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave("Env3Boxplots.svg",dpi=300,height=8,width=12)
  ggsave("Env3Boxplots.jpeg",dpi=300,height=8,width=12)
  
  
  #save data from that run
  setwd("/home1/02443/bw4sz/DimDiv/DimDivResults")
  setwd(paste(Tax,Phylo,Func,sep="_"))
  
#Maps need rgdal, can't be done here
  save.image("plotting.Rdata")
}

#Run the plotting function for all sets of hypothesis

system.time(Hyplist.func(Tax="Sorenson_Null",Phylo="Phylosor.Phylo_Null",Func="MNTD_Null"))

save.image("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")
stopCluster(cluster)
q()
n
exit
