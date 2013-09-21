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
require(FD)
require(fields)
require(GGally)
require(stringr)
require(scales)
###########################
###############Read in data
###########################

#sink output for overnight runs so we can see it tomorrow
#sink("C:/Users/Ben/Dropbox/Shared Ben and Catherine/DimDivRevision/Results/OvernightOutput.txt")
#load data if desired
load("C:/Users/Ben/Dropbox/Shared Ben and Catherine/DimDivRevision/Results/DimDivRevision.RData")

###Define Source Functions

source("C:/Users/Ben/Dropbox/Scripts/DimDiv/Scripts/geb12021-sup-0004-si.R.txt")
#Function is called beta tf


lappend <- function(lst, obj) {
  lst[[length(lst)+1]] <- obj
  return(lst)
}

#define a helpful function get species list commid MUST be in quotes!

Getsplist<-function(commID){
  names(siteXspp[commID,which(siteXspp[commID,]==1)])
}

##########################################################
#Read in data
##########################################################

#Load in the data
#load("C:/Users/Ben/Dropbox/Shared Ben and Catherine/DimDivEntire/Output Data/Workspace.RData")

##set correct working directory to the dropbox Files for Analysis folder, whereever it is
setwd("C:\\Users\\Ben\\Dropbox\\Shared Ben and Catherine\\DimDivEntire\\Files for Analysis")  ###Change this to the Files for Analysis folder in your dropbox, no need to move it. 

#Read in species matrix
siteXspp <- read.csv("siteXspp_Oct20_2011.csv", row.names=1)

#Species probabilities
sp.prob<-apply(siteXspp,2,sum)/sum(siteXspp)

#Get entire species list
splist<-colnames(siteXspp)

#Read in phylogeny
tree<-read.nexus("ColombiaPhylogenyUM.tre")

#Read in names file to replace names in Nexis file
spnames<-read.table(file="SpNameTree.txt" , sep = "\t", header = TRUE)

#Replace tip.label with Spnames#
tree$tip.label<-as.character(spnames$SpName) 

#Run PCD and split out into functional compenents, 
tree.func<-read.tree("func.tre")

#Color the func phylo by clade # this only works the 2nd time through
#col.clade<-as.vector(sapply(tree.func$tip.label,function(x) clades[clades$Dash==x,"Clade"]))
#plot.phylo(tree.func, cex=.8, tip.color=as.numeric(as.factor(col.clade)))

#bring in traits
morph <- read.csv("C:\\Users\\Ben\\Dropbox\\Lab paper 1 Predicted vs observed assemblages\\MorphologyShort.csv",na.strings="9999")

#just get males
morph.male<-morph[morph$Sex=="Macho",c("SpID","ExpC","Peso","AlCdo")]
morph.complete<-morph.male[complete.cases(morph.male),]

#aggregate for species
agg.morph<-aggregate(morph.complete,list(morph.complete$SpID),mean)
mon<-agg.morph[,-2]
colnames(mon)<-c("Species","Bill","Mass","WingChord")
rownames(mon)<-mon[,1]
mon<-mon[,-1]

#Replace spaces with underscore
rownames(mon)<-gsub(" ","_",rownames(mon))

#######################################################
#Compute Environmental Dissimilarity

##############################################################
#Step 2 Bring in Env Info and Dissimilarity 
##############################################################
comm<-siteXspp


#Extract env information from each locality
#Import Localities
loc<-readShapePoints("Locs_proj.shp")
head(loc)

#Put all raster layers into one folder, point towards that folder.
envL<-list.files(pattern=".tif$","D:\\Ben\\GIS\\Data\\DimDiv\\EnvLayers",full.names=T)
env.list<-lapply(envL,raster)
names(env.list)<-c("AnnualPrecip","Elev","H_mean","AnnualTemp","Tree")

#Get any of the shapefiles too
envP<-list.files("D:\\Ben\\GIS\\Data\\DimDiv\\EnvLayers",".shp$",full.names=T)
poly.list<-lapply(envP,readShapePoly)

#Extract raster info to points
extract.list<-lapply(env.list, function(x) extract(x,loc))
poly.df<-over(loc,poly.list[[1]])

#THIS MUST BE DONE MANUALLY and must be Correct order compared to envL! 
names(extract.list)<-c("AnnualPrecip","Elev","H_mean","AnnualTemp","Tree")
Envtable<-data.frame(loc@data,as.data.frame(extract.list))
Envtable$Biome<-poly.df$BIOME

#remember to divide temp by 10
Envtable$AnnualTemp<-Envtable$AnnualTemp/10

#check for any NA's
IDerror<-Envtable[!complete.cases(Envtable),]

###Ive checked in arcmap, and i have to believe we arent getting data from mangroves...
#so set the 293 as biome 2
#and set 495 as biome 1 
Envtable[Envtable$ID_Comm==293,"Biome"]<-2
Envtable[Envtable$ID_Comm==495,"Biome"]<-1

#Hmean error
#nearest cells looks to be ~3100, estimate by two nearest cells, by coast
Envtable[Envtable$ID_Comm==495,"H_mean"]<-3100

#tree error, its in a river! 85 on both sides of the cell. 
Envtable[Envtable$ID_Comm==194,"Tree"]<-85

#also in a river, nearest cell
Envtable[Envtable$ID_Comm==525,"Tree"]<-45
write.csv(Envtable,"Envtable.csv")

#Create distance matrices of each Environmental variable, store it in a list. 
Evar<-lapply(Envtable[7:12],function(x)(dist(x, method='maximum', diag=TRUE, upper=TRUE,)))

#Create euclidean distance matrix from the long/lat coordinates
longlat<-as.matrix(cbind(Envtable$LongDecDeg,Envtable$LatDecDeg))
rownames(longlat)<-Envtable$ID_Comm
edist<-rdist.earth(longlat,longlat,miles=FALSE)
Evar<-lappend(Evar,edist)

#Bring in Cost Path from a R gdistance script
#############iF you are not changing ANY of the above steps uncomment the next 2 lines, and SKIP this section of data
CostPathMatrix<-read.csv("CostPathCost.csv", row.names=1)
rownames(CostPathMatrix)<-Envtable$ID_Comm
colnames(CostPathMatrix)<-Envtable$ID_Comm

#Import Friction layer
elevr<-raster("studyarea_1km.tif")

#Try different elevation layer
elev.test<-crop(elevr,extent(loc)*1.2)
elev.test<-aggregate(elev.test,5)

#Find shortest cost path
cl<-makeCluster(8,"SOCK")
registerDoSNOW(cl)
costPath.list<-foreach(x = 1:length(loc),.packages=c("raster","gdistance")) %dopar% {
  
  #pick the original site. 
  orig<-loc[x,]
  #What elevation is the origin
  elev_origin<-extract(elev.test,orig)[[1]]
  
  #Get the difference between the origin elevation and every cell in the raster
  elev_diff<-abs(elev_origin-elev.test)
  
  #create a the transition matrix where permeability is high where elev difference is low
  trans_diff<-transition(elev_diff,function(x) 1/mean(x),8)
  
  #Diagonal Cell Correction, diagonals are father away than adjacent cells. 
  slope <- geoCorrection(trans_diff)
  
  #Remember this cost surface is only valid for this site as the origen, ie. we want to create just this column in the cost distance matrix
  #Cost Distance
  cdist<-costDistance(slope,orig,loc)
  #labelling is key here.
  
  return(list(cdist))}
stopCluster(cl)

#reshape cost path for each locality pair
m.costlist<-melt(costPath.list)
CostPathMatrix<-cast(m.costlist,X2~L1)[,-1]
rownames(CostPathMatrix)<-loc$ID_Comm
colnames(CostPathMatrix)<-loc$ID_Comm

#Write to file
setwd("C:/Users/Ben/Dropbox/Shared Ben and Catherine/DimDivEntire/Output Data")
write.csv(CostPathMatrix,"CostPathCost.csv")

##############################If you did skip making the costpath, start here again. 
Evar<-lappend(Evar,as.matrix(CostPathMatrix))
names(Evar)<-c(names(Envtable[7:ncol(Envtable)]),"Euclid", "CostPathCost")

#Same Ecosystem or Across Ecosystem
names(Evar)
Evar[["Biome"]][Evar[["Biome"]]>0]<-1

####Straight Line Cost: What is the total elevation accumulation between two sites?
#Warning, this function takes awhile to run. 
cl<-makeCluster(8,"SOCK")
registerDoSNOW(cl)
system.time(Euclid_DeltaElev<-foreach(x=1:length(loc),.packages=c("raster","rgdal")) %dopar% {
  lapply(1:length(loc),function(y){
    
    #Get a test lines class and make a spatial line
    test.line<-Line(rbind(loc[x,],loc[y,]))
    sp.test.line<-SpatialLines(list(Lines(test.line,ID="A")))
    
    #Extract the elevation
    #Order the elevation at get profile 
    EuclidProfile<-extract(elev.test,sp.test.line,cellnumbers=T)[[1]]
    rc <- rowColFromCell(elev.test, EuclidProfile[,1]) 
    b <- as.data.frame(cbind(rc, value=EuclidProfile[,2]) )
    
    #Get the change in elevation
    Euclid_Elev<-sum(abs(diff(b$value,1)),na.rm=T)
    return(Euclid_Elev)})
})
stopCluster(cl)

#melt and cast into a full matrix
m.test<-melt(Euclid_DeltaElev)

#easier to write this to file?
write.csv(m.test,"C:\\Users\\Ben\\Dropbox\\Shared Ben and Catherine\\DimDivEntire\\Output Data\\Euclid_Delta_Elevmelt.csv")
Euclid_Elev_matrix<-cast(m.test,L1~L2)[,-1]

#Set rownmaes and column names in the CORRECT ORDER.
rownames(Euclid_Elev_matrix)<-loc$ID_Comm
colnames(Euclid_Elev_matrix)<-loc$ID_Comm

#Or just read it in from file
#Euclid_Elev_matrix.tocast<-read.csv("C:\\Users\\Ben\\Dropbox\\Shared Ben and Catherine\\DimDivEntire\\Output Data\\Euclid_Delta_Elevmelt.csv")
#Euclid_Elev_matrix<-as.matrix(cast(Euclid_Elev_matrix.tocast,L1~L2)[,-1])
#rownames(Euclid_Elev_matrix)<-loc$ID_Comm
#colnames(Euclid_Elev_matrix)<-loc$ID_Comm

#append it to Evar and name it
Evar<-lappend(Evar,as.matrix(Euclid_Elev_matrix))
names(Evar)[9]<-"DeltaElevEuclid"

#Write Environmental Layers Matrices to File
setwd("C:\\Users\\Ben\\Dropbox\\Shared Ben and Catherine\\DimDivEntire\\Output Data")
for(x in 1:length(Evar))
  write.csv(as.matrix(Evar[[x]]),file=paste(names(Evar[x]),".csv",""))

#create giant data frame using the comparisons between all sites
Evar.m<-lapply(Evar,function(x) {
  y<-as.matrix(x)
  colnames(y)<-Envtable$ID_Comm
  rownames(y)<-Envtable$ID_Comm
  return(y)})

test.Evar<-melt(Evar.m)
compare.env<-cast(test.Evar,X1+X2~L1)
colnames(compare.env)<-c("To","From",names(compare.env[-c(1,2)]))

#save to file
save.image("C:/Users/Ben/Dropbox/Shared Ben and Catherine/DimDivRevision/Results/DimDivRevision.RData")

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
tax_nullM<-rbind.fill(lapply(1:1000,function(x){
null_comm<-commsimulator(comm,method="r1")
d<-as.matrix(vegdist(null_comm,binary=TRUE,upper=FALSE,diag=FALSE))
sorenson<-data.frame(melt(d)[melt(upper.tri(d))$value,],Iteration=x)
colnames(sorenson)<-c("To","From","Sorenson","Iteration")
return(sorenson)}))

#For each pair of assemblage compare the observed metrics to the null distribution. 

cl<-makeCluster(8,"SOCK")
registerDoSNOW(cl)
null_taxlists<-foreach(x=1:nrow(sorenson),.combine=rbind) %dopar% {
 
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
  }
stopCluster(cl)

#####################################################
#Perform Phylogenetic and Trait Betadiversity
#####################################################

beta_all<-function(comm,tree,traits){
  
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
  
  #get all pairwise combinations of sites, depending if you want a full matrix (null) or sparse matrix (observed values)
  pair.w<-combn(rownames(siteXspp_traits),2,simplify=FALSE)
    
  #loop through all pairs of assemblages and get the functional overlap
  pairwise.beta<-foreach(x=pair.w,.packages=c("vegan","reshape"),.errorhandling="pass") %do%{
    source("C:/Users/Ben/Dropbox/Scripts/DimDiv/Scripts/geb12021-sup-0004-si.R.txt")
    Villeger<-beta_TF(siteXspp_traits[rownames(siteXspp_traits) %in% x,] ,as.matrix(mon_cut))$beta
    Villeger<-melt(Villeger)
    cast(Villeger,~X1+X2)
  }
  
  #get rid of the NA rows
  toremove<-sapply(pairwise.beta,function(x) is.character(x[[1]]))
  pairwise.beta.removed<-rbind.fill(pairwise.beta[!toremove])
  
  #Get the order of inputs
  pairwise.order<-t(sapply(pair.w,function(x) {matrix(nrow=1,ncol=2,x)}))
  pairwise.order.removed<-pairwise.order[!toremove,]
  
  #If there are just two rows, it needs to be turned back to a matrix
  if(class(pairwise.order.removed)=="character"){pairwise.order.removed<-t(matrix(pairwise.order.removed))}
  func.beta<-data.frame(pairwise.order.removed,pairwise.beta.removed)[,-3]
  
  #Combine the dataframes
  colnames(func.beta)[1:2]<-c("To","From")
  
  #If functional is empty, 
  if(nrow(func.beta)==0){func.beta<-data.frame(pairwise.order,beta_taxonomic=NA,nestedness_functional=NA,nestedness_taxonomic=NA,turnover_functional=NA,turnover_taxonomic=NA)
                         colnames(func.beta)[1:2]<-c("To","From")}
  
  Allmetrics<-merge(func.beta,Phylosor.phylo,by=c("To","From"))

  return(Allmetrics)}

system.time(beta_metrics<-beta_all(comm=comm,tree=tree,traits=mon))

#Visualizations of the beta metrics
head(beta_metrics)
beta_metrics<-merge(beta_metrics,sorenson,c("To","From"))

#Get rid of the nestedness components
beta_metrics<-beta_metrics[,colnames(beta_metrics) %in% c("To","From","beta_functional","Phylosor.Phylo","Sorenson") ]

#####################################################
#Merge Betadiversity and Environmnetal Dissimilairity
#####################################################
data.merge<-merge(compare.env,beta_metrics,by=c("To","From"))

save.image("C:/Users/Ben/Dropbox/Shared Ben and Catherine/DimDivRevision/Results/DimDivRevision.RData")

#################################################################################
#The Null model on 219*219 assemblages is 23871 comparisons, this is too extreme to do in a reasonable time
#I propose, we just randomize those species comparisons which have greater or less than taxonomic betadiversity given the null model of richness
table(null_taxlists$Sorenson_Null)

#Get just the high and low comparisons, the number of null comparisons is equal to the number of rows in this matrix
taxHL<-null_taxdf[null_taxlists$Sorenson_Null=="High"|null_taxlists$Sorenson_Null=="Low",]
nrow(taxHL)

#Richness of the assemblages
richness_sites<-apply(siteXspp,1,sum)

#Composition of each site
sp.list<-apply(siteXspp,1,function(x){
  names(x[which(x==1)])
})

############################
#Null model for phylogenetic and trait
############################

cl<-makeCluster(8,"SOCK")
registerDoSNOW(cl)

Null_dataframe<-foreach(j=1:1000,.combine=rbind,.packages=c("foreach","reshape","picante")) %dopar% {
  
nullPT<-lapply(1:nrow(taxHL),function(x){
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
}

stopCluster(cl)

#Keep desired columns, ignoring nestedness components
Null_dataframe<-Null_dataframe[,colnames(Null_dataframe) %in% c("To","From","beta_functional","Phylosor.Phylo","Sorenson","Iteration") ]

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
null_lists<-list()
for (x in 1:nrow(data.cNull)){
   
  #Select Row
  rowS<-data.cNull[x,]
  
  #Grab all the iteration rows that match these richness 
  null_rows<-Null_dataframe[Null_dataframe$To==rowS$To & Null_dataframe$From==rowS$From,]
  
  null_stats<-sapply(colnames(null_rows)[!colnames(null_rows) %in% c("To","From","Iteration")],function(y){
        
    if(!is.finite(rowS[,y])) return(NA)
    #Create a distribution of null values
    test_stat<-ecdf(null_rows[,y]) (rowS[,y])
    
        if(test_stat >= .95) return(answer<-"High")
    if(test_stat <= .05) return(answer<-"Low")
    if( test_stat <= .95 & test_stat >= .05 ) return(answer<-"Random")
    return(answer)
  })
  
  null_lists[[x]]<-data.frame(t(c(To=rowS$To,From=rowS$From,null_stats)))
}

#Bind together the null model outputs
rownames(null_lists)<-NULL
null_PhyloTrait<-rbind.fill(null_lists)
colnames(null_PhyloTrait)<-c("To","From",paste(colnames(Null_dataframe)[!colnames(Null_dataframe) %in% c("To","From","Iteration")],"Null",sep="_"))

#Combine the environmental, observed and null metrics into a huge dataframe
data.df.null<-merge(data.merge,null_PhyloTrait,by=c("To","From"))
data.df.null<-merge(data.df.null,null_taxdf,by=c("To","From"))

#Legacy correction, data.merge is data.df, sorry
data.df<-data.merge

#Write to file
write.csv(data.df,"FinalData.csv")
write.csv(data.df.null,"FinalDataNull.csv")


#Or save data
save.image("C:/Users/Ben/Dropbox/Shared Ben and Catherine/DimDivRevision/Results/DimDivRevision.RData")

#Data Generation Complete
##########################################################################################
##########################################################################################

##########################################################################################
#Tables and Statistics
#########################################################################################

setwd("C:\\Users\\Ben\\Dropbox\\Shared Ben and Catherine\\DimDivRevision\\Results")

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

#PCDp Phylo and hulls
#PCD func and PCD phylo
p<-ggplot(data.df,aes(y=beta_functional,x=Phylosor.Phylo,col=Sorenson)) + geom_point() 
p<-p+ theme_bw() + scale_color_gradient("Sorenson",low="gray90",high="black")
p<-p+ ylab("Trait Betadiversity") + xlab("Phylogenetic Betadiversity") + coord_equal()
p
ggsave("PhylosorPhylovConvexHull_Taxonomic.svg",height=7,width=7.5,dpi=300)

#Phylo phylosor and Hulls
p<-ggplot(data.df,aes(y=beta_functional,x=Phylosor.Phylo,col=Elev)) + geom_point() 
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

setwd("C:\\Users\\Ben\\Dropbox\\Shared Ben and Catherine\\DimDivRevision\\Results")

Hyplist.func<-function(Tax,Phylo,Func){
  setwd("C:\\Users\\Ben\\Dropbox\\Shared Ben and Catherine\\DimDivRevision\\Results")
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
  boot.run<-foreach(x=1:1000,.export="data.df") %do% {
    
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
  
  setwd("C:\\Users\\Ben\\Dropbox\\Shared Ben and Catherine\\DimDivRevision\\Results")
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
  
  #Draw Lines between all hypothesis one sites
  elevr<-raster("C:\\Users\\Ben\\Dropbox\\Shared Ben and Catherine\\DimDivEntire\\Files for Analysis\\studyarea_1km.tif")
  
  #Function for spatial lines for all hypothesis, overlayed on an Ecuador Map
  
  #ugly function to get to the directory level, sorry. 
  setwd("C:\\Users\\Ben\\Dropbox\\Shared Ben and Catherine\\DimDivRevision\\Results")
  setwd(paste(Tax,Phylo,Func,sep="_"))
  dir.create("Maps")
  setwd("Maps")
  
  #only Create maps of three way comparisons
 
  #Lines for each hypothesis, plotted individually
  #plot all combinations
  map_names<-c("Low.Low.Low","High.Low.Low","High.High.Low","Low.High.Low","High.High.High","Low.High.High","Low.Low.High","High.Low.High")
  
  cl<-makeCluster(8,"SOCK")
  maps.Hyp<-foreach(f=1:length(map_names),.packages=c("reshape","raster","rgdal")) %do% {
    jpeg(paste(map_names[f],"RowSwap.jpeg",sep=""),height= 10, width = 10, unit="in", res=300)
    if(length(Hyplist[[map_names[f]]])>0){
      coords<-apply(Hyplist[[map_names[f]]],1,function(x){
        x_stat<-Envtable[Envtable$CommID %in% as.numeric(x["To"]),c("LatDecDeg","LongDecDeg")]
        y_stat<-Envtable[Envtable$CommID %in% as.numeric(x["From"]),c("LatDecDeg","LongDecDeg")]
        comb<-data.frame(x_stat,y_stat)
        colnames(comb)<-c("xmin","ymin","xmax","ymax")
        return(comb)
      })
      cordmatrix<-rbind.fill(coords)
    }
    plot(elevr, axes=FALSE,main="",legend=F,col=colorRampPalette(c("grey90", "grey6"))(255))
    if(length(Hyplist[[map_names[f]]])>0){
      for (j in 1:nrow(cordmatrix)) {
        arrows(y0=cordmatrix[j,1],x0=cordmatrix[j,2],y1=cordmatrix[j,3],x1=cordmatrix[j,4],length=0, lwd=1,col="grey20")}
    }
    points(loc, col='grey20', pch=15,cex=.5)
    dev.off()
  }
  stopCluster(cl)
  
  
  # Together on one graph, if needed.
  jpeg("AllhypothesisRowSwap.jpeg",quality=100,res=300,units="in",height=12,width=20)
  par(mfrow=c(2,4))
  
  #Plot lines
  
  for (x in 1:length(Hyplist)){
    coords<-apply(Hyplist[[x]],1,function(x){
      x_stat<-Envtable[Envtable$CommID %in% as.numeric(x["To"]),c("LatDecDeg","LongDecDeg")]
      y_stat<-Envtable[Envtable$CommID %in% as.numeric(x["From"]),c("LatDecDeg","LongDecDeg")]
      comb<-cbind(x_stat,y_stat)
      colnames(comb)<-c("xmin","xmax","ymin","ymax")
      return(comb)
    })
    cordmatrix<-rbind.fill(coords)
    
    #Plot map and create arrows between assemblages
    plot(elevr, axes=FALSE,main=names(Hyplist)[x],legend=F,col=colorRampPalette(c("grey90", "grey6"))(255))
    points(loc, col='red', pch=10,cex=.5)
    for (x in 1:nrow(cordmatrix)) {
      arrows(y0=cordmatrix[x,1],x0=cordmatrix[x,2],y1=cordmatrix[x,3],x1=cordmatrix[x,4],length=0, lwd=1)
      
    }
  }
  dev.off()
  
  
  ###Done
  #save data from that run
  setwd("C:\\Users\\Ben\\Dropbox\\Shared Ben and Catherine\\DimDivRevision\\Results")
  setwd(paste(Tax,Phylo,Func,sep="_"))
  
  save.image("plotting.Rdata")
}

#Run the plotting function for all sets of hypothesis

system.time(Hyplist.func(Tax="Sorenson_Null",Phylo="Phylosor.Phylo_Null",Func="beta_functional_Null"))

