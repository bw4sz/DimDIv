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
#sink("C:/Users/Jorge/Dropbox/Shared Ben and Catherine/DimDivRevision/Results/OvernightOutput.txt")
#load data if desired

load("C:/Users/Ben/Dropbox/Shared Ben and Catherine/DimDivRevision/Results/DimDivRevision.RData")

###Define Source Functions

source("C:/Users/Jorge/Dropbox/Scripts/DimDiv/Scripts/geb12021-sup-0004-si.R.txt")
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
#load("C:/Users/Jorge/Dropbox/Shared Ben and Catherine/DimDivEntire/Output Data/Workspace.RData")

##set correct working directory to the dropbox Files for Analysis folder, whereever it is
setwd("C:\\Users\\Jorge\\Dropbox\\Shared Ben and Catherine\\DimDivEntire\\Files for Analysis")  ###Change this to the Files for Analysis folder in your dropbox, no need to move it. 

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
morph <- read.csv("C:\\Users\\Jorge\\Dropbox\\Lab paper 1 Predicted vs observed assemblages\\MorphologyShort.csv",na.strings="9999")

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
setwd("C:/Users/Jorge/Dropbox/Shared Ben and Catherine/DimDivRevision/Results/")
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
write.csv(m.test,"C:\\Users\\Jorge\\Dropbox\\Shared Ben and Catherine\\DimDivEntire\\Output Data\\Euclid_Delta_Elevmelt.csv")
Euclid_Elev_matrix<-cast(m.test,L1~L2)[,-1]

#Set rownmaes and column names in the CORRECT ORDER.
rownames(Euclid_Elev_matrix)<-loc$ID_Comm
colnames(Euclid_Elev_matrix)<-loc$ID_Comm

#Or just read it in from file
#Euclid_Elev_matrix.tocast<-read.csv("C:\\Users\\Jorge\\Dropbox\\Shared Ben and Catherine\\DimDivEntire\\Output Data\\Euclid_Delta_Elevmelt.csv")
#Euclid_Elev_matrix<-as.matrix(cast(Euclid_Elev_matrix.tocast,L1~L2)[,-1])
#rownames(Euclid_Elev_matrix)<-loc$ID_Comm
#colnames(Euclid_Elev_matrix)<-loc$ID_Comm

#append it to Evar and name it
Evar<-lappend(Evar,as.matrix(Euclid_Elev_matrix))
names(Evar)[9]<-"DeltaElevEuclid"

#Write Environmental Layers Matrices to File
setwd("C:/Users/Jorge/Dropbox/Shared Ben and Catherine/DimDivRevision/Results/")
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
save.image("C:/Users/Jorge/Dropbox/Shared Ben and Catherine/DimDivRevision/Results/DimDivRevision.RData")

#Send it to the Cluster for null model Analysis