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

###########################
###############Read in data
###########################

#load data if desired
load("C:/Users/Jorge/Dropbox/Shared Ben and Catherine/DimDivRevision/Results/DimDivRevision.rData.RData")

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

######################################################
#Create a function for computing betadiversity metrics
#######################################################

##########################
# test matrix and traits #
comm<-siteXspp[1:100,]
#traits<-mon     
##########################

beta_all<-function(comm,tree,traits,cores,FullMatrix){

#Phylogenetic PCD
PCD.phylo<-pcd(comm,tree)
PCD.phylo<-melt(lapply(PCD.phylo,as.matrix))
colnames(PCD.phylo)<-c("To","From","PCD.value.phylo","PCD.phylo")

#Phylosor Calculation see Bryant 2008
Phylosor.phylo<-melt(as.matrix(phylosor(comm,tree)))
colnames(Phylosor.phylo)<-c("To","From","Phylosor.Phylo")

#Merge phylometrics
phylometrics<-merge(Phylosor.phylo,PCD.phylo,c("To","From"))

#####################################
##Taxonomic Betadiversity
sorenson<-melt(as.matrix(vegdist(comm,binary=TRUE)))
colnames(sorenson)<-c("To","From","Sorenson")
######################################

#Phylogenetic PCD
PCD.func<-pcd(comm,tree=tree.func)
PCD.func<-melt(lapply(PCD.func,as.matrix))
colnames(PCD.func)<-c("To","From","PCD.value.func","PCD.func")

#Phylosor Calculation see Bryant 2008
Phylosor.func<-melt(as.matrix(phylosor(comm,tree.func)))
colnames(Phylosor.func)<-c("To","From","Phylosor.Func")

#Merge phylometrics
funcmetrics<-merge(Phylosor.func,PCD.func,c("To","From"))

#merge this into the phylometrics
Phylo_Tax2<-merge(phylometrics,funcmetrics,c("To","From"))

#################
#Trait Metrics
#################
Phylo_Tax<-merge(Phylo_Tax2,sorenson,c("To","From"))

#############################################
#Non-dendrogram approach, functional approach employed by villeger 2013
#############################################

#Trait frame needs to match siteXSpp table
mon_cut<-traits[rownames(traits) %in% colnames(comm),]

#There are eight species without traits, take them out for just this portion of the analysis, keep the assemblage lsit
siteXspp_traits<-comm[,colnames(comm) %in% rownames(mon_cut)]
                  
#get all pairwise combinations of sites, depending if you want a full matrix (null) or sparse matrix (observed values)
if(FullMatrix==FALSE){pair.w<-combn(rownames(siteXspp_traits),2,simplify=FALSE)}
if(FullMatrix==TRUE){
  #Get the combinations of the null model, we only want  "A" compared to "B"
  pair.w<-expand.grid(rownames(null.siteXspp.matrix)[1:length(richness_levels)],rownames(null.siteXspp.matrix)[(length(richness_levels)+1):(length(richness_levels)*2)])
  pair.w<- as.list(as.data.frame(t(pair.w)))
}

#loop through all pairs of assemblages and get the functional overlap
pairwise.beta<-foreach(x=pair.w,.packages=c("vegan","reshape"),.errorhandling="pass") %do%{
  source("C:/Users/Jorge/Dropbox/Scripts/DimDiv/Scripts/geb12021-sup-0004-si.R.txt")
  Villeger<-beta_TF(siteXspp_traits[rownames(siteXspp_traits) %in% x,] ,as.matrix(mon_cut))$beta
  Villeger<-melt(Villeger)
  cast(Villeger,~X1+X2)
  }

toremove<-sapply(pairwise.beta,function(x) is.character(x[[1]]))

#get rid of the NA rows
toremove<-sapply(pairwise.beta,function(x) is.character(x[[1]]))
pairwise.beta.removed<-rbind.fill(pairwise.beta[!toremove])

#Get the order of inputs
pairwise.order<-t(sapply(pair.w,function(x) {matrix(nrow=1,ncol=2,x)}))
pairwise.order.removed<-pairwise.order[!toremove,]

#Combine the dataframes
func.beta<-data.frame(pairwise.order.removed,pairwise.beta.removed)[,-3]
colnames(func.beta)[1:2]<-c("To","From")
Allmetrics<-merge(func.beta,Phylo_Tax,by=c("To","From"))

#Cast out the full frame
require(reshape2)

Allmetrics1<-dcast(Allmetrics,...~PCD.phylo,value.var="PCD.value.phylo")
colnames(Allmetrics1)[14:16]<-c("PCD.phylo","PCDc.phylo","PCDp.phylo")

Allmetrics2<-dcast(Allmetrics1,...~PCD.func,value.var="PCD.value.func")
colnames(Allmetrics2)[15:17]<-c("PCD.func","PCDc.func","PCDp.func")

return(Allmetrics2)}

system.time(beta_metrics<-beta_all(comm=comm,tree=tree,traits=mon,FullMatrix=FALSE))
                            
#Visualizations of the beta metrics
head(beta_metrics)

#################################################################################
#################################################################################
#This is just an idea, but it seems duplicatous to draw this huge matrix with replacment 219*219
#What makes more sense is just visualize the type of assemblages we have

#Any assemblage drawing from a comparison of richness = 6 in A and richness =8 in B is drawn from the same distribution
#Therefore just simulate one distribution for each unique type of assemblage comparison
richness_sites<-apply(siteXspp,1,sum)
richness_levels<-as.numeric(names(table(apply(siteXspp,1,sum))))

paste("Iteration is:",x)

null.siteXspp.matrix<-matrix(nrow=length(richness_levels)*2,ncol=length(splist))
rownames(null.siteXspp.matrix)<-rep(richness_levels,2)
colnames(null.siteXspp.matrix)<-splist

#Insert the correct numbers of 0,1 for each row
for(x in 1:nrow(null.siteXspp.matrix)){
  rich<-as.numeric(rownames(null.siteXspp.matrix)[x])
  null.siteXspp.matrix[x,]<-sample(c(rep(0,length(splist)-rich),rep(1,rich)),replace=FALSE)
}

#for the foreach function there can't be duplicate rownmaes, set the first as A and the 2nd as B
rownames(null.siteXspp.matrix)[1:length(richness_levels)]<-paste(rownames(null.siteXspp.matrix)[1:length(richness_levels)],"A")
rownames(null.siteXspp.matrix)[(length(richness_levels)+1):(length(richness_levels)*2)]<-paste(rownames(null.siteXspp.matrix)[(length(richness_levels)+1):(length(richness_levels)*2)],"B")

#Compute null distributions for each combination of taxonomic diversity, use FullMatrix=TRUE
cl<-makeCluster(8,"SOCK")
registerDoSNOW(cl)


system.time(null_models<-foreach(x=1:50,.combine=rbind,.packages=c("vegan","picante","reshape","foreach")) %dopar% { 
  print(x)
  null.matrix<-commsimulator(null.siteXspp.matrix,"r0")
  null_beta_metrics<-beta_all(comm=null.matrix,tree=tree,traits=mon,FullMatrix=TRUE)
  return(data.frame(null_beta_metrics,Iteration=x))
})
stopCluster(cl)

#Null model for beta metrics, this is just here for sample, the real code will need a cluster!

##############################################################
#Step 2 Bring in Env Info and Dissimilarity 
##############################################################

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
elev.test<-aggregate(elev.test,10)

#Find shortest cost path
cl<-makeCluster(8,"SOCK")
registerDoSNOW(cl)
costPath.list<-foreach(x = 1:length(loc),.packages=c("raster","gdistance")) %dopar% {
  print(x)
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
setwd("C:/Users/Jorge/Dropbox/Shared Ben and Catherine/DimDivEntire/Output Data")
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
setwd("C:\\Users\\Jorge\\Dropbox\\Shared Ben and Catherine\\DimDivEntire\\Output Data")
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

#####################################################
#Merge Betadiversity and Environmnetal Dissimilairity
#####################################################
data.merge<-merge(compare.env,beta_metrics,by=c("To","From"))

######################################################################
#######################Null Model Analysis
######################################################################

#################################################
#Perform Null Model on Clusters
#################################################
#See ParallelRandomization.R for this script, contact the author
#Not included since it is specific the particular supercomputing cluster used, and not transferable.

#For each pair of assemblage compare the observed metrics to the null distribution. 
null_lists<-list()
for (x in 1:nrow(data.merge)){
print(x)
#Select Row
rowS<-data.merge[x,]
  
#What is the richness of each assemblage in this comparison?
#Paste A and B on that to get the row names of the null model dataframe
richness_To<-paste(richness_sites[names(richness_sites) %in% rowS$To],"A")
richness_From<-paste(richness_sites[names(richness_sites) %in% rowS$From],"B")

#Grab all the iteration rows that match these richness 
null_rows<-null_models[null_models$To==richness_To & null_models$From==richness_From,]

null_stats<-sapply(colnames(null_rows)[-c(1,2,18)],function(y){

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
null_dataframe<-rbind.fill(null_lists)
colnames(null_dataframe)<-c("To","From",paste(colnames(null_dataframe)[-c(1,2)],"Null",sep="_"))

#Combine the environmental, observed and null metrics into a huge dataframe
data.df<-merge(data.merge,null_dataframe,by=c("To","From"))

#Or save data
save.image("C:/Users/Jorge/Dropbox/Shared Ben and Catherine/DimDivEntire/Output Data/Workspace.RData")

#Data Generation Complete
##########################################################################################
##########################################################################################

##########################################################################################
#Tables and Statistics
#########################################################################################

setwd("C:\\Users\\Jorge\\Dropbox\\Shared Ben and Catherine\\DimDivRevision\\Results")

#Get the bounds of each 
range_metrics<-list()

for(x in 0:14){
  print(x)
range_min<-aggregate(data.df[,12+x],list(data.df[,27+x]),min,na.rm=TRUE)
range_max<-aggregate(data.df[,12+x],list(data.df[,27+x]),max,na.rm=TRUE)

range_val<-data.frame(Index=range_min[,1],Min=range_min[,2],Max=range_max[,2])
range_metrics[[x+1]]<-range_val
}

names(range_metrics)<-colnames(data.df)[12:26]
range_metrics<-melt(range_metrics)

write.csv(range_metrics,"Range_Metrics.csv")


#Find Prevalence of each combination
data_prev<-lapply(colnames(data.df)[27:40],function(x){
range_prev<-table(data.df[,x])/nrow(data.df)})

names(data_prev)<-colnames(data.df)[27:40]
data_prev<-melt(data_prev)
data_prev<-cast(data_prev,L1~Var.1)
rownames(data_prev)<-data_prev[,1]
data_prev<-data_prev[,-1]

write.csv(round(data_prev,3)*100,"NullPrevalence.csv")

###################################################
#correlation and comparisons

#Taxonomic
Tax_plot<-ggpairs(data.df[,c("beta_taxonomic","PCDc.phylo","Sorenson")]) 
Tax_plot
ggsave("Tax_plot.svg")

#Phylogenetic
Phylo_plot<-ggpairs(data.df[,c("Phylosor.Phylo","PCDp.phylo")]) 
Phylo_plot
ggsave("Phylo_plot.svg")

#Functional
Func_plot<-ggpairs(data.df[,c("Phylosor.Func","PCDp.func","beta_functional")]) 
Func_plot
ggsave("Func_plot.svg")

#Plot each of the metrics with their ranges

range_plots<-lapply(12:26,function(x){
  print(colnames(data.df)[x])
print(colnames(data.df)[x+15])
p<-ggplot(data.df,aes(y=data.df[,colnames(data.df)[x]],x=data.df[,colnames(data.df)[x+15]])) + geom_boxplot()
p<-p+labs(y=colnames(data.df)[x],x=colnames(data.df)[x+15])
filnam<-paste(colnames(data.df)[x],"_range.jpeg")
ggsave(filnam,plot=p,height=7,width=5,dpi=300)
return(p)})

###################################################
#Combinations of the dimensions of betadiversity
###################################################

#Create a function that takes the input of which null metrics you want to use to create output lists

#Create multiple options for the hyplist, hold them in a list and spit them to file seperately

#For testing

Tax<-"Sorenson_Null"
Phylo<-"PCDp.phylo_Null"
Func<-"beta_functional_Null"

#Input as names of columns

Hyplist.func<-function(Tax,Phylo,Func){
  
  #Create directory
  dir.store<-dir.create(paste(Tax,Phylo,Func,sep="_"))
  setwd(paste(Tax,Phylo,Func,sep="_"))
  
  #Delinate the combinations of betadiversity
  Hyplist<-split(data.df,list(data.df[,Tax],data.df[,Phylo],data.df[,Func]))
  
  #Proportion of combinations for each null model, later make a recursive object that grabs this output for plotting
  prop.Hyp<-melt(sapply(Hyplist,nrow)/nrow(data.df))*100
  prop.Hyp$Hyp<-as.character(rownames( prop.Hyp))
  write.csv(prop.Hyp,"ProportionHypotheses.csv")
  
  Hyplist<-Hyplist[sapply(Hyplist,nrow) > 1]
  HypBox<-melt(Hyplist, id.vars=c("To","From"),measure.vars=c("Euclid","CostPathCost","AnnualPrecip","Elev"))
  
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
cl<-makeCluster(8,"SOCK")
registerDoSNOW(cl)

#run 1000 iterations
boot.run<-times(1000) %dopar% {
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
stopCluster(cl)

boot.run<-melt(boot.run)

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
ggsave("\Appendix2.jpeg",dpi=300,height=6,width=8)

#Now you have a data frame with the difference in means for all reps and the true stat at the end
#Our goal is to figure out if its within the 95 CI interval and then cast it as a matrix

######################################
#Figure Creation and Export
######################################

#Write hypothesis list to files 
setwd("C:\\Users\\Jorge\\Dropbox\\Shared Ben and Catherine\\DimDivEntire\\Output Data")

for(x in 1:length(Hyplist.raw)){
  n<-names(Hyplist.raw[x])
  write.csv(Hyplist.raw[[x]],paste(n,"rowSwap.csv"),row.names=FALSE)}

#Create Boxplots for all variables across all hypothesis and entire dataset
setwd("C:\\Users\\Jorge\\Dropbox\\Shared Ben and Catherine\\DimDivEntire\\BoxplotsEntire")
for (x in 1:length(Evar)){
  qplot(data=HypBox[HypBox$variable==names(Evar[x]),],x=" ",y=value, geom="boxplot") + facet_grid(.~L1,scale="free_x") + theme_bw() + ylab(paste("Dissim:",names(Evar[x]))) + xlab("") + geom_hline(aes(yintercept=median(data.d[,names(Evar[x])],na.rm=TRUE)),col='red')
  ggsave(width=12,height=7,paste(names(Evar[x]),"jpeg",sep="."))}

#Just for the CostPath we want log
qplot(data=HypBox[HypBox$variable=="CostPathCost",],x=L1,y=log(value), geom="boxplot") + theme_bw() +ylab("log(CostPathCost)") + geom_hline(aes(yintercept=log(median(data.d[,"CostPathCost"],na.rm=TRUE))),col='red')
ggsave(width=12,height=7,paste(names(Evar["CostPathCost"]),"jpeg",sep="."))

#remove duplicates from cross comparison on everything but Cost Path, since thats not bidirectional
remove.bi<-function(y){
  y.split<-split(y,y$L1)
  tofill<-lapply(y.split,function(z){
    dub<-apply(z,1,function(x){
      paste(min(x["To"],x["From"]),max(x["To"],x["From"]),sep="_")
    })
    out<-z[!duplicated(dub),]
    return(out)})
  return(rbind.fill(tofill))}

AnPrecip<-remove.bi(HypBox[HypBox$variable=="AnnualPrecip",])
AnPrecip$variable<-"Precipitation (mm)"
Elev.<-remove.bi(HypBox[HypBox$variable=="Elev",])
Elev.$variable<-"Elevation (m)"
Euclid.<-remove.bi(HypBox[HypBox$variable=="Euclid",])
Euclid.$variable="Euclid (km)"


library(gridExtra)
Precip<-qplot(data=AnPrecip,x=L1,y=value, geom="boxplot") + theme_bw() + ylab("") + geom_hline(aes(yintercept=median(data.d[,"AnnualPrecip"],na.rm=TRUE)),col='grey40',linetype="dashed") + opts( axis.text.y = theme_blank()) + xlab("") + coord_flip()
ggsave("Precip.jpeg")
Elev<-qplot(data=Elev.,x=L1,y=value, geom="boxplot") + theme_bw() + ylab("") + geom_hline(aes(yintercept=median(data.d[,"Elev"],na.rm=TRUE)), linetype="dashed",col='grey40') + opts( axis.text.y = theme_blank()) + xlab("") + coord_flip()
ggsave("Elev.jpeg")
Euclid<-qplot(data=Euclid.,x=L1,y=value, geom="boxplot") + theme_bw() + ylab("") + geom_hline(aes(yintercept=median(data.d[,"Euclid"],na.rm=TRUE)), linetype="dashed",col='grey40') + opts( axis.text.y = theme_blank()) + xlab("") + coord_flip()
ggsave("Euclid.jpeg",height=3,width=10,dpi=300)
#note outliers
CostC<-qplot(data=HypBox[HypBox$variable=="CostPathCost",],x=L1,y=value, geom="boxplot") + theme_bw() + ylab("") + geom_hline(aes(yintercept=median(data.d[,"CostPathCost"],na.rm=TRUE)), linetype="dashed",col='grey40') + opts( axis.text.y = theme_blank()) + xlab("") + coord_flip() + ylim(0,7.5e8)
ggsave("CostC.jpeg")

#Or try it rbind
#I want to remove those extreme CostPath Values, everyhitng above 1e9
#Only two data points out there, make sure to put in figure legend.
CostC.<-HypBox[HypBox$variable=="CostPathCost",]
CostC.<-CostC.[CostC.$value < 1e9,]
CostC.$variable<-"Cost Path"

#THe names need to match with the entire dataset
levels(HypBox$variable)[8]<-"Cost Path"
levels(HypBox$variable)[7]<-"Euclid (km)"
levels(HypBox$variable)[1]<-"Precipitation (mm)"
levels(HypBox$variable)[2]<-"Elevation (m)"

toplot<-rbind(AnPrecip,Elev.,Euclid.,CostC.)
colnames(toplot)<-c("To","From","env","val","hyp")
toplot$env<-as.factor(toplot$env)
intercepts<-aggregate(toplot$val,list(toplot$env),median,na.rm=T)
colnames(intercepts)<-c("env","int")
ggplot(toplot) + geom_boxplot(aes(hyp,val)) + theme_bw() + theme(axis.text.x=element_blank()) + xlab("") + geom_hline(data=intercepts,aes(yintercept=int,group=env), linetype="dashed",col='grey40') + ylab("") + facet_wrap(~env,scales="free")
ggsave("FacetBoxRowSwap10quant.pdf",dpi=300,height=8,width=12)


#Draw Lines between all hypothesis one sites
setwd("C:\\Users\\Jorge\\Dropbox\\Shared Ben and Catherine\\DimDivEntire\\Files for Analysis")  ###Change this to the Files for Analysis folder in your dropbox, no need to move it. 
elevr<-raster("studyarea_1km.tif")

#Function for spatial lines for all hypothesis, overlayed on an Ecuador Map
setwd("C:\\Users\\Jorge\\Dropbox\\Shared Ben and Catherine\\DimDivEntire\\LineMaps90")

#Lines for each hypothesis, plotted individually
for (f in 1:length(Hyplist.raw)){
  jpeg(paste(names(Hyplist.raw)[f],"RowSwap.jpeg",sep=""),height= 10, width = 10, unit="in", res=300)
  coords<-apply(Hyplist.raw[[f]],1,function(x){
    x_stat<-Envtable[Envtable$CommID %in% as.numeric(x["To"]),c("LatDecDeg","LongDecDeg")]
    y_stat<-Envtable[Envtable$CommID %in% as.numeric(x["From"]),c("LatDecDeg","LongDecDeg")]
    comb<-data.frame(x_stat,y_stat)
    colnames(comb)<-c("xmin","ymin","xmax","ymax")
    return(comb)
  })
  cordmatrix<-rbind.fill(coords)
  plot(elevr, axes=FALSE,main="",legend=F,col=colorRampPalette(c("grey90", "grey6"))(255))
  for (j in 1:nrow(cordmatrix)) {
    arrows(y0=cordmatrix[j,1],x0=cordmatrix[j,2],y1=cordmatrix[j,3],x1=cordmatrix[j,4],length=0, lwd=1,col="grey20")}
  points(loc, col='grey20', pch=15,cex=.5)
  dev.off()
}

Hyplist.raw<-Hyplist.raw[sapply(Hyplist.raw,nrow) > 2]

# Together on one graph, if needed.
jpeg("AllhypothesisRowSwap.jpeg",quality=100,res=300,units="in",height=12,width=20)
par(mfrow=c(2,4))

#Plot lines

for (x in 1:length(Hyplist.raw)){
  coords<-apply(Hyplist.raw[[x]],1,function(x){
    x_stat<-Envtable[Envtable$CommID %in% as.numeric(x["To"]),c("LatDecDeg","LongDecDeg")]
    y_stat<-Envtable[Envtable$CommID %in% as.numeric(x["From"]),c("LatDecDeg","LongDecDeg")]
    comb<-cbind(x_stat,y_stat)
    colnames(comb)<-c("xmin","xmax","ymin","ymax")
    return(comb)
  })
  cordmatrix<-rbind.fill(coords)
  plot(elevr, axes=TRUE, main=names(Hyplist.raw)[x])
  points(loc, col='red', pch=10,cex=.5)
  for (x in 1:nrow(cordmatrix)) {
    arrows(y0=cordmatrix[x,1],x0=cordmatrix[x,2],y1=cordmatrix[x,3],x1=cordmatrix[x,4],length=0, lwd=1)
    
  }
}
dev.off()


#Plot residuals of PCD
head(data.d)

#Function for spatial lines for all hypothesis
setwd("C:\\Users\\Jorge\\Dropbox\\Shared Ben and Catherine\\DimDivEntire\\Output Data")

ggplot(data.d,aes(x=PCDc,y=PCDp,col=Elev)) + geom_point() + theme_bw() + scale_color_gradient("Elevation",low="gray90",high="black") + xlab("Taxonomic") + ylab("Phylogenetic") + ylim(0,1.8) + xlim(0,1.8)
ggsave("CvP.pdf",height=7,width=7,dpi=300)

ggplot(data.d,aes(x=PCDc,y=PCDp.f,col=Elev)) + geom_point() + theme_bw() + scale_color_gradient("Elevation",low="gray90",high="black",guide="none") + xlab("Taxonomic") + ylab("Trait") + ylim(0,1.8) + xlim(0,1.8)
ggsave("CvF.pdf",height=7,width=7,dpi=300)

ggplot(data.d,aes(x=PCDp,y=PCDp.f,col=Elev)) + geom_point() + theme_bw() + scale_color_gradient("Elevation",low="gray90",high="black") + xlab("Phylogenetic") + ylab("Trait") + ylim(0,1.8) + xlim(0,1.8) + geom_abline()
ggsave("PvT.pdf",height=7,width=7.5,dpi=300)

}

###Done