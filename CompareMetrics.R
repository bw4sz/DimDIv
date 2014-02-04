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

##########################################################
#Read in data
##########################################################


#Set dropbox path
droppath<-"C:/Users/Jorge/Dropbox/"

#Set git path
gitpath<-"C:/Users/Jorge/Documents/DimDiv/"


#load data if desired
load(paste(droppath,"Shared Ben and Catherine/DimDivRevision/Results/CompareMetrics.RData",sep=""))

###Define Source Functions

source(paste(gitpath,"geb12021-sup-0004-si.R.txt",sep=""))
#Function is called beta tf

##########################################################
#Read in data
##########################################################

#Load in the data
#load("C:/Users/Ben/Dropbox/Shared Ben and Catherine/DimDivEntire/Output Data/Workspace.RData")

#Read in species matrix
siteXspp <- read.csv(paste(gitpath,"InputData/siteXspp_Oct20_2011.csv", sep=""), row.names=1)

#Get entire species list
splist<-colnames(siteXspp)

#Read in phylogeny
tree<-read.nexus(paste(gitpath,"InputData/ColombiaPhylogenyUM.tre",sep=""))

#Read in names file to replace names in Nexis file
spnames<-read.table(file=paste(gitpath,"InputData/SpNameTree.txt",sep="") , sep = "\t", header = TRUE)

#Replace tip.label with Spnames#
tree$tip.label<-as.character(spnames$SpName) 


#bring in traits
morph <- read.csv(paste(gitpath,"InputData/MorphologyShort.csv",sep=""),na.strings="9999")

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

#Zscores, standardized by sd and subtracted means
means<-apply(mon,2,mean)

Bill<-mon$Bill - means["Bill"]/sd(mon$Bill)
Mass<-mon$Mass - means["Mass"]/sd(mon$Mass)
WingChord<-(mon$WingChord - means["WingChord"])/sd(mon$WingChord)

z.scores<-data.frame(Bill,Mass,WingChord)
rownames(z.scores)<-rownames(mon)

#Functional dendrogram
tree.func<-as.phylo(hclust(dist(z.scores)))

clades<-read.csv(paste(droppath,"Shared Ben and Catherine\\DimDivEntire\\Files for Analysis\\Cladelist.txt",sep=""),header=FALSE)
head(clades)
colnames(clades)<-c("row","Clade","Genus","Species","double","Engish")

head(clades)
clades$Dash<-gsub(" ","_",clades$double)
#Color the func phylo by clade # this only works the 2nd time through
col.clade<-as.vector(sapply(tree.func$tip.label,function(x) clades[clades$Dash==x,"Clade"]))
plot.phylo(tree.func, cex=.35, tip.color=as.numeric(as.factor(col.clade)))

##################################
#Define Function to Compare Metrics
###################################
comm<-siteXspp


beta_all<-function(comm,tree,traits,cores){

#Phylogenetic PCD
PCD.phylo<-pcd(comm,tree)

#turn to matrix, and set the diagonal to NA
PCD.phylo<-melt(
  lapply(PCD.phylo,function(x){
  m<-as.matrix(x)
  diag(m)<-NA
  return(m)})
  )

colnames(PCD.phylo)<-c("To","From","PCD.value.phylo","PCD.phylo")

#Phylosor Calculation see Bryant 2008
phylo.matrix<-as.matrix(phylosor(comm,tree))
diag(phylo.matrix)<-NA
Phylosor.phylo<-melt(phylo.matrix)
colnames(Phylosor.phylo)<-c("To","From","Phylosor.Phylo")
Phylosor.phylo$Phylosor.Phylo<-1-Phylosor.phylo$Phylosor.Phylo

#Merge phylometrics
phylometrics<-merge(Phylosor.phylo,PCD.phylo,c("To","From"))

#####################################
##Taxonomic Betadiversity
sorenson<-melt(as.matrix(vegdist(comm,binary=TRUE)))
colnames(sorenson)<-c("To","From","Sorenson")
######################################

#Functional dendorgram PCD
PCD.func<-pcd(comm,tree=tree.func)
PCD.func<-melt(
  lapply(PCD.func,function(x){
    m<-as.matrix(x)
    diag(m)<-NA
    return(m)})
)

colnames(PCD.func)<-c("To","From","PCD.value.func","PCD.func")

#Phylosor Calculation see Bryant 2008
Phylosor.func<-melt(as.matrix(phylosor(comm,tree.func)))
colnames(Phylosor.func)<-c("To","From","Phylosor.Func")
Phylosor.func$Phylosor.Func<-1-Phylosor.func$Phylosor.Func

#Merge phylometrics
funcmetrics<-merge(Phylosor.func,PCD.func,c("To","From"))

#merge this into the phylometrics
Phylo_Tax2<-merge(phylometrics,funcmetrics,c("To","From"))

#################
#Trait Metrics
#################
Phylo_Tax<-merge(Phylo_Tax2,sorenson,c("To","From"))

Phylo_Tax3<-dcast(Phylo_Tax,...~PCD.phylo,value.var="PCD.value.phylo")
Phylo_Tax4<-dcast(Phylo_Tax3,...~PCD.func,value.var="PCD.value.func")

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
cl<-makeCluster(cores,"SOCK")
registerDoSNOW(cl)
system.time(pairwise.beta<-foreach(x=pair.w,.packages=c("vegan","reshape"),.errorhandling="pass") %dopar%{
  source(paste(gitpath,"geb12021-sup-0004-si.R.txt",sep=""))
  Villeger<-beta_TF(siteXspp_traits[rownames(siteXspp_traits) %in% x,] ,as.matrix(mon_cut))$beta
  Villeger<-melt(Villeger)
  cast(Villeger,~X1+X2)
  })
stopCluster(cl)

toremove<-sapply(pairwise.beta,function(x) is.character(x[[1]]))

#get rid of the NA rows
toremove<-sapply(pairwise.beta,function(x) is.character(x[[1]]))
pairwise.beta.removed<-rbind.fill(pairwise.beta[!toremove])

#Get the order of inputs
pairwise.order<-t(sapply(pair.w,function(x) {matrix(nrow=1,ncol=2,x)}))
pairwise.order.removed<-pairwise.order[!toremove,]
func.beta<-data.frame(pairwise.order.removed,pairwise.beta.removed)[,-3]

#Combine the dataframes
colnames(func.beta)[1:2]<-c("To","From")
Allmetrics<-merge(func.beta,Phylo_Tax,by=c("To","From"))


#Trait distance, sum of the trait distance, first use principal components and directly uses euclidian space
#Courtesy of ben holt, imperial college

## siteXspp_traits = community matrix (grid cell = rows, species = cols)
## mon_cut = trait data
## sp.list is a list with the names of the species found in each grid cell
# the name of each entry in sp.list is the name of the grid cell it refers to


source(paste(gitpath,"BenHolttraitDiversity.R",sep=""))

#create sp.list
sp.list<-apply(siteXspp_traits,1,function(x){
  names(x[which(x==1)])
})

dists<-as.matrix(dist(traits,method="euclidean"))

rownames(dists) <- rownames(traits)
colnames(dists) <- rownames(traits)

sgtraitMNTD <- sapply(rownames(siteXspp_traits),function(i){
  
  #Iterator count
  print(round(which(rownames(siteXspp_traits)==i)/nrow(siteXspp_traits),3))
  
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
Allmetrics0<-merge(Allmetrics,melt.MNTD,by=c("To","From"))

#Cast out the full frame
Allmetrics1<-dcast(Allmetrics0,...~PCD.phylo,value.var="PCD.value.phylo")
colnames(Allmetrics1)[15:17]<-c("PCD.phylo","PCDc.phylo","PCDp.phylo")

Allmetrics2<-dcast(Allmetrics1,...~PCD.func,value.var="PCD.value.func")
colnames(Allmetrics2)[16:18]<-c("PCD.func","PCDc.func","PCDp.func")

return(Allmetrics2)}

system.time(beta_metrics<-beta_all(comm=comm,tree=tree,traits=z.scores,cores=8))
                            
#Visualizations of the beta metrics
head(beta_metrics)

save.image(paste(droppath,"Shared Ben and Catherine/DimDivRevision/Results/CompareMetrics.RData",sep=""))

##########################################################################################
#Tables and Statistics
#########################################################################################
data.df<-beta_metrics

setwd("C:\\Users\\Jorge\\Dropbox\\Shared Ben and Catherine\\DimDivRevision\\Results")

dir.create("CompareMetrics")
setwd("CompareMetrics")

###################################################
#correlation and comparisons

#Taxonomic
svg("Tax_plot.svg")
Tax_plot<-ggpairs(data.df[,c("beta_taxonomic","PCDc.phylo","Sorenson")]) 
Tax_plot
dev.off()

#Phylogenetic
svg("Phylo_plot.svg")
Phylo_plot<-ggpairs(data.df[is.finite(data.df$PCDp.phylo),][,c("Phylosor.Phylo","PCDp.phylo")]) 
Phylo_plot
dev.off()

#Functional
svg("Func_plot.svg",height=8,width=8)
Func_plot<-ggpairs(data.df[is.finite(data.df$PCDp.func),][,c("Phylosor.Func","PCDp.func","beta_functional","MNTD")]) 
Func_plot
dev.off()

###Other Scatter plots of interest
#Function for spatial lines for all hypothesis

#PCD phylo and PCD Taxonomic
ggplot(data.df,aes(x=PCDc.phylo,y=PCDp.phylo,col=Elev)) + geom_point() + theme_bw() + scale_color_gradient("Elevation",low="gray90",high="black") + xlab("Taxonomic PCDc") + ylab("Phylogenetic PCDp") + xlim(0,1.8) + ylim(0,1.8) + theme(aspect.ratio=1)
ggsave("PCDcvPCDpPhylo.svg",height=7,width=7,dpi=300)

#PCD func and PCD Taxonomic
ggplot(data.df,aes(x=PCDc.func,y=PCDp.func,col=Elev)) + geom_point() + theme_bw() + scale_color_gradient("Elevation",low="gray90",high="black",guide="none") + xlab("Taxonomic PCDc") + ylab("Trait PCDp") + xlim(0,1.8) + ylim(0,1.8) + theme(aspect.ratio=1)
ggsave("PCDcvPCDpFunc.svg",height=7,width=7,dpi=300)

#PCD func and PCD phylo
p<-ggplot(data.df,aes(x=PCDp.phylo,y=PCDp.func,col=Elev)) + geom_point() 
p<-p+ theme_bw() + scale_color_gradient("Elevation",low="gray90",high="black")
p<-p+ xlab("Phylogenetic PCDp") + ylab("Trait PCDp") + xlim(0,1.8) + ylim(0,1.8) + theme(aspect.ratio=1)
p
ggsave("PCDpPhylovPCDpFunc.svg",height=7,width=7.5,dpi=300)

#Phylosor v Phylosor
#PCD func and PCD phylo
p<-ggplot(data.df,aes(y=Phylosor.Func,x=Phylosor.Phylo,col=Elev)) + geom_point() 
p<-p+ theme_bw() + scale_color_gradient("Elevation",low="gray90",high="black")
p<-p+ ylab("Phylogenetic Func") + xlab("Phylogenetic Phylosor") + coord_equal()
p
ggsave("Phylosor_Elevation.svg",height=7,width=7.5,dpi=300)

#Phylosor v Phylosor
#PCD func and PCD phylo
p<-ggplot(data.df,aes(y=Phylosor.Func,x=Phylosor.Phylo,col=Sorenson)) + geom_point() 
p<-p+ theme_bw() + scale_color_gradient("Sorenson",low="gray90",high="black")
p<-p+ ylab("Taxonomic Sorenson") + xlab("Phylogenetic Phylosor") + coord_equal()
p
ggsave("Phylosor_Taxonomic.svg",height=7,width=7.5,dpi=300)

#MNTD and hulls
p<-ggplot(data.df,aes(y=beta_functional,x=MNTD,col=Sorenson)) + geom_point() 
p<-p+ theme_bw() + scale_color_gradient("Sorenson",low="gray90",high="black")
p<-p+ ylab("Trait Convex Hull") + xlab("MNTD") + coord_equal()
p
ggsave("MNTDvConvexHull_Taxonomic.svg",height=7,width=7.5,dpi=300)

#MNTD and hulls
p<-ggplot(data.df,aes(y=beta_functional,x=MNTD,col=Elev)) + geom_point() 
p<-p+ theme_bw() + scale_color_gradient("Elevation",low="gray90",high="black")
p<-p+ ylab("Trait Convex Hull") + xlab("MNTD") + coord_equal()
p
ggsave("MNTDvConvexHull_Elevation.svg",height=7,width=7.5,dpi=300)

#
#PCDp Phylo and hulls
#PCD func and PCD phylo
p<-ggplot(data.df,aes(y=beta_functional,x=Phylosor.Phylo,col=Sorenson)) + geom_point() 
p<-p+ theme_bw() + scale_color_gradient("Sorenson",low="gray90",high="black")
p<-p+ ylab("Trait Convex Hull") + xlab("Phylogenetic Phylosor") + coord_equal()
p
ggsave("PhylosorPhylovConvexHull_Taxonomic.svg",height=7,width=7.5,dpi=300)

#Phylo phylosor and Hulls
p<-ggplot(data.df,aes(y=beta_functional,x=Phylosor.Phylo,col=Elev)) + geom_point() 
p<-p+ theme_bw() + scale_color_gradient("Elevation",low="gray90",high="black")
p<-p+ ylab("Trait Convex Hull") + xlab("Phylogenetic Phylosor") + coord_equal()
p
ggsave("PhylosorPhylovConvexHull_Elevation.svg",height=7,width=7.5,dpi=300)

#Phylo PCDp and Hulls
p<-ggplot(data.df,aes(y=beta_functional,x=PCDp.phylo,col=Elev)) + geom_point() 
p<-p+ theme_bw() + scale_color_gradient("Sorenson",low="gray90",high="black")
p<-p+ ylab("Trait Convex Hull") + xlab("Phylogenetic PCDp") + coord_equal()
p
ggsave("PhylosorPhylovConvexHull_Elevation.svg",height=7,width=7.5,dpi=300)

#Phylo phylosor and Hulls
p<-ggplot(data.df,aes(y=Phylosor.Phylo,x=PCDc.phylo,col=Elev)) + geom_point() 
p<-p+ theme_bw() + scale_color_gradient("Elevation",low="gray90",high="black")
p<-p+ ylab("Phylosor") + xlab("PCDc") + coord_equal()
p

#Phylo PCDp and Hulls
p<-ggplot(data.df,aes(y=beta_functional,x=PCDp.phylo,col=Sorenson)) + geom_point() 
p<-p+ theme_bw() + scale_color_gradient("Sorenson",low="gray90",high="black")
p<-p+ ylab("Trait Convex Hull") + xlab("Phylogenetic PCDp") + coord_equal()
p
ggsave("PhylosorPhylovConvexHull_Taxonomic.svg",height=7,width=7.5,dpi=300)




###############################
#Correlation
###############################
#Just get the columns we want to run in the metrics
data.d<-data.df[,colnames(data.df) %in% c("beta_functional","Phylosor.Phylo","Phylosor.Func","Sorenson","PCDp.phylo","PCDp.func","MNTD")]

metric_cor<-sapply(colnames(data.d), function(x){ sapply(colnames(data.d),function(y){
  cor(data.d[,x],data.df[,y],method="spearman",use="complete.obs")
})})

write.csv(metric_cor,"Metric_Cor.csv")

###########Compare metrics of trait 
#Dendrogram
dendogram.trait<-cophenetic(tree.func)
prc_traits<-prcomp(mon)

#PCA dist
newSGdist <- as.matrix(dist(prc_traits$x))

#Euclid Dist
euclid.trait<-as.matrix(dist(mon))

#Combine
m.trait<-melt(list(dendrogram=dendogram.trait,MNTD=newSGdist,Euclid=euclid.trait))
colnames(m.trait)<-c("To","From","dist","Metric")

trait.compare<-dcast(m.trait,...~Metric,value.var="dist")

svg("Trait.comparison.svg",height=9,width=8)
ggplot(data=trait.compare[complete.cases(trait.compare),-c(1,2,4)],aes(x=dendrogram,y=MNTD)) + geom_point() + theme_bw()
dev.off()

###########Catherine NSF Grant#####################
data.df2<-dcast(data.df,...~PCD.phylo,value.var="PCD.value.phylo")
colnames(data.df2)[c(15,16,17)]<-paste(colnames(data.df2)[c(15,16,17)],"Phylo",sep="_")
data.df2<-dcast(data.df2,...~PCD.func,value.var="PCD.value.func")
colnames(data.df2)[c(15,16,17)]<-paste(colnames(data.df2)[c(15,16,17)],"Phylo",sep="_")

ggplot(data.df2,aes(x=PCD_Phylo,y=PCDp)) + geom_point() + xlim(c(0,2)) + ylim(c(0,2)) + 

ggplot(data.df,aes(x=))
