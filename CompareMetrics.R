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
sink("C:/Users/Jorge/Dropbox/Shared Ben and Catherine/DimDivRevision/Results/OvernightOutput.txt")
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



##################################
#Define Function to Compare Metrics
###################################
set.seed(100)
comm<-siteXspp[sample(1:nrow(siteXspp)),]


beta_all<-function(comm,tree,traits,FullMatrix){

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
  pair.w<-expand.grid(rownames(comm)[1:length(richness_levels)],rownames(comm)[(length(richness_levels)+1):(length(richness_levels)*2)])
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

prc_traits<-prcomp(mon_cut)
newSGdist <- dist(prc_traits$x)
source("C:/Users/Jorge/Documents/DimDiv/BenHolttraitDiversity.R")

#create sp.list
sp.list<-apply(siteXspp_traits,1,function(x){
  names(x[which(x==1)])
})

dists <- as.matrix(newSGdist)

rownames(dists) <- rownames(mon_cut)
colnames(dists) <- rownames(mon_cut)

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
require(reshape2)

Allmetrics1<-dcast(Allmetrics0,...~PCD.phylo,value.var="PCD.value.phylo")
colnames(Allmetrics1)[15:17]<-c("PCD.phylo","PCDc.phylo","PCDp.phylo")

Allmetrics2<-dcast(Allmetrics1,...~PCD.func,value.var="PCD.value.func")
colnames(Allmetrics2)[16:18]<-c("PCD.func","PCDc.func","PCDp.func")

return(Allmetrics2)}

system.time(beta_metrics<-beta_all(comm=comm,tree=tree,traits=mon,FullMatrix=FALSE))
                            
#Visualizations of the beta metrics
head(beta_metrics)

save.image("C:/Users/Jorge/Dropbox/Shared Ben and Catherine/DimDivEntire/Output Data/Workspace.RData")

##########################################################################################
#Tables and Statistics
#########################################################################################
data.df<-beta_metrics

setwd("C:\\Users\\Jorge\\Dropbox\\Shared Ben and Catherine\\DimDivRevision\\Results")

#Get the bounds of each 
range_metrics<-list()

for(x in 0:15){
  print(x)
range_min<-aggregate(data.df[,12+x],list(data.df[,28+x]),min,na.rm=TRUE)
range_max<-aggregate(data.df[,12+x],list(data.df[,28+x]),max,na.rm=TRUE)

range_val<-data.frame(Index=range_min[,1],Min=range_min[,2],Max=range_max[,2])
range_metrics[[x+1]]<-range_val
}

names(range_metrics)<-colnames(data.df)[12:27]
range_metrics<-melt(range_metrics)

write.csv(range_metrics,"Range_Metrics.csv")

#Find Prevalence of each combination
data_prev<-lapply(colnames(data.df)[28:42],function(x){
range_prev<-table(data.df[,x])/nrow(data.df)})

names(data_prev)<-colnames(data.df)[28:42]
data_prev<-melt(data_prev)
data_prev<-cast(data_prev,L1~Var.1)
rownames(data_prev)<-data_prev[,1]
data_prev<-data_prev[,-1]

write.csv(round(data_prev,3)*100,"NullPrevalence.csv")

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

##########################################
#Plot each of the metrics with their ranges
##########################################
range_plots<-lapply(12:26,function(x){
print(colnames(data.df)[x])
print(colnames(data.df)[x+15])
p<-ggplot(data.df,aes(y=data.df[,colnames(data.df)[x]],x=data.df[,colnames(data.df)[x+15]])) + geom_boxplot()
p<-p+labs(y=colnames(data.df)[x],x=colnames(data.df)[x+15])
filnam<-paste(colnames(data.df)[x],"_range.jpeg")
ggsave(filnam,plot=p,height=7,width=5,dpi=300)
return(p)})


###############################
#Correlation
###############################
#Just get the columns we want to run in the metrics
data.d<-data.df[,colnames(data.df) %in% c("beta_functional","Phylosor.Phylo","Phylosor.Func","Sorenson","PCDp.phylo","PCDp.func")]

#Get the colums we want to run in env
#Set 1: Get the correlations from the dataset 
data.env<-data.df[,colnames(data.df) %in% c("AnnualPrecip","Elev", "CostPathCost","DeltaElevEuclid","H_mean","Tree","Euclid")]
true_cor<-sapply(colnames(data.d), function(x){ sapply(colnames(data.env),function(y){
  cor(data.d[,x],data.env[,y],method="spearman",use="complete.obs")
})})

write.csv(true_cor,"Env_Cor.csv")


#################################
#Compare Metrics
#################################
require(stringr)
#
require(reshape)
setwd("C:\\Users/Jorge/Dropbox/Shared Ben and Catherine/DimDivRevision/Results/")
#get the prevalances for each hypothesis for each species

Hyp.all<-list.files(full.names=TRUE,pattern="ProportionHypotheses",recursive=TRUE)

Hyp.dat<-lapply(Hyp.all,read.csv)

#I'm no good at regex
names(Hyp.dat)<-sapply(Hyp.all,function(x){
  strsplit(x,"/")[[1]][2]
})

m.dat<-melt(Hyp.dat,variable_name="ignore")
colnames(m.dat)[4]<-"Prevalence"

#Split into component pieces
m.dat<-data.frame(m.dat,colsplit(m.dat$Hyp,"\\.",c("Taxonomic","Phylogenetic","Trait")))
#melt those pieces
m.dat<-melt(m.dat,measure.vars=c("Taxonomic","Phylogenetic","Trait"),variable_name="Betadiversity")
colnames(m.dat)[7]<-"Combination"

#drop the NULL in the word labels
m.dat$L1<-gsub("_Null","",m.dat$L1)

#all possible combinations
ggplot(m.dat,aes(fill=L1,y=Prevalence,x=Hyp)) +geom_bar(position="dodge") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_wrap(~Hyp,scales="free_x") + scale_x_discrete(labels="") + labs(fill="Metrics") + theme_bw()
ggsave("Metric_compare.jpeg",height=8,width=12)

#Without random combinations
m.datNoRandom<-m.dat[!m.dat$Hyp %in% levels(m.dat$Hyp)[str_detect(levels(m.dat$Hyp),"Random")],]

ggplot(m.datNoRandom,aes(fill=L1,y=Prevalence,x=Hyp)) +geom_bar(position="dodge") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_wrap(~Hyp,scales="free_x") + scale_x_discrete(labels="") + labs(fill="Metrics") + theme_bw()
ggsave("Metric_compareNorandom.jpeg",height=8,width=11)

#plot dimensions of betadiversity
dim_beta<-aggregate(m.dat$Prevalence,list(m.dat$Combination,m.dat$Betadiversity,m.dat$L1),sum)
colnames(dim_beta)<-c("Combination","Betadiversity","Metrics","Prevalence")
ggplot(dim_beta,aes(fill=Metrics,y=Prevalence,x=Combination)) + geom_bar(position="dodge")+ facet_wrap(~Betadiversity
) + theme_bw() + labs(fill="Metrics")
ggsave("MetricComponents.jpeg",height=8,width=11)

#Another way to look at components
#Taxonomic
TaxC<-melt(data.df,id.var=c("To","From"),measure.vars=c("Sorenson_Null","PCDc.phylo_Null"))
p<-ggplot(TaxC,aes(x=value,fill=variable)) + geom_bar(col="black",position="dodge",aes(y=(..count..)/23871)) + labs(fill="Metric") + scale_y_continuous(label = percent) + theme_bw() + scale_fill_brewer(palette="Greys")
p + scale_color_continuous(guide='none') + ylab("Prevalence") + xlab("")
ggsave("Taxmetrics.svg",dpi=300,height=8,width=11)

PhyloC<-melt(data.df,id.var=c("To","From"),measure.vars=c("PCDp.phylo_Null","Phylosor.Phylo_Null"))
p<-ggplot(na.omit(PhyloC),aes(x=value,fill=variable)) + geom_bar(col="black",position="dodge",aes(y=(..count..)/23871)) + labs(fill="Metric") + scale_y_continuous(label = percent) + theme_bw() + scale_fill_brewer(palette="Greys")
p + scale_color_continuous(guide='none') + ylab("Prevalence") + xlab("")
ggsave("Phylometrics.svg",dpi=300,height=8,width=11)

TraitC<-melt(data.df,id.var=c("To","From"),measure.vars=c("PCDp.func_Null","Phylosor.Func_Null","beta_functional_Null",'MNTD_Null'))
p<-ggplot(na.omit(TraitC),aes(x=value,fill=variable)) + geom_bar(col="black",position="dodge",aes(y=(..count..)/23871)) + labs(fill="Metric") + scale_y_continuous(label = percent) + theme_bw() + scale_fill_brewer(palette="Greys")
p + scale_color_continuous(guide='none') + ylab("Prevalence") + xlab("")
ggsave("Traitmetrics.svg",dpi=300,height=8,width=11)

save.image("C:/Users/Jorge/Dropbox/Shared Ben and Catherine/DimDivRevision/Results/DimDivRevision.RData")

#prevalence of any combination
round(ftable(data.df$Sorenson_Null,data.df$Phylosor.Phylo_Null,data.df$beta_functional_Null)/23871,digits=3)*100

#detach our overnight sink
sink()
