#Post Cluster Figure Creation
#DimDivScript

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


#Set dropbox path
droppath<-"C:/Users/Ben//Dropbox/"

#load data from cluster

load(paste(droppath,"Shared Ben and Catherine/DimDivRevision/500Iterations/DimDivRevisionCluster.RData",sep=""))

data.df<-read.csv(paste(droppath,"Shared Ben and Catherine/DimDivRevision/500Iterations/FinalData.csv",sep=""))

data.df.null<-read.csv(paste(droppath,"Shared Ben and Catherine/DimDivRevision/500Iterations/FinalDataNull.csv",sep=""))


#There is an anomaly in the null cluster, if the assemblages are identical, there is no quantile of the distribution, so the cumulative distribution is =1, thus making the null high
data.df.null[data.df.null$Sorenson==0,]

data.df.null[data.df.null$Sorenson==0,"Phylosor.Phylo_Null"]<-"Low"
data.df.null[data.df.null$Sorenson==0,"MNTD_Null"]<-"Low"

#set to low!


##########################################################################################
#Tables and Statistics
#########################################################################################

setwd(paste(droppath,"Shared Ben and Catherine\\DimDivRevision\\Results",sep=""))

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
  range_prev<-table(data.df.null[,x])/(nrow(comm)*(nrow(comm)-1)/2)
  })

names(data_prev)<-colnames(data.df.null)[15:17]
data_prev<-melt(data_prev)
data_prev<-cast(data_prev,L1~Var.1)
rownames(data_prev)<-data_prev[,1]
data_prev<-data_prev[,-1]

write.csv(round(data_prev,3)*100,"NullPrevalence.csv")

##################################################################
#######################ScatterPlots###############################
##################################################################

#Trait versus Phylogenetic, 
p<-ggplot(data.df,aes(y=MNTD,x=Phylosor.Phylo,col=Sorenson)) + geom_point() 
p<-p+ theme_bw() + scale_color_gradient("Sorenson",low="gray90",high="black")
p<-p+ ylab("Trait Betadiversity") + xlab("Phylogenetic Betadiversity")
p
ggsave("PhylosorPhylovConvexHull_Taxonomic.svg",height=7,width=7.5,dpi=300)

#Phylo phylosor and Hulls
p<-ggplot(data.df,aes(y=MNTD,x=Phylosor.Phylo,col=Elev)) + geom_point() 
p<-p+ theme_bw() + scale_color_gradient("Elevation",low="gray90",high="black")
p<-p+ ylab("Trait Convex Hull") + xlab("Phylogenetic Phylosor") 
p
ggsave("PhylosorPhylovConvexHull_Elevation.svg",height=7,width=7.5,dpi=300)

#Show ranges and plots
#Phylogenetic
ggplot(data.df.null,aes(y=Phylosor.Phylo,fill=Sorenson_Null,x=Phylosor.Phylo_Null)) + geom_boxplot() + theme_bw()
ggsave("PhylogeneticRanges.svg",width=9,height=9)

#Trait
ggplot(data.df.null,aes(y=MNTD,fill=Sorenson_Null,x=MNTD_Null)) + geom_boxplot() + theme_bw()
ggsave("TraitRanges.svg",width=9,height=9)

ggplot(data.df.null,aes(y=MNTD,fill=Phylosor.Phylo_Null,x=MNTD_Null)) + geom_boxplot() + theme_bw() + facet_wrap(~Sorenson_Null)

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
  setwd(paste(droppath,"Shared Ben and Catherine\\DimDivRevision\\Results",sep=""))
  
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
  prop.Hyp<-melt(sapply(Hyplist,nrow)/nrow(data.df)*100)
  prop.Hyp$Hyp<-as.character(rownames( prop.Hyp))
  write.csv(prop.Hyp,"ProportionHypotheses.csv")
  
  ###########################As i see it we need two randomization test, to say the difference in elevation is "high" we should compare that to the dataset as a whole.
  #To say that env differences are between hypothesis, we need to swap those labels.
  #1) Are the environmental differences "high" compared to the dataset?
  #Step 1, get the actual mean values for the env variables for each hypothesis
  
  true.median<-sapply(Hyplist,function(x){
    apply(x[,colnames(x) %in% c("Euclid","CostPathCost","AnnualPrecip","Elev")],2,mean,na.rm=TRUE)})
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
    
    CI<-boot.ci(boot(x$value,function(j,i) mean(j[i],na.rm=T), R=1000)) 
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
    write.csv(Hyplist[[x]],paste(n,".csv"),row.names=FALSE)}
  
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
  HypBox[HypBox$variable=="CostPathCost" & HypBox$value > quantile(data.df$CostPathCost,.95),"value"]<-NA
  
  #Get the medians
  colnames(HypBox)<-c("To","From","Diss","value","Hyp")
  intercepts<-data.frame(Diss=levels(HypBox$Diss),value=round(apply(data.df[,levels(HypBox$Diss)],2,median,na.rm=TRUE),2))
  
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
  elevr<-raster(paste(droppath,"Shared Ben and Catherine\\DimDivEntire\\Files for Analysis\\studyarea_1km.tif",sep=""))
  
  #Function for spatial lines for all hypothesis, overlayed on an Ecuador Map
  
  setwd(paste(droppath,"Shared Ben and Catherine\\DimDivRevision\\Results",sep=""))
  setwd(paste(Tax,Phylo,Func,sep="_"))
  dir.create("Maps")
  setwd("Maps")
  
  #only Create maps of three way comparisons
  
  #Lines for each hypothesis, plotted individually
  #plot all combinations
  map_names<-c("Low.Low.Low","High.Low.Low","High.High.Low","Low.High.Low","High.High.High","Low.High.High","Low.Low.High","High.Low.High")
  
  cl<-makeCluster(7,"SOCK")
  maps.Hyp<-foreach(f=1:length(map_names),.packages=c("reshape","raster","rgdal")) %do% {
    jpeg(paste(map_names[f],".jpeg",sep=""),height= 10, width = 10, unit="in", res=300)
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
    plot(elevr, axes=FALSE,main="",legend=F,col=colorRampPalette(c("grey95", "grey10"))(255))
    if(length(Hyplist[[map_names[f]]])>0){
      for (j in 1:nrow(cordmatrix)) {
        arrows(y0=cordmatrix[j,1],x0=cordmatrix[j,2],y1=cordmatrix[j,3],x1=cordmatrix[j,4],length=0, lwd=1,col="black")}
    }
    points(loc, col='grey10', pch=15,cex=.5)
    dev.off()
  }
  stopCluster(cl)
  
  ############################################################
  #Univaraiate Results
  ############################################################
  
  #Split hypothesis
  SorensonHL<-split(data.df.null,data.df.null$Sorenson_Null)
  PhylosorHL<-split(data.df.null,data.df.null$Phylosor.Phylo_Null)
  MNTDHL<-split(data.df.null,data.df.null$MNTD_Null)
  
  one_way<-c(SorensonHL,PhylosorHL,MNTDHL)
  
  names(one_way)<-paste(rep(c("Taxonomic","Phylogenetic","Trait"),each=3),names(one_way))
  
  #Delinate the combinations of betadiversity
  Hyplist<-one_way[sapply(one_way,nrow) > 0]
  remove.level<-names(Hyplist)[str_detect(names(Hyplist),"Random")]
  Hyplist<-Hyplist[!names(Hyplist) %in% remove.level]
  HypBox<-melt(Hyplist, id.vars=c("To","From"),measure.vars=c("Euclid","CostPathCost","AnnualPrecip","Elev"))
  
  #Proportion of combinations for each null model, later make a recursive object that grabs this output for plotting
  prop.Hyp<-melt(sapply(Hyplist,nrow)/nrow(data.df)*100)
  prop.Hyp$Hyp<-as.character(rownames( prop.Hyp))
  write.csv(prop.Hyp,"UnivariateProportionHypotheses.csv")
  
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
  ggsave("Appendix2oneway.jpeg",dpi=300,height=6,width=8)
  
  #Now you have a data frame with the difference in means for all reps and the true stat at the end
  #Our goal is to figure out if its within the 95 CI interval and then cast it as a matrix
  
  ######################################
  #Figure Creation and Export
  ######################################
  
  for(x in 1:length(Hyplist)){
    n<-names(Hyplist[x])
    write.csv(Hyplist[[x]],paste(n,".csv"),row.names=FALSE)}
  
  #Show three way combinations only.
 
  setwd(paste(droppath,"Shared Ben and Catherine\\DimDivRevision\\Results",sep=""))
  setwd(paste(Tax,Phylo,Func,sep="_"))
  dir.create("1WayBoxplots")
  setwd("1WayBoxplots")
  
  ###############################################
  #To do, remove outlier values from CostPath?
  ###############################################
  
  head(HypBox)
  
  
  colnames(HypBox)<-c("To","From","Diss","value","Hyp")
  
  #split the dimension and the direction out 
  HypBox[,c("Dimension","Direction")]<-colsplit(HypBox$Hyp," ",c("Dimension","Direction"))
  
  #Get the true medians across the entire dataset
  intercepts<-data.frame(Diss=levels(HypBox$Diss),value=round(apply(data.df[,levels(HypBox$Diss)],2,median,na.rm=TRUE),2))
  
  #remove cost path outliers greater than 95th quartile.
  HypBox[HypBox$Diss=="CostPathCost" & HypBox$value > quantile(data.df$CostPathCost,.95),"value"]<-NA

  p<-ggplot(HypBox,aes(x=Direction,y=value,fill=Dimension)) + geom_boxplot()+ facet_wrap(~Diss,scales="free",drop=FALSE) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(data=intercepts,aes(yintercept=value,group=Diss), linetype="dashed",col='grey40') + theme_bw() + scale_x_discrete(drop=FALSE)
  p + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_brewer(palette="Greys")
  
  #remove outliers greater than the 95th quartile
  
  
  ggsave("Env1Boxplots.svg",dpi=300,height=8,width=12)
  ggsave("Env1Boxplots.jpeg",dpi=300,height=8,width=12)
  
  #Draw Lines between all hypothesis one sites
  elevr<-raster(paste(droppath,"Shared Ben and Catherine\\DimDivEntire\\Files for Analysis\\studyarea_1km.tif",sep=""))
  
  #Plot spatial combinations for each High low of tax, phylo, trait.
  for(y in 1:length(one_way)){
    one<-one_way[[y]]
    for (x in 1:length(one)){
      Hyplist<-one[[x]]
      coords<-apply(Hyplist,1,function(x){
        x_stat<-Envtable[Envtable$CommID %in% as.numeric(x["To"]),c("LatDecDeg","LongDecDeg")]
        y_stat<-Envtable[Envtable$CommID %in% as.numeric(x["From"]),c("LatDecDeg","LongDecDeg")]
        comb<-cbind(x_stat,y_stat)
        colnames(comb)<-c("xmin","xmax","ymin","ymax")
        return(comb)
      })
      cordmatrix<-rbind.fill(coords)
      
      #Plot map and create arrows between assemblages
      jpeg(paste(paste(names(one_way)[y],names(one)[x]),".jpeg",sep=""),height= 10, width = 10, unit="in", res=300)
      plot(elevr, axes=FALSE,main=paste(names(one_way)[y],names(one)[x]),legend=F,col=colorRampPalette(c("grey90", "grey6"))(255))
      points(loc, col='red', pch=10,cex=.5)
      for (x in 1:nrow(cordmatrix)) {
        arrows(y0=cordmatrix[x,1],x0=cordmatrix[x,2],y1=cordmatrix[x,3],x1=cordmatrix[x,4],length=0, lwd=1)
      }
      dev.off()
    }
  }
  
  #Read the env analysis for just the univariate results
  
  
  
  ###Done
  #save data from that run
  setwd(paste(droppath,"Shared Ben and Catherine\\DimDivRevision\\Results",sep=""))
  setwd(paste(Tax,Phylo,Func,sep="_"))
  
  save.image("plotting.Rdata")
}

#Run the plotting function for all sets of hypothesis

system.time(Hyplist.func(Tax="Sorenson_Null",Phylo="Phylosor.Phylo_Null",Func="MNTD_Null"))


################################################################
#Reviewer 1 would like to see a better explanation of trait diversity and intra specific variation
#######################################################################

morph<-read.csv(paste(droppath,"Shared Ben and Catherine\\DimDivRevision\\Results\\CompareMetrics\\Morphology.csv",sep=""),na.strings="9999",dec=".")

#just use the species in the analysis
morph<-morph[morph$SpID %in% gsub("_"," ",colnames(siteXspp)),]

#bring in clades

clades<-read.csv(paste(droppath,"Shared Ben and Catherine\\DimDivEntire\\Files for Analysis\\Cladelist.txt",sep=""),header=FALSE)

head(clades)
colnames(clades)<-c("row","Clade","Genus","Species","double","Engish")

morph<-merge(morph,clades[,c("double","Clade")],by.x="SpID",by.y="double")

head(morph)

sd(morph$ExpC,na.rm=TRUE)

m.sp<-aggregate(morph$ExpC,list(morph$SpID),function(x){
  sd(x,na.rm=TRUE)})

mean(m.sp$x,na.rm=TRUE)

#Start with bill length
ggplot(data=morph,aes(SpID,ExpC)) + theme_bw() + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_wrap(~Clade,nrow=3,scales="free_x") + ylab("Bill Length(mm)")

anova(morph$SpID,morph$ExpC)

ggsave(paste(droppath,"Shared Ben and Catherine\\DimDivRevision\\Results\\CompareMetrics\\BillLength.jpeg",sep=""),height=12,width=20,dpi=300,units="in")   
      
#Eliminate any rounding errors for the weight, a couple bad sites from Aves Data Base
#The heaviest hummingbird in the world is around 22g, eliminate anything larger, clearly some data cleaning issue

morph<-morph[morph$Peso < 24 & !is.na(morph$Peso),]

ggplot(data=morph,aes(SpID,Peso)) + theme_bw() + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_wrap(~Clade,nrow=3,scales="free") + ylab("Bill Length(mm)")
ggsave(paste(droppath,"Shared Ben and Catherine\\DimDivRevision\\Results\\CompareMetrics\\BillLength.jpeg",sep=""),height=12,width=20,dpi=300,units="in")   

<<<<<<< HEAD
#############################################
#Multiple Regression Coefficients
###############################################

require(ppcor)
pcor.test(data.df$MNTD,data.df$Phylosor.Phylo,data.df$Sorenson)




=======

#
#################################################################
#Inspect individual points, something is wrong with low low low
#################################################################

#Plot points
plot(loc,cex=.1)
text(loc,labels=loc$CommID,cex=.5)

#review specific comparisons
insp<-function(a,b){
  e<-c(a,b)
  co<-data.df.null[data.df.null$To %in% e & data.df.null$From %in% e,]
  co.sp<-loc[loc$CommID %in% e,]
  plot(elev.test)
  points(co.sp,col="red",pch=16)
  a<-Getsplist(e[[1]])
  b<-Getsplist(e[[2]])
  shared<-sum(a%in%b)
  
  richnessA<-length(a)
  richnessB<-length(b)
  unsharedA<-sum(!a%in%b)
  unsharedB<-sum(!b%in%a)
  
  print(data.frame(shared,richnessA,richnessB,unsharedA,unsharedB))
  
  print(paste(co$Sorenson_Null,co$Phylosor.Phylo_Null,co$MNTD_Null))
  return(co)
  
  }

#WHatkind do fit into low low low

lll<-data.df.null[data.df.null$Sorenson_Null %in% "Low" & data.df.null$Phylosor.Phylo_Null %in% "Low" & data.df.null$MNTD_Null %in% "Low",]

#For an example of a suprising inclusion in lll
#see  129 223

ob<-insp("251","301")

#what is the mean Null distribution for this comparison
e<-c("251","301")

Nulls<-Null_dataframe[Null_dataframe$From %in% e & Null_dataframe$To %in% e,]

head(Nulls)

ecdf(Nulls$Phylosor.Phylo) (ob$Phylosor.Phylo)

ggplot(Nulls,aes(x=Phylosor.Phylo)) + geom_histogram() + geom_vline(aes(xintercept=ob$Phylosor.Phylo),col="Red") + theme_bw()

a<-Getsplist("301")
b<-Getsplist("251")

b[!b %in% a]
a[!a %in% b]


#What is the mean null phyloenetic betadiversiy
aggN<-aggregate(Null_dataframe$Phylosor.Phylo,list(Null_dataframe$To,Null_dataframe$From),mean,na.rm=TRUE)

colnames(aggN)<-c("To","From","meanNull")

#What is the number of shared species
aggN$unshared<-apply(aggN,1,function(x){
  a<-Getsplist(as.character(x[[1]]))
  b<-Getsplist(as.character(x[[2]]))
  unshared<-sum(!b%in%a) + sum(!a%in%b)
})

aggN$totalN<-apply(aggN,1,function(x){
  a<-Getsplist(as.character(x[[1]]))
  b<-Getsplist(as.character(x[[2]]))
  return(length(a)+length(b))
})

#plot ratio
ggplot(aggN,aes(x=unshared/totalN,y=meanNull)) + geom_point() + geom_smooth(method="lm")

# WHen species list are more similiar it is harder to get low phylogenetic diversity

ggplot(data.df.null,aes(y=Sorenson,x=Phylosor.Phylo_Null)) + geom_boxplot()

#Simulation test

doN<-function(richness_To,richness_From,overlapping){
  
  #Create slots for each community
  nullTo<-matrix(nrow=richness_To,ncol=1)
  nullFrom<-matrix(nrow=richness_From,ncol=1)
  
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
  print(vegdist(null_siteXspp,binary=TRUE))
  Full=1-as.numeric(phylosor(null_siteXspp,tree))
  return(Full)}


T_20204<-replicate(500,doN(richness_To=20,richness_From=20,overlapping=4))

T_551<-replicate(500,doN(richness_To=5,richness_From=5,overlapping=1))

p<-ggplot() + geom_density(aes(T_551),col="Red") + geom_density(aes(T_20204),col="blue",add=TRUE,alpha=.8) + xlab("PhyloBeta") + theme_bw()
p<-p + geom_vline(aes(xintercept=quantile(T_20204,.95)),col="Blue",linetype="dashed")
p<-p+ geom_vline(aes(xintercept=quantile(T_551,.95)),col="Red",linetype="dashed")
p<-p + geom_vline(aes(xintercept=quantile(T_20204,.05)),col="Blue",linetype="dashed")
p<-p+ geom_vline(aes(xintercept=quantile(T_551,.05)),col="Red",linetype="dashed")
p + ggtitle("Sorenson=.8,Red<-n=5 shared = 1,Blue<-n=20,shared=4")
ggsave(paste(droppath,"Shared Ben and Catherine\\DimDivRevision\\Results\\TaxNullTest.svg",sep=""),height=9,width=8,dpi=300)

quantile(T_20204,.95)
quantile(T_551,.95)

#Compute phylogenetic and trait metrics on null assemblages

#give a brief biogeographic viewpoint
data.df.null$meanE<-apply(data.df.null,1,function(x){
  a<-Envtable[Envtable$CommID %in% as.numeric(x["To"]),"Elev"]
  b<-Envtable[Envtable$CommID %in% as.numeric(x["From"]),"Elev"]
  return(mean(a,b))
})

ggplot(data.df.null,aes(y=meanE,x=Phylosor.Phylo_Null)) + geom_boxplot()

ggplot(data.df.null,aes(y=meanE,x=Sorenson)) + geom_point() + geom_smooth()

aggN$meanE<-apply(aggN,1,function(x){
  a<-Envtable[Envtable$CommID %in% as.numeric(x["To"]),"Elev"]
  b<-Envtable[Envtable$CommID %in% as.numeric(x["From"]),"Elev"]
  return(mean(a,b))
})

ggplot(aggN,aes(x=meanE,y=shared)) + geom_point() + geom_smooth(method="lm")

ggplot(data.df.null,aes(x=Phylosor.Phylo_Null,y=Phylosor.Phylo)) + geom_boxplot()

data.df.null[data.df.null$Phylosor.Phylo_Null %in% "High",][which.min(data.df.null[data.df.null$Phylosor.Phylo_Null %in% "High",]$Phylosor.Phylo),]
>>>>>>> 6bbdc2e926be5a888630ceb6f34364b6c511b5c7
