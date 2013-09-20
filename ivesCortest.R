#Ives 2010 correlation for Randomization and PCD method
#in #### are notes directly from ives

setwd("C:\\Users\\Jorge\\Dropbox\\Shared Ben and Catherine\\DimDivRevision\\Results")

#Written by Ben Weinstein Stony Brook University

#Read in data

data.df<-read.csv("FinalData.csv")[,-1]

#Just get the columns we want to run in the metrics
data.d<-data.df[,colnames(data.df) %in% c("beta_functional","Phylosor.Phylo","Phylosor.Func","Sorenson","PCDp.phylo","PCDp.func")]

#Get the colums we want to run in env
#Set 1: Get the correlations from the dataset 
data.env<-data.df[,colnames(data.df) %in% c("AnnualPrecip","Elev", "CostPathCost","DeltaElevEuclid","H_mean","Tree","Euclid")]
true_cor<-sapply(colnames(data.d), function(x){ sapply(colnames(data.env),function(y){
  cor(data.d[,x],data.env[,y],method="spearman",use="complete.obs")
})})

write.csv(true_cor,"Env_Cor.csv")

metric_cor<-sapply(colnames(data.d), function(x){ sapply(colnames(data.d),function(y){
  cor(data.d[,x],data.df[,y],method="spearman",use="complete.obs")
})})

write.csv(metric_cor,"Metric_Cor.csv")

true_cor<-abs(true_cor)


#### Ives 2010: To do this, we first created 10,000 data sets by randomly permutingmenvironmental variables among communities and correlating these to the observed PCD values.
#For all environmental variables
cl<-makeCluster(8,"SOCK")
registerDoSNOW(cl)

  full<-list()
  for (g in 7:2){
  print(g)  
  #get the highest correlation from 10000 sampled datasets.
 
  ives_cor<-times(1000) %dopar% {
    
    topick<-sample(colnames(data.env),g)
    #Resample desired columns
    dat.test<-apply(data.env[,colnames(data.env)[colnames(data.env) %in% topick]],2,function(x)sample(x))
    
    #get correlation between all env variables and PCD values
    cor_all<-sapply(colnames(data.d), function(x){ sapply(colnames(dat.test),function(y){
      cor(data.d[,x],dat.test[,y],method="spearman",use="complete.obs")
    })
    })
    
    #for each column, find the max value
    max_col<-rbind(apply(cor_all,2,max))
    return(list(max_col))}
  
  m.ives_cor<-melt(ives_cor)
  #Create null distributions, each column below is the expected distribution is there no relationship among variables
  observed_cor<-cast(m.ives_cor,L1~X2)[,-1]
  
  #Get the cumulative distribution and the pvalue for every combination of metrics and env variable
  #For the highest correlation, sorted value is the last position [27], the 2nd highest correlation is 2nd to last [26] etc.
  
  cor_rank<-lapply(colnames(data.d),function(x){ 
    final_p<-1-ecdf(observed_cor[,x]) (true_cor[which(rank(true_cor[,x])==g),x])
    m.cor<-melt(true_cor[,x])
    m.cor$rank<-rank(m.cor[,1])
    m.cor<-m.cor[order(m.cor$rank),]
    cbind(rownames(m.cor[g,]),m.cor[g,],final_p)
  })
  names(cor_rank)<-colnames(data.d)
  
  #Format the output
  m.cor_rank<-melt(cor_rank)
  colnames(m.cor_rank)<-c("Env","Var","p","dis")
  toreturn<-cast(m.cor_rank,dis+Env~Var,value="p")
  print(toreturn)
  full[[which(g %in% 7:2)]]<-toreturn
 
  }

stopCluster(cl)

#cast out dataframe

