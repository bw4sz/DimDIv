#Richness of the two assemblages
richness_To<-richness_sites[names(richness_sites) %in% data.df[x,]$To]
richness_From<-richness_sites[names(richness_sites) %in% data.df[x,]$From]

#Create slots for each community
nullTo<-matrix(nrow=richness_To,ncol=1)
nullFrom<-matrix(nrow=richness_From,ncol=1)

#number of species overlapping
species.To<-sp.list[names(sp.list) %in% data.df[x,]$To][[1]]
species.From<-sp.list[names(sp.list) %in% data.df[x,]$From][[1]]
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
  nullTo[(overlapping+1):richness_To,]<-sample(splist[!splist %in% overlap.species],richness_To-overlapping,replace=FALSE)
  }

#Fill the From community, but cannot pick from any species in the filled To community

speciestopick<-splist[!splist %in% overlap.species & !splist %in% nullTo]

if(!overlapping==richness_From){
  nullFrom[(overlapping+1):richness_From,]<-sample(speciestopick,richness_From-overlapping,replace=FALSE)}

#Create siteXspp table for new assemblages
null_a<-melt(list(nullTo,nullFrom))

#create binary matrix
null_siteXspp<-(!is.na(cast(null_a,L1~value)[,-1]))*1
rownames(null_siteXspp)<-c(as.character(data.df[x,]$To),as.character(data.df[x,]$From))
