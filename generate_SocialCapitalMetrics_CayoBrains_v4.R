#Generate pre-hurricane social capital metrics
# This script generates pre-hurricane social capital metrics with in mind the following objectives:
#(1) determining which factors pre-hurricane may predict the change in social rates post-hurricane. Why do some 
#individual change their social rate  by a lot and others only by a little, or even decrease?
#(2) determining which factors pre-hurricane may predict death post-hurricane (did some factor made individuals 
#more likely to die in the year following the hurricane?
#(3) Used as a dataframe for other analyses requiring pre-hurricane social capital factors.
#Functions called: functions_GlobalNetworkMetrics
#Input: GroomingEvents.txt, AgonisticActions.txt, FocalData.txt, GroupByYear.txt, 
#Output Social Capital Metrics: "SocialCapital.RData"
#GROOMING:
#- GroomIN and GroomOUT: based on focal data (so duration/hrs followed). Not standardized. 
#I am  not using the usual "network approach" here but rather using the grooming.txt file durations directly
#- CSIGroom: (GroomIN + GroomOUT). std.CSIGroom = CSIGroom/mean(CSIGroom) [non-zero mean, for that group and year]
#- Network measures: degree, betweenness, eigenvector centrality and clustering coeff. IMPORTANT NOTE: These are based 
#on standardized weights (i.e. divided by mean)!
#I also add standardized degree, divided by the non-zero degree mean, for that group and year.
#PROXIMITY: 
#- CSIprox = proximity rate = number of proximity partners per scans = total numPartners/numScans
#- Network measures: degree, betweenness, eigenvector centrality and clustering coeff. 
#IMPORTANT NOTE: weights = #proximity events between dyad/hrs followed. Weights are standardized (i.e. divided by non-zero mean)! 
#COMBINATION: numPartners & CSI combining grooming and proximity (taking the average of the two measures)
#AGGRESSION:
#- AggIN, AggOUT: # aggressive interactions/hrs followed. All aggressive interactions are considered here.
#- CSIagg = AggIN + AggOut. Also standardized measure, dividing by group mean for that year.
#SDB & Vigiliance rates:
# sdb or vigiance events / #hrs followed 
#(both standardized and non-standardized measures).

library(stringr)
library(igraph)
library(ggcorrplot)
library(ggplot2)
library(ggridges)
library(reshape)
library(reshape2)
# library(writexl)
library(CePa)
library(doBy)

#Load scan data and population info
setwd("/Users/camilletestard/Documents/GitHub/CayoBrains") 
source("Functions/functions_GlobalNetworkMetrics.R")
source("Functions/KinshipPedigree.R")

setwd("/Users/camilletestard/Dropbox/Cleaned Cayo Data/Raw data") 
bigped <- read.delim("PEDIGREE_2021.txt", sep="\t")

#Get and preprocess biobanking data
setwd("/Users/camilletestard/Desktop/CayoBrains/Biobanking_info/")
biobanking_data = read.csv("2016_Biobanking_metadata.csv");names(biobanking_data)[1]="id"
biobanking_data$id = as.character(biobanking_data$id)
biobanking_data$id[biobanking_data$id=="2.00E+09"]="2E9"
biobanking_data$id[biobanking_data$id=="7.00E+00"]="7E0"
biobanking_data$id[biobanking_data$id=="7.00E+03"]="7E3"
biobanking_data$id[biobanking_data$id=="8.00E+02"]="8E2"
biobanking_data = biobanking_data[, c("id","age","days_between","weight_kg", "brain_weight_grams", "trapping_date")]
biobanking_data$id = as.character(biobanking_data$id)
biobanking_data$MOM = bigped$BehaviorMom[match(biobanking_data$id, bigped$AnimalId)]

brainweights <- read.csv("CBRU_Brain_Weights.csv")

#Add HH dominance for juveniles. Unfortunately these were not calculated for KK or S. 
hh.dominance <- read.csv("HH_Dominance.csv");names(hh.dominance)[1]="id"
hh.dominance$id = as.character(hh.dominance$id)
hh.dominance$id[hh.dominance$id=="2.00E+09"]="2E9"
hh.dominance$id[hh.dominance$id=="2.00E+08"]="2E8"
hh.dominance$id[hh.dominance$id=="7.00E+00"]="7E0"
hh.dominance$id[hh.dominance$id=="7.00E+03"]="7E3"
hh.dominance$id[hh.dominance$id=="8.00E+02"]="8E2"

# biobanking_data1 = read.csv("2016_Biobanking_metadata.csv");names(biobanking_data1)[1]="id"
# biobanking_data2 = read.csv("2018_Biobanking_metadata.csv");names(biobanking_data2)[1]="id"
# matched_col1 = match(names(biobanking_data2), names(biobanking_data1));
# matched_col1 = matched_col1[-which(is.na(matched_col1))]
# matched_col2 = match(names(biobanking_data1), names(biobanking_data2));
# matched_col2 = matched_col2[-which(is.na(matched_col2))]
# biobanking_data = rbind(biobanking_data1[matched_col1], biobanking_data2[matched_col2])

gene.expr.IDs = read.csv("GeneExpression_ID.csv");names(gene.expr.IDs)[1]="id"
gene.expr.IDs$id = as.character(gene.expr.IDs$id)
gene.expr.IDs$id[gene.expr.IDs$id=="2.00E+09"]= "2E9"

group = c("HH")
years = c(2016)
groupyears = c("HH2016")
SocialCapital.ALL = data.frame()

#####################################################################
# Compute Social Capital Metrics, per individual, per year
#####################################################################

gy=1
# for (gy in 1:length(groupyears)){

print(paste("%%%%%%%%%%%%%%%%%% ",groupyears[gy], "%%%%%%%%%%%%%%%%%%"))

#Load data
setwd("/Users/camilletestard/Dropbox/Cleaned Cayo Data/Output") 
groom_data = read.csv(paste("Group",groupyears[gy],"_GroomingEvents.txt", sep = ""))
agg_data = read.csv(paste("Group",groupyears[gy],"_AgonisticActions.txt", sep = ""))
focal_data = read.csv(paste("Group",groupyears[gy],"_FocalData.txt", sep = ""))
meta_data = read.csv(paste("Group",groupyears[gy],"_GroupByYear.txt", sep = ""))
prox_data = read.csv(paste("Group",groupyears[gy],"_ProximityGroups.txt", sep = ""))

#Make sure all IDs are in character
groom_data$groom_giver = as.character(groom_data$groom_giver)
groom_data$groom_reciever = as.character(groom_data$groom_reciever)
# agg_data$agonsim.loser = as.character(agg_data$agonsim.loser)
# agg_data$agonism_winner = as.character(agg_data$agonism_winner)

if (groupyears[gy] == "HH2016"){
  #Get rank for HH subadults
  meta_data$percofsex.dominanted <- hh.dominance$percofsex.domianted[match(meta_data$id, hh.dominance$id)]
  meta_data$ordinal.rank <- hh.dominance$ordinal.rank[match(meta_data$id, hh.dominance$id)]
}

#Create Social Capital Data frame & add Sex, Age, Rank, Group and Year
SocialCapitalData= meta_data[,c("id","sex","age","ordinal.rank","percofsex.dominanted")]
names(SocialCapitalData)=c("id","sex","age","ordrank","percentrank")
SocialCapitalData$group = group[gy]
SocialCapitalData$year = years[gy]
# SocialCapitalData$percentrank = SocialCapitalData$percentrank/100 #Make sure rank is in decimals

#Add categorical column in aggressive data file
agg_data$agg_categ[agg_data$agonism_type=="FearGrm"|
                     agg_data$agonism_type=="avoid"|
                     agg_data$agonism_type=="Submit"] = "submissive"
agg_data$agg_categ[agg_data$agonism_type=="threat"|
                     agg_data$agonism_type=="noncontactAgg"|
                     agg_data$agonism_type=="contactAgg"|
                     agg_data$agonism_type=="displace"] = "aggressive"

#####################################################################
## For GROOMING DATA
#####################################################################

##############################
#Using a non-network approach

# 1. Output weighted edgelist from the groom data.
#x = as.character(groom_data$observation_name)
#groom_data$focalID = toupper(as.character(substr(x,12, 14))) #find focal ID from observation session name
groom.ID = as.character(meta_data$id)
groom.give = data.frame(); groom.receive = data.frame(); groom.partner = data.frame(); id=1 #Initialize
for (id in 1:length(groom.ID)){ #For all IDs
  groom.give[id,"id"] = groom.ID[id]; groom.receive[id,"id"] = groom.ID[id]; #Initialize groom give and groom receive df for ID "id"
  groom.give[id,"duration"] = sum(groom_data$constrained_duration[groom_data$groom_giver == groom.ID[id]], na.rm =T) # & is a groom giver
  groom.receive[id,"duration"] = sum(groom_data$constrained_duration[groom_data$groom_reciever == groom.ID[id]], na.rm =T)# & is a groom receiver
  groom.partner[id, "id"] = groom.ID[id]
  groom.partner[id, "numPartners"] = length(unique(c(groom_data$groom_reciever[groom_data$groom_giver == groom.ID[id]], 
                                                     groom_data$groom_giver[groom_data$groom_reciever == groom.ID[id]])))
}#Outputs 3 structures: groom.give and groom receive - all IDs and duration engaged in both states. And groom.partner which gets the number of unique partners per individual.
#(Note that I may be underestiamting the number of partners since unknowns "juvenile females" are considered all 1 partner)

#IMPORTANT NOTE: I also include gooming done outside of focal sampling. Otherwise there is a mismatch between time spent grooming and other social network measures

#IMPORTANT NOTE 2: There is grooming with unidentified partners or partners outside of the group. This leads to some 
#discrepancies with social network measures. To get a more accurate picture of number of partners (i.e. including unidentified individuals), I compute it manually instead 
# of through igraph. 

#GROOM GIVE
hrs.followed.giver = meta_data$hrs.focalfollowed[match(groom.give$id, meta_data$id)] #find the number of hours followed for each groom giver ID
groom.give$weight <- round(groom.give$duration / hrs.followed.giver, 5) #add weight information by dividing by the #hrs spent observing --> this yields rate

#GROOM RECEIVE
hrs.followed.reciever = meta_data$hrs.focalfollowed[match(groom.receive$id, meta_data$id)]
groom.receive$weight <- round(groom.receive$duration / hrs.followed.reciever, 5) #add weight information by dividing by the #hrs spent observing

#NUMBER UNIQUE PARTNERS
groom.partner$adjusted.numPartners <- round(groom.partner$numPartners/ meta_data$hrs.focalfollowed, 5)

# 2. Add Groom weighted in-degree and out-degree (weighted)
SocialCapitalData$GroomIN = groom.receive$weight[match(meta_data$id,groom.receive$id)]
SocialCapitalData$GroomOUT = groom.give$weight[match(meta_data$id,groom.give$id)]
SocialCapitalData$GroomIN[is.na(SocialCapitalData$GroomIN)]=0; SocialCapitalData$GroomOUT[is.na(SocialCapitalData$GroomOUT)]=0
SocialCapitalData$std.GroomIN = SocialCapitalData$GroomIN/mean(SocialCapitalData$GroomIN)
SocialCapitalData$std.GroomOUT = SocialCapitalData$GroomOUT/mean(SocialCapitalData$GroomOUT)

SocialCapitalData$numPartnersGroom = groom.partner$numPartners
mean_numP = mean(groom.partner$adjusted.numPartners)
SocialCapitalData$std.numPartnersGroom = groom.partner$adjusted.numPartners/mean_numP

#grooming reciprocity index
IN=SocialCapitalData$GroomIN; OUT=SocialCapitalData$GroomOUT
SocialCapitalData$Groom.Recip = 1-(OUT-IN)/(OUT+IN)
SocialCapitalData$Groom.Recip[SocialCapitalData$Groom.Recip=="NaN"]=1
#Output index is between 2 (All groom IN, no groom OUT) and 0 (all groom OUT, no IN). 
# 1-> Groom IN=Groom OUT
#NaN is when there is not grooming recieved or given.

#Grooming index or CSI
SocialCapitalData$CSIgroom = SocialCapitalData$GroomIN + SocialCapitalData$GroomOUT
meanGroomRate = mean(SocialCapitalData$CSIgroom) #for standardization. If we wish to!
SocialCapitalData$std.CSIgroom = SocialCapitalData$CSIgroom/meanGroomRate

SocialCapitalData$isIsolated=0; SocialCapitalData$isIsolated[SocialCapitalData$numPartnersGroom<=1]=1
length(which(SocialCapitalData$isIsolated==1))

####################################  
# 4. Create grooming network metrics. This uses a network approach.
#Find all unique IDs
unqIDs = groom.ID

# Output the Master Edgelist of all possible pairs given the unique IDs.
masterEL = calcMasterEL_groom(unqIDs)
groom_data$conc = paste(groom_data$groom_giver,groom_data$groom_reciever,sep=".")
# table(groom_data$conc)
for (ii in 1:length(masterEL$conc)){
  masterEL$count[ii] = sum(groom_data$constrained_duration[which(groom_data$conc == masterEL$conc[ii])], na.rm=T)
}
#Compute grooming rates using hours followed
masterEL$hrs.followed = rowSums(cbind(meta_data$hrs.focalfollowed[match(masterEL$givingID, meta_data$id)], #hours followed is the mean num hrs 
                                      meta_data$hrs.focalfollowed[match(masterEL$receivingID, meta_data$id)]))/2 #followed between the pair of each edge
masterEL$weight =  masterEL$count/masterEL$hrs.followed #groom rate
#Get a standard weight
mean_weight = mean(masterEL$weight[which(masterEL$weight != 0)])
masterEL$stdWeight = masterEL$weight/mean_weight;
# NAidx = which(is.na(masterEL$stdWeight)); if(length(NAidx)!=0){break} #Check to make sure we don't have NAs

#Compute # weak connections and strong connections
EdgeList = masterEL[which(masterEL$weight != 0),c("givingID","receivingID","stdWeight")]
partners=data.frame(); id=1
for (id in 1:length(groom.ID)){ #For all IDs
  idx=c(which(EdgeList$givingID==groom.ID[id]), which(EdgeList$receivingID==groom.ID[id]))
  partners[id,"id"]=groom.ID[id]; 
  partners[id,"weak.P"] = length(which(EdgeList$stdWeight[idx]<0.5))
  partners[id,"strong.P"] = length(which(EdgeList$stdWeight[idx]>1.5))
}

#Create an adjacency matrix to generate igraph object for social network measures
weightedEL = masterEL[,c("givingID","receivingID","stdWeight")] #only keep columns of interest
adjMat = dils::AdjacencyFromEdgelist(weightedEL)
data = adjMat[["adjacency"]]; rownames(data) = adjMat[["nodelist"]]; colnames(data) = adjMat[["nodelist"]]

#read adjacency matrix
m=as.matrix(data) # coerces the data set as a matrix
am.g=graph.adjacency(m,mode="directed",weighted=T) # this will create an directed 'igraph object'. Change qualifiers to make "undirected" or unweighted (null)
graph = am.g 

#Get the network measures
NetworkMetrics = data.frame(matrix(NA, nrow = length(V(graph)), ncol = 7)); names(NetworkMetrics)=c("id","deg","between","eig.cent", "clusterCoeff")
NetworkMetrics$id = as_ids(V(graph))

# #Unweighted degree (undirected)
# NetworkMetrics$deg<-igraph::degree(graph, mode = "total")
# NetworkMetrics$adj.deg<-NetworkMetrics$deg/hrs.followed #account for # hrs followed!
# mean_numP = mean(NetworkMetrics$deg, na.rm=T)
# NetworkMetrics$std.deg = NetworkMetrics$deg/mean_numP #standardize partner number by group average

#Weighted betweenness
NetworkMetrics$between<-igraph::betweenness(graph, v=V(graph), directed=T,normalized=T)
mean_between = mean(NetworkMetrics$between, na.rm=T)
NetworkMetrics$std.between<-NetworkMetrics$between/mean_between

#Weighted eigenvector centrality
A <-igraph::eigen_centrality(graph, directed=T, scale=T)
eig.cen = as.data.frame(A["vector"])
NetworkMetrics$eig.cent = eig.cen$vector
mean_EVC = mean(NetworkMetrics$eig.cent, na.rm=T)
NetworkMetrics$std.eig.cent<-NetworkMetrics$eig.cent/mean_EVC

#Weighted clustering coeff
NetworkMetrics$clusterCoeff = transitivity(graph, type = "weighted", isolates='zero')
mean_clusterCoeff = mean(NetworkMetrics$clusterCoeff, na.rm=T)
NetworkMetrics$std.clusterCoeff<-NetworkMetrics$clusterCoeff/mean_clusterCoeff

#Closeness.
NetworkMetrics$closeness<-igraph::closeness(graph, v=V(graph), mode = "all", normalized=T)
mean_closeness = mean(NetworkMetrics$closeness, na.rm=T)
NetworkMetrics$std.closeness<-NetworkMetrics$closeness/mean_closeness

SocialCapitalData[,c("between.groom","std.between.groom","eig.cent.groom","std.eig.cent.groom","clusterCoeff.groom","std.clusterCoeff.groom",
                     "closeness.groom","std.closeness.groom")] = 
  NetworkMetrics[match(meta_data$id, NetworkMetrics$id), c("between","std.between","eig.cent","std.eig.cent","clusterCoeff","std.clusterCoeff",
                                                           "closeness","std.closeness")] # Not including number of partners...
## SocialCapitalData[,c("numPartnersGroom","std.numPartnersGroom","between.groom","eig.cent.groom","clusterCoeff.groom","closeness.groom")] = 
#   NetworkMetrics[match(meta_data$id, NetworkMetrics$id), c("deg","std.deg","between","eig.cent","clusterCoeff","closeness")]

##################################################################
#Visualize network
#increase space between nodes if overlapping. Choose graph layout.
l <- layout_nicely(am.g,niter=500,area=vcount(am.g)^10,repulserad=vcount(am.g)^8)
V(am.g)$label.cex <- 1  #changes size of labels of vertices

#link up attributes file with network
V(am.g)$sex=as.factor(SocialCapitalData$sex[match(V(am.g)$name,SocialCapitalData$id)]) #sex attribute

#set colour of sexes
V(am.g)$color=V(am.g)$sex #assign the "Sex" attribute as the vertex color
V(am.g)$color=gsub("1","plum1",V(am.g)$color) #Females will be orange
V(am.g)$color=gsub("2","seagreen2",V(am.g)$color) #Males will be lightblue
V(am.g)$color=gsub("0","white",V(am.g)$color) #unknown sex will be white

#set degree attribute
V(am.g)$degree=degree(am.g); V(am.g)$degree=V(am.g)$degree/mean(V(am.g)$degree)
V(am.g)$CSI=SocialCapitalData$std.CSIgroom[match(V(am.g)$name,SocialCapitalData$id)]; V(am.g)$CSI=V(am.g)$CSI/mean(V(am.g)$CSI)
V(am.g)$eigcent= NetworkMetrics$eig.cent; V(am.g)$eigcent[V(am.g)$eigcent>0.9] = 0.4; V(am.g)$eigcent=V(am.g)$eigcent/mean(V(am.g)$eigcent)
V(am.g)$clustercoeff= NetworkMetrics$clusterCoeff; V(am.g)$clustercoeff=V(am.g)$clustercoeff/mean(V(am.g)$clustercoeff)
V(am.g)$closeness= NetworkMetrics$closeness; V(am.g)$closeness=V(am.g)$closeness/mean(V(am.g)$closeness)
V(am.g)$between= NetworkMetrics$between; V(am.g)$between=V(am.g)$between/mean(V(am.g)$between)

#set path for saving and plot graph
setwd("/Users/camilletestard/Desktop/CayoBrains/Results")
postscript(paste(groupyears[gy],"NetworkGraph_Degree.eps"))
#tiff(paste(groupyears[gy],"NetworkGraph_Degree.tiff"), units="in", width=10, height=8, res=300, compression = 'lzw')
plot.igraph(am.g,layout=l, vertex.label=NA, vertex.color=V(am.g)$color, vertex.size=2+3*V(am.g)$degree,edge.color="grey20",
            edge.width=E(am.g)$weight,edge.arrow.size = 0.4, edge.curved=0.5, main = paste("Social Network", groupyears[gy]))
dev.off()

##################################################################
#Kinship Ratio amongst grooming partners

#Compute pedigree for all IDs
pedigree = bigped[,c("AnimalId","Dam","Sire")]
ped <- KinshipPedigree(pedigree)

#Find Kin relationship for existing relationships
el=weightedEL[weightedEL$stdWeight!=0,]
KC      <- NULL; for(i in 1:length(el[,1])){ 
  KC[i] <-  ped[which(rownames(ped)==as.character(el$givingID[i])) , which(colnames(ped)==as.character(el$receivingID[i]))]
}
el$KC   <- round(KC, 4)
el$KinPairClass <- "unrelated"
# el$KinPairClass[which(el$KC >= .125 & el$KC < .25)] <- "dRel"
el$KinPairClass[which(el$KC >= .125)] <- "rel"

#Find number of kin amonsgt all partners
kin=data.frame()
for (id in 1:length(groom.ID)){ #For all IDs
  kin[id,"id"]=groom.ID[id]
  kin[id,"numKin"]=length(which(el$givingID==groom.ID[id]|el$receivingID==groom.ID[id] & el$KinPairClass=="rel"))
  kin[id,"numTotal"]=length(which(el$givingID==groom.ID[id]|el$receivingID==groom.ID[id]))
  kin[id,"strengthKin"]=sum(el$stdWeight[which(el$givingID==groom.ID[id]|el$receivingID==groom.ID[id] & el$KinPairClass=="rel")])
  kin[id,"strengthTotal"]=sum(el$stdWeight[which(el$givingID==groom.ID[id]|el$receivingID==groom.ID[id])])
}
kin$propKin = kin$strengthKin/kin$strengthTotal

SocialCapitalData$Kin = kin$strengthKin[match(meta_data$id,kin$id)]
SocialCapitalData$std.Kin = SocialCapitalData$Kin/mean(SocialCapitalData$Kin)
SocialCapitalData$propKin = kin$propKin[match(meta_data$id,kin$id)]
SocialCapitalData$numKin = kin$numKin[match(meta_data$id,kin$id)]

#####################################################################
## For PROXIMITY DATA
#####################################################################

##############################
#Using a non-network approach

prox_partners = as.data.frame(str_split(prox_data$in.proximity, c(","), simplify = TRUE))
colnames(prox_partners)[1]="focalID"

unqIDs = as.character(meta_data$id)
#Find the number of scans per focal ID
proxRate = data.frame()
for (id in 1:length(unqIDs)){
  proxRate[id, "id"]= unqIDs[id]
  scans = which(prox_partners$focalID == unqIDs[id]) #find the number of scans where focal ID = id
  proxRate[id, "numScans"] = length(scans) #num scans
  proxRate[id, "numPartners"] = 0
  #Find the number of partners
  for (i in 1:length(scans)){ #for all scans of id
    numProxPartners = length(which(prox_partners[scans[i],2:length(prox_partners)] != "")) #find the number of partners at that scan
    proxRate[id, "numPartners"] = proxRate[id, "numPartners"] + numProxPartners #add #partners through scans
  }
  proxRate[id, "proxRate"] = proxRate$numPartners[id]/proxRate$numScans[id] #rate is the average number of partner per proximity scan
}
meanProxRate = mean(proxRate$proxRate, na.rm = T)#to standardize proximity rate later if we wish to.
proxRate$CSIprox = proxRate$proxRate; proxRate$std.CSIprox = proxRate$proxRate/meanProxRate

SocialCapitalData$CSIprox = proxRate$CSIprox[match(meta_data$id,proxRate$id)]
SocialCapitalData$std.CSIprox = proxRate$std.CSIprox[match(meta_data$id,proxRate$id)]

#####################################################################
## Combining proximity and grooming data

SocialCapitalData$CSI = (SocialCapitalData$CSIprox + SocialCapitalData$CSIgroom)/2
SocialCapitalData$std.CSI = (SocialCapitalData$std.CSIprox + SocialCapitalData$std.CSIgroom)/2

#####################################################################
## For AGGRESSION DATA
#####################################################################
#I am including all aggression types and only focal data to get accurate weights

# 1. Output weighted edgelist from the aggression data.
#AGGRESSION GIVE = aggression "winner"
agg.give = as.data.frame(table(agg_data$agonism_winner[which(agg_data$focal_individual=="agonism.winner")]));names(agg.give)[1]='id'
agg.give = agg.give[-which(nchar(as.character(agg.give$id))>3),] #only include known IDs (i.e. 3 characters)
hrs.followed.giver = meta_data$hrs.focalfollowed[match(agg.give$id, meta_data$id)] #find #hours followed
#Exclude individuals not in meta data file
if (length(which(is.na(hrs.followed.giver))) !=0) {
  agg.give = agg.give[-which(is.na(hrs.followed.giver)),]
  hrs.followed.giver = hrs.followed.giver[-which(is.na(hrs.followed.giver))]
}
#Add weights
agg.give$weight = agg.give$Freq/hrs.followed.giver

#AGGRESSION RECEIVE = aggression "loser"
agg.receive = as.data.frame(table(agg_data$agonism_loser[which(agg_data$focal_individual=="agonism.loser")]));names(agg.receive)[1]='id'
agg.receive = agg.receive[-which(nchar(as.character(agg.receive$id))>3),]
hrs.followed.receive= meta_data$hrs.focalfollowed[match(agg.receive$id, meta_data$id)]
#Exclude individuals not in meta data file
if (length(which(is.na(hrs.followed.receive))) !=0) {
  agg.receive = agg.receive[-which(is.na(hrs.followed.receive)),]
  hrs.followed.receive = hrs.followed.receive[-which(is.na(hrs.followed.receive))]
}
#Add weights
agg.receive$weight = agg.receive$Freq/hrs.followed.receive

# 2. Add Aggression weighted in-degree and out-degree (weighted)
SocialCapitalData$AggOUT = agg.give$weight[match(meta_data$id,agg.give$id)]
SocialCapitalData$AggIN = agg.receive$weight[match(meta_data$id,agg.receive$id)]
SocialCapitalData$AggIN[is.na(SocialCapitalData$AggIN)]=0; SocialCapitalData$AggOUT[is.na(SocialCapitalData$AggOUT)]=0
#Standerdize
SocialCapitalData$std.AggIN=SocialCapitalData$AggIN/mean(SocialCapitalData$AggIN)
SocialCapitalData$std.AggOUT=SocialCapitalData$AggOUT/mean(SocialCapitalData$AggOUT)

#aggression reciprocity index
IN=SocialCapitalData$AggIN; OUT=SocialCapitalData$AggOUT
SocialCapitalData$Agg.Recip = 1-(OUT-IN)/(OUT+IN)
SocialCapitalData$Agg.Recip[SocialCapitalData$Agg.Recip=="NaN"]=1
#Output index is between 2 (All aggression IN, no aggression OUT) and 0 (all aggression OUT, no IN). 
# 1-> Agg IN=Agg OUT
#NaN is when there is not aggression recieved or given. I decided to set them to 1 for now (i.e. equal amount given & received). 

#Dyadic Aggression Index
SocialCapitalData$CSIAgg = SocialCapitalData$AggIN + SocialCapitalData$AggOUT
meanAggRate = mean(SocialCapitalData$CSIAgg)
SocialCapitalData$std.CSIAgg = SocialCapitalData$CSIAgg/meanAggRate

#####################################################################
## Relationship tenor
Agg=SocialCapitalData$CSIAgg
Aff=SocialCapitalData$CSIgroom

SocialCapitalData$tenor = 1-(Agg-Aff)/(Aff+Agg)
#Output index is between 2 (All affiliation, no aggression) and 0 (all aggression, no affiliation). 
# 1-> Aggression = affiliation
#NaN is when there is not grooming recieved or given.

#####################################################################
## For vigilance and sdb rates
#####################################################################

unqIDs = as.character(meta_data$id); vig.sdb=data.frame(matrix(NA,length(unqIDs), 3)); colnames(vig.sdb)=c("id","vig.freq","sdb.freq")
for (id in 1:length(unqIDs)){
  vig.sdb$id[id] = unqIDs[id]
  vig.sdb$vig.freq[id] = length(which(focal_data$focal_id == unqIDs[id] & (focal_data$behaviour == "Vigilnce" | focal_data$behaviour == "Vigilnce_overtime")))
  vig.sdb$sdb.freq[id] = length(which(focal_data$focal_id == unqIDs[id] & (focal_data$behaviour == "SelfGrm" | focal_data$behaviour == "SelfGrm_overtime" | focal_data$behaviour == "Scratch" | focal_data$behaviour == "Scratch_overtime")))
}
hrs.followed = meta_data$hrs.focalfollowed[match(unqIDs, meta_data$id)]
vig.sdb$vig.ra = vig.sdb$vig.freq/hrs.followed
vig.sdb$vig.I = vig.sdb$vig.ra/mean(vig.sdb$vig.ra)
vig.sdb$sdb.ra = vig.sdb$sdb.freq/hrs.followed
vig.sdb$sdb.I = vig.sdb$sdb.ra/mean(vig.sdb$sdb.ra)

#Keep the standerdized version
SocialCapitalData$vig.ra = vig.sdb$vig.I[match(meta_data$id,vig.sdb$id)]
SocialCapitalData$sdb.ra = vig.sdb$sdb.I[match(meta_data$id,vig.sdb$id)]


##################################################################
#Measure of loneliness - Mismatch between social interest and social attainment
##################################################################

#Measure social interest:
Loneliness=data.frame(); id=1
for (id in 1:length(unqIDs)){
  Loneliness[id,"id"]=unqIDs[id]
  Loneliness[id,"GrmPrsnt"] = nrow(focal_data[focal_data$focal_id==unqIDs[id] & focal_data$behaviour=="GrmPrsnt",])
  Loneliness[id,"initiator"] = nrow(focal_data[focal_data$focal_id==unqIDs[id] & 
                                                 focal_data$initiator=="focal"& focal_data$behaviour!="Leave",])
  passcont = focal_data$constrained_duration[focal_data$focal_id==unqIDs[id] & focal_data$behaviour=="passcont"]
  Loneliness[id,"passcont"] = ifelse(length(passcont)==0,0,passcont)
}
Loneliness$social.interest=Loneliness$GrmPrsnt+Loneliness$initiator
Loneliness$social.attainment=groom.receive$duration+groom.give$duration+Loneliness$passcont

#Account for differences in #hrs observed
hrs.followed = meta_data$hrs.focalfollowed[match(unqIDs, meta_data$id)]
Loneliness$social.interest=Loneliness$social.interest/hrs.followed
Loneliness$social.attainment=Loneliness$social.attainment/hrs.followed

#Standerdize by group mean
Loneliness$std.social_interest=Loneliness$social.interest/mean(Loneliness$social.interest)
Loneliness$std.social_attainment=Loneliness$social.attainment/mean(Loneliness$social.attainment)

#Get loneliness measure
Loneliness$loneliness = 1-((Loneliness$std.social_attainment-Loneliness$std.social_interest)/
                             (Loneliness$std.social_interest+Loneliness$std.social_attainment))

SocialCapitalData$loneliness = Loneliness$loneliness
SocialCapitalData$social.interest = Loneliness$social.interest
SocialCapitalData$social.attainment = Loneliness$social.attainment
SocialCapitalData$std.social.interest = Loneliness$std.social_interest
SocialCapitalData$std.social.attainment = Loneliness$std.social_attainment
#Notes: Output measure is between 0 and 2, with 2 being the most lonely (all social interest, no social attainment)
# and 0 being the least lonely (all social )


# ##################################################################
# #Measure of social complexity based on Fisher et al 2017
# 
# unqIDs = as.character(meta_data$id); 
# df_submissive=data.frame(matrix(nrow=length(unqIDs), ncol=length(unqIDs))); colnames(df_submssive)=unqIDs; rownames(df_submssive)=unqIDs; 
# df_aggressive=data.frame(matrix(nrow=length(unqIDs), ncol=length(unqIDs))); colnames(df_aggressive)=unqIDs; rownames(df_aggressive)=unqIDs; 
# df_proximity=data.frame(matrix(nrow=length(unqIDs), ncol=length(unqIDs))); colnames(df_proximity)=unqIDs; rownames(df_proximity)=unqIDs; 
# df_grooming=data.frame(matrix(nrow=length(unqIDs), ncol=length(unqIDs))); colnames(df_grooming)=unqIDs; rownames(df_grooming)=unqIDs; 
# df_all=data.frame(matrix(nrow=length(unqIDs), ncol=length(unqIDs))); colnames(df_all)=unqIDs; rownames(df_all)=unqIDs; 
# df_diversity=data.frame(matrix(nrow=length(unqIDs), ncol=length(unqIDs)));   colnames(df_diversity)=unqIDs; rownames(df_diversity)=unqIDs;
# df_CSI=data.frame(matrix(nrow=length(unqIDs), ncol=length(unqIDs))); colnames(df_CSI)=unqIDs; rownames(df_CSI)=unqIDs; 
# df_freq=data.frame(matrix(nrow=length(unqIDs), ncol=length(unqIDs))); colnames(df_freq)=unqIDs; rownames(df_freq)=unqIDs; 
# df_tenor=data.frame(matrix(nrow=length(unqIDs), ncol=length(unqIDs))); colnames(df_tenor)=unqIDs; rownames(df_tenor)=unqIDs; 
# 
# prox_data$in.proximity=as.character(prox_data$in.proximity)
# prox_partners = as.data.frame(str_split(prox_data$in.proximity, c(","), simplify = TRUE))
# colnames(prox_partners)[1]="focalID"
# 
# id1=1; id2=2
# for (id1 in 1:length(unqIDs)){
#   for (id2 in 1:length(unqIDs)){
#     #submission events
#     df_submissive[id1,id2]=length(which((agg_data$agonism_winner==unqIDs[id1] & agg_data$agonism_loser==unqIDs[id2]|
#                         agg_data$agonism_winner==unqIDs[id2] & agg_data$agonism_loser==unqIDs[id1])&
#                           agg_data$agg_categ=="submissive"))
#     #aggression events
#     df_aggressive[id1,id2]=length(which((agg_data$agonism_winner==unqIDs[id1] & agg_data$agonism_loser==unqIDs[id2]|
#                                            agg_data$agonism_winner==unqIDs[id2] & agg_data$agonism_loser==unqIDs[id1])&
#                                           agg_data$agg_categ=="aggressive"))
#     #grooming events
#     df_grooming[id1,id2]=length(which(groom_data$groom_giver==unqIDs[id1] & groom_data$groom_reciever==unqIDs[id2]|
#                                          groom_data$groom_giver==unqIDs[id2] & groom_data$groom_reciever==unqIDs[id1]))
#     #proximity events
#     scans=prox_partners[prox_partners$focalID==unqIDs[id1],2:length(prox_partners)]
#     prox.ids = trimws(scans[scans!=""])
#     df_proximity[id1,id2]=df_proximity[id1,id2]+length(which(str_detect(prox.ids,unqIDs[id2])))
#     df_proximity[id2,id1]=df_proximity[id2,id1]+length(which(str_detect(prox.ids,unqIDs[id2])))
#   }
# }
# df_all = df_submissive+df_aggressive+df_grooming+df_proximity; 
# colnames(df_all)=unqIDs; rownames(df_all)=unqIDs; 
# df_all[df_all==0]=NA
# colSums(df_proximity,na.rm =T) 
# 
# data.all = data.frame(matrix(ncol=5,nrow=70*70)); colnames(data.all)=c("dyad","diversity","CSI","freq","tenor")
# 
# #Behavioral diversity index
# df_diversity = 1/((df_submissive/df_all)^2 + (df_aggressive/df_all)^2 + 
#   (df_grooming/df_all)^2 + (df_proximity/df_all)^2)
# diversity = dils::EdgelistFromAdjacency(df_diversity, nodelist=make.names(unqIDs))
# diversity$conc=paste(diversity$fromnode, diversity$tonode,sep='.')
# 
# data.all[,c("dyad","diversity")]=diversity[,c("conc","weight")]
# data.all$dyad=str_replace(data.all$dyad,"X",""); data.all$dyad=str_replace(data.all$dyad,"X","");
# 
# #Dyadic sociality index
# groom = melt(df_grooming/meta_data$hrs.focalfollowed)
# groom_mean=mean(groom$value[groom$value!=0])
# prox = melt(df_proximity/meta_data$hrs.focalfollowed)
# prox_mean=mean(prox$value[prox$value!=0])
# df_CSI = ((df_grooming/meta_data$hrs.focalfollowed)/groom_mean +
#   (df_proximity/meta_data$hrs.focalfollowed)/prox_mean)/2
# CSI = dils::EdgelistFromAdjacency(df_CSI)
# data.all$CSI=CSI$weight
# 
# #Interaction frequency index
# aggression = melt(df_aggressive/meta_data$hrs.focalfollowed)
# agg_mean=mean(aggression$value[aggression$value!=0])
# submission = melt(df_submissive/meta_data$hrs.focalfollowed)
# sub_mean=mean(submission$value[submission$value!=0])
# 
# df_freq = ((df_grooming/meta_data$hrs.focalfollowed)/groom_mean +
#              (df_proximity/meta_data$hrs.focalfollowed)/prox_mean +
#              (df_aggressive/meta_data$hrs.focalfollowed)/agg_mean +
#              (df_submissive/meta_data$hrs.focalfollowed)/sub_mean)/4
# freq = dils::EdgelistFromAdjacency(df_freq)
# data.all$freq=freq$weight
# 
# #Tenor
# df_tenor = (df_grooming + df_proximity)/(df_aggressive+df_submissive+df_grooming+df_proximity)
# tenor = dils::EdgelistFromAdjacency(df_tenor)
# data.all$tenor=tenor$weight
# 
# #Total number of interactions
# num_interac = dils::EdgelistFromAdjacency(df_all)
# data.all$num = num_interac$weight
# 
# #Only keep dyads with no NA
# data.all=data.all[!is.na(data.all$num),]
# data.minObs = data.all[data.all$num>2,]
# 
# #visualize dyads
# corel_matrix=as.data.frame(cor(data.minObs[,c(2,4,5)]))
# ggcorrplot(corel_matrix, method='circle', lab=T,legend.title=c("corr. coeff."))
# 
# x=data.minObs$diversity; y=data.minObs$freq; z=data.minObs$tenor
# scatter3D(x,y,z, colvar=NULL, col="blue",cex=1)
# plotrgl()
#  
# data.minObs$dyad

###################################################################
# Merge and save data
SocialCapital.ALL = rbind(SocialCapital.ALL, SocialCapitalData)
# }

#Combine biobanking and social metric data
#SocialCapital.ALL[,names(biobanking_data)[-1]] = biobanking_data[match(SocialCapital.ALL$id,biobanking_data$id),c(2:4)]

#Combine with brain weight data
SocialCapital.ALL[,names(brainweights)[c(5,6)]] = brainweights[match(SocialCapital.ALL$id,brainweights$animal_ID),c(5,6)]

# #Z-score regressors.
# SocialCapital.ALL[,c(3,5,8:ncol(SocialCapital.ALL))]=scale(SocialCapital.ALL[,c(3,5,8:ncol(SocialCapital.ALL))], scale = FALSE)


###################################################################
# Select individuals
# Depending on who we have biological data for, we might want to select a subset of individuals.

#Only consider monkeys in the genomic analysis ID list.
#SocialCapital.ALL=merge(gene.expr.IDs,SocialCapital.ALL, by="id")

#Only consider monkeys NOT in the genomic analysis ID list.
# outersect <- function(x, y, ...) {
#   big.vec <- c(x, y, ...)
#   duplicates <- big.vec[duplicated(big.vec)]
#   setdiff(big.vec, unique(duplicates))
# }
# id_list = as.data.frame(outersect(as.character(gene.expr.IDs$id),as.character(SocialCapital.ALL$id)))
# names(id_list)[1]="id"
# SocialCapital.ALL=merge(SocialCapital.ALL,id_list, by="id")

###################################################################
# Select regressors

# #ALL COMPUTED PARAMETERS
# selected_regressors = SocialCapital.ALL[,1:43]
# # names(selected_regressors)=c("id","sex","age","rank","groom OUT","groom IN","#partners","between","EVC","clusterCoeff",
# #                              "closeness","#kin","aggression give","aggression rec.","tenor","vigilance","SDB","loneliness")
# setwd('/Users/camilletestard/Desktop/CayoBrains/')
# write.csv(selected_regressors,'social_metrics_kennyOnly.csv')

#SOCIAL INTEGRATION MODEL
#One measure of "objective" social integration measure and one proxy of "perceived social isolation"
# selected_regressors = SocialCapital.ALL[,c("id","sex","age","percentrank","std.numPartnersGroom","numPartnersGroom",
#                                            "std.AggIN","brain_wt")]
# names(selected_regressors)=c("id","sex","age","rank","#partners","#partners_abs","aggression rec.",
#                              "brain weight")

# #POSITION IN SOCIAL NETWORK MODEL
# #Select a subset which recapitulate the social profile of an individual.
# selected_regressors = SocialCapital.ALL[,c("id","sex","age","percentrank","between.groom",
#                                             "eig.cent.groom","closeness.groom","brain_wt")]
# names(selected_regressors)=c("id","sex","age","rank","betweenness","eig. centrality","closeness","brain weight")

#COMBINATION OF MODELS
# Combined social integration and network position models
selected_regressors = SocialCapital.ALL[,c("id","sex","age","percentrank","std.numPartnersGroom","between.groom",
                                           "eig.cent.groom","closeness.groom","brain_wt")]
names(selected_regressors)=c("id","sex","age","rank","#partners",
                             "betweenness","eig. centrality","closeness","brain weight")

# #SOCIALLY ISOLATED VS. WELL-CONNECTED MALES
# selected_regressors = SocialCapital.ALL[,c("id","sex","age","percentrank","isIsolated","brain_wt")]
# names(selected_regressors)[4]="rank"
# selected_regressors=selected_regressors[selected_regressors$sex=="M",] #select only males
# selected_regressors$sex <-NULL

# #SEX AND AGE MODEL
# selected_regressors = SocialCapital.ALL[,c("id","sex","age","ordrank","days_between")]
# names(selected_regressors)=c("id","sex","age","rank","days_between")

#Make sure regressors are in the correct order of subject
setwd('/Users/camilletestard/Desktop/CayoBrains/DBM_scripts_2021/monkeys')
monkey_list = list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
monkey_list=substr(monkey_list[2:length(monkey_list)],3,5)
monkey_list_df = data.frame(monkey_list); names(monkey_list_df)="id"
output = merge(selected_regressors,monkey_list_df,by="id")
#output<-selected_regressors

#Format for full glm design matrices:
#Format sex column
output$sex=as.character(output$sex)
output$sex[output$sex=="F"]=1; output$sex[output$sex=="M"]=0; output$sex=as.numeric(output$sex)
# #Format ordinal rank column
# output$rank=as.character(output$rank)
# output$rank[output$rank=="L"]=1;output$rank[output$rank=="M"]=2;output$rank[output$rank=="H"]=3;
# output$rank=as.numeric(output$rank)
#Remove column we will not use
output$id=NULL

#Scale
output[,-1]=scale(output[,-1])

# #Sex & age model - Format sex column
# output$sex=as.character(output$sex)
# output$female=0; output$female[output$sex=="F"]=1;
# output$male=0; output$male[output$sex=="M"]=1;
# output$sex=NULL
# #Format ordinal rank column
# output$rank=as.character(output$rank)
# output$rank[output$rank=="L"]=1;output$rank[output$rank=="M"]=2;output$rank[output$rank=="H"]=3;
# output$rank=as.numeric(output$rank)
# #Remove column we will not use
# #output$id=NULL; output=output[,c("female","male","age","rank")]
# #Scale
# output[,c("age","rank")]=scale(output[,c("age","rank")])

# # Social isolated vs. well-connected males
# output$sex=NULL;#output$id=NULL
# output$isIsolated[output$isIsolated<0]=0; output$isIsolated[output$isIsolated>1]=1
# output$isIsolated=as.factor(output$isIsolated)
# ggplot(output, aes(x=isIsolated,y=age, fill=isIsolated))+
#   geom_boxplot()+theme_classic(20)+
#   stat_summary(fun.y=mean, geom="point", shape=23, size=4, colour="black")
# ggsave("isolated_males_age.tiff")

#Add random regressors
#set.seed(0.001) # just to make it reproducible
# #integration model
# output$shuffled_partners = sample(output$`#partners`)
# output$shuffled_loneliness = sample(output$loneliness)
#netpos model
# output$shuffled_ECV = sample(output$`eig. centrality`)
# output$shuffled_closeness = sample(output$closeness)

setwd('/Users/camilletestard/Desktop/CayoBrains/Models')
write.csv(output,'integration_model_notScaled_HH2016.csv',row.names = F)

##################################################################################################
# #Visualize values for each regressor, for each individual. Are there obvious patterns across IDs?:

#Set directory to save results
setwd("/Users/camilletestard/Desktop/CayoBrains/Results")

#Remove outliers: 
output$`eig. centrality`[output$`eig. centrality`>(6)]<-NA; 

# 1. Heatmap
matrix_plot = as.matrix(output[,-1])
row.names(matrix_plot)=output$id
order_id_partners <- as.character(output$id[order(matrix_plot[,4])])

long.form.plot = melt(matrix_plot)
long.form.plot$Var1<- factor(long.form.plot$Var1, levels = order_id_partners)
long.form.plot$Var2<- factor(long.form.plot$Var2, levels = c("sex","age","rank","aggression rec.","betweenness",
                                                             "eig. centrality","closeness","#partners"))

ggplot(long.form.plot, aes(x=factor(Var1),y=Var2,fill=value))+
  geom_tile()+
  scale_fill_gradient(low="yellow", high="blue")+
  #scale_y_discrete(limits = rev(levels(Var2)))+
  xlab("Individual Monkeys")+
  theme_classic(base_size = 20)+
  theme(axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave('SocialParameterDiveristy.eps')
ggsave('SocialParameterDiversity.png')

#Ridgline plot
mdata <- melt(output[,-1]); #mdata[is.na(mdata$value),2]=0
ggplot(mdata, aes(x=value, y=variable, fill=variable))+
  geom_density_ridges(jittered_points = TRUE, position = "raincloud",
                      alpha = 0.7, scale = 0.9)+
  theme_classic(base_size=25)+
  theme(legend.position="none")+
  # ggtitle('Social Brain Model')+
  ylab('')+ xlab('')#+xlim(c(-3, 3))
ggsave('SocialParameterDistribution.eps')
ggsave('SocialParameterDistribution.png')

#Compute and visualize correlation matrix
corel_matrix=as.data.frame(cor(output[,-1]))
#Get p-value for correlation
p_matrix=as.data.frame(cor_pmat(output[,-1]))

#Mask correlation matrix by p-value:
corr = corel_matrix; 
#corr[p_matrix > 0.01]=0
corr[corr == 1]=0

#Statistical test for multicolinearity
# test_model=lm(age~sex+rank+`days between`+`#partners`+
#                 +loneliness+`aggression rec.`+`aggression give`+ shuffled_partners+shuffled_loneliness, data=output)
# test_model=lm(age~sex+rank+Grooming+`days between`+betweenness 
#               +`eig. centrality`+`cluster coeff`+closeness+ shuffled_closeness+shuffled_ECV, data=output)
# performance::check_model(test_model)˜

#Visualize
ggcorrplot(corr, method='circle', lab=T,legend.title=c("corr. coeff."))
ggsave('correlation_matrix_all_Variables.png')
