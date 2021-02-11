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
#- DSIGroom: (GroomIN + GroomOUT). std.DSIGroom = DSIGroom/mean(DSIGroom) [non-zero mean, for that group and year]
#- Network measures: degree, betweenness, eigenvector centrality and clustering coeff. IMPORTANT NOTE: These are based 
#on standardized weights (i.e. divided by mean)!
#I also add standardized degree, divided by the non-zero degree mean, for that group and year.
#PROXIMITY: 
#- DSIprox = proximity rate = number of proximity partners per scans = total numPartners/numScans
#- Network measures: degree, betweenness, eigenvector centrality and clustering coeff. 
#IMPORTANT NOTE: weights = #proximity events between dyad/hrs followed. Weights are standardized (i.e. divided by non-zero mean)! 
#COMBINATION: numPartners & DSI combining grooming and proximity (taking the average of the two measures)
#AGGRESSION:
#- AggIN, AggOUT: # aggressive interactions/hrs followed. All aggressive interactions are considered here.
#- DSIagg = AggIN + AggOut. Also standardized measure, dividing by group mean for that year.
#SDB & Vigiliance rates:
# sdb or vigiance events / #hrs followed 
#(both standardized and non-standardized measures).

library(stringr)
library(igraph)
library(ggcorrplot)

#Load scan data and population info
setwd("C:/Users/Camille Testard/Documents/GitHub/CayoBrains") 
source("Functions/functions_GlobalNetworkMetrics.R")
setwd("C:/Users/Camille Testard/Desktop/Desktop_CayoBrains/Biobanking_info/")
biobanking_data = read.csv("2016_Biobanking_metadata.csv")
names(biobanking_data)[1]="id"
#Temporarily load dominance for HH here:
dominance = read.csv("HH_Dominance.csv");names(dominance)[1]="id"

biobanking_data$id = as.character(biobanking_data$id)
biobanking_data$id[biobanking_data$id=="2.00E+09"]="2E9"
biobanking_data$id[biobanking_data$id=="7.00E+00"]="7E0"
biobanking_data$id[biobanking_data$id=="7.00E+03"]="7E3"
biobanking_data$id[biobanking_data$id=="8.00E+02"]="8E2"
biobanking_data = biobanking_data[,c(1,9,18,19)]

group = c("HH")
years = c(2016)
groupyears = c("HH2016")
SocialCapital.ALL = data.frame()

#####################################################################
# Compute Social Capital Metrics, per individual, per year
#####################################################################

gy=1
for (gy in 1:length(groupyears)){
  
  print(paste("%%%%%%%%%%%%%%%%%% ",groupyears[gy], "%%%%%%%%%%%%%%%%%%"))
  
  #Load data
  setwd("C:/Users/Camille Testard/Desktop/Desktop-Cayo-Maria/Behavioral_Data/Data All Cleaned") 
  groom_data = read.csv(paste("Group",groupyears[gy],"_GroomingEvents.txt", sep = ""))
  agg_data = read.csv(paste("Group",groupyears[gy],"_AgonisticActions.txt", sep = ""))
  focal_data = read.csv(paste("Group",groupyears[gy],"_FocalData.txt", sep = ""))
  meta_data = read.csv(paste("Group",groupyears[gy],"_GroupByYear.txt", sep = ""))
  prox_data = read.csv(paste("Group",groupyears[gy],"_ProximityGroups.txt", sep = ""))
  
  ## Format data
  
  #Make sure all IDs are in character
  groom_data$groom_giver = as.character(groom_data$groom_giver)
  groom_data$groom_reciever = as.character(groom_data$groom_reciever)
  # agg_data$agonsim.loser = as.character(agg_data$agonsim.loser)
  # agg_data$agonism_winner = as.character(agg_data$agonism_winner)
  
  #Create Social Capital Data frame & add Sex, Age, Rank, Group and Year
  SocialCapitalData= meta_data[,c("id","sex","age","ordinal.rank","percofsex.dominanted")]
  names(SocialCapitalData)=c("id","sex","age","ordrank","percentrank")
  SocialCapitalData$group = group[gy]
  SocialCapitalData$year = years[gy]
  
  #Temporarily use HH_Dominance.csv file to add ranks (until group-by-year files are updated):
  SocialCapitalData[,c("ordrank","percentrank")]=dominance[match(SocialCapitalData$id,dominance$id),
                                                           c("ordinal.rank","percofsex.domianted")]
  
  #####################################################################
  ## For GROOMING DATA
  #####################################################################
  
  ##############################
  #Using a non-network approach
  
  # 1. Output weighted edgelist from the groom data.
  x = as.character(groom_data$observation_name)
  groom_data$focalID = toupper(as.character(substr(x,12, 14))) #find focal ID from observation session name
  groom.ID = as.character(meta_data$id)
  groom.give = data.frame(); groom.receive = data.frame(); id=1 #Initialize
  for (id in 1:length(groom.ID)){ #For all IDs
    groom.give[id,"id"] = groom.ID[id]; groom.receive[id,"id"] = groom.ID[id]; #Initialize groom give and groom receive df for ID "id"
    groom.give[id,"duration"] = sum(groom_data$constrained_duration[which(groom_data$focalID == groom.ID[id] #add groom.give duration if ID is focal
                                                              & groom_data$groom_giver == groom.ID[id])], na.rm =T) # & is a groom giver
    groom.receive[id,"duration"] = sum(groom_data$constrained_duration[which(groom_data$focalID == groom.ID[id] #add groom.give duration if ID is focal
                                                                 & groom_data$groom_reciever == groom.ID[id])], na.rm =T)# & is a groom receiver
  }#Outputs 2 structures: groom.give and groom receive - all IDs and duration engaged in both states.
  
  #GROOM GIVE
  hrs.followed.giver = meta_data$hrs.focalfollowed[match(groom.give$id, meta_data$id)] #find the number of hours followed for each groom giver ID
  groom.give$weight <- round(groom.give$duration / hrs.followed.giver, 5) #add weight information by dividing by the #hrs spent observing --> this yields rate
  
  #GROOM RECEIVE
  hrs.followed.reciever = meta_data$hrs.focalfollowed[match(groom.receive$id, meta_data$id)]
  groom.receive$weight <- round(groom.receive$duration / hrs.followed.reciever, 5) #add weight information by dividing by the #hrs spent observing
  
  # 2. Add Groom weighted in-degree and out-degree (weighted)
  SocialCapitalData$GroomIN = groom.receive$weight[match(meta_data$id,groom.receive$id)]
  SocialCapitalData$GroomOUT = groom.give$weight[match(meta_data$id,groom.give$id)]
  SocialCapitalData$GroomIN[is.na(SocialCapitalData$GroomIN)]=0; SocialCapitalData$GroomOUT[is.na(SocialCapitalData$GroomOUT)]=0
  SocialCapitalData$std.GroomIN=SocialCapitalData$GroomIN/mean(SocialCapitalData$GroomIN)
  SocialCapitalData$std.GroomOUT=SocialCapitalData$GroomOUT/mean(SocialCapitalData$GroomOUT)
  
  #grooming reciprocity index
  IN=SocialCapitalData$GroomIN; OUT=SocialCapitalData$GroomOUT
  SocialCapitalData$Groom.Recip = 1-(OUT-IN)/(OUT+IN)
  SocialCapitalData$Groom.Recip[SocialCapitalData$Groom.Recip=="NaN"]=1
  #Output index is between 2 (All groom IN, no groom OUT) and 0 (all groom OUT, no IN). 
  # 1-> Groom IN=Groom OUT
  #NaN is when there is not grooming recieved or given.
  
  #Grooming index or DSI
  TotalGroom = SocialCapitalData$GroomIN + SocialCapitalData$GroomOUT
  meanGroomRate = mean(TotalGroom) #for standardization. If we wish to!
  SocialCapitalData$DSIgroom = TotalGroom
  SocialCapitalData$std.DSIgroom = TotalGroom/meanGroomRate

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
  mean_weight = mean(masterEL$weight[which(masterEL$weight != 0)])
  masterEL$stdWeight = masterEL$weight/mean_weight;
  NAidx = which(is.na(masterEL$stdWeight)); if(length(NAidx)!=0){break} #Check to make sure we don't have NAs
  
  #Compute # weak connections and strong connections
  EdgeList = masterEL[which(masterEL$weight != 0),c("givingID","receivingID","stdWeight")]
  for (id in 1:length(groom.ID)){ #For all IDs
  }
  
  #Create an adjacency matrix to generate igraph object for social network measures
  weightedEL = masterEL[,c("givingID","receivingID","weight")] #only keep columns of interest
  adjMat = dils::AdjacencyFromEdgelist(weightedEL)
  data = adjMat[["adjacency"]]; rownames(data) = adjMat[["nodelist"]]; colnames(data) = adjMat[["nodelist"]]
  
  #read adjacency matrix
  m=as.matrix(data) # coerces the data set as a matrix
  am.g=graph.adjacency(m,mode="directed",weighted=T) # this will create an directed 'igraph object'. Change qualifiers to make "undirected" or unweighted (null)
  graph = am.g 
  
  #Get the network measures
  NetworkMetrics = data.frame(matrix(NA, nrow = length(V(graph)), ncol = 7)); names(NetworkMetrics)=c("id","deg","INdeg","OUTdeg","between","eig.cent", "clusterCoeff")
  NetworkMetrics$id = as_ids(V(graph))
  #Unweighted degree (undirected)
  NetworkMetrics$deg<-igraph::degree(graph)
  mean_numP = mean(NetworkMetrics$deg, na.rm=T)
  NetworkMetrics$std.deg = NetworkMetrics$deg/mean_numP #standardize partner number
  #Weighted betweenness
  NetworkMetrics$between<-igraph::betweenness(graph, v=V(graph), directed=T)
  #Unweighted eigenvector centrality
  A <-igraph::eigen_centrality(graph, directed=T, scale=T)
  eig.cen = as.data.frame(A["vector"])
  NetworkMetrics$eig.cent = eig.cen$vector
  #Weighted clustering coeff
  NetworkMetrics$clusterCoeff = transitivity(graph, type = "localundirected")
  NetworkMetrics$clusterCoeff[is.nan(NetworkMetrics$clusterCoeff)]=0

  SocialCapitalData[,c("numPartnersGroom","std.numPartnersGroom","between.groom","eig.cent.groom","clusterCoeff.groom")] = 
    NetworkMetrics[match(meta_data$id, NetworkMetrics$id), c("deg","std.deg","between","eig.cent","clusterCoeff")]

  #####################################################################
  ## For PROXIMITY DATA
  #####################################################################
  
  ####################################  
  # Create proximity network metrics. This uses a network approach.
  rscans = prox_data
  
  # Output the Master Edgelist of all possible pairs given the unique IDs.
  masterEL = calcMasterEL(unqIDs)
  
  # 4. Output weighted edgelist from the Master Edgelist.
  options(warn = -1) #set options to ignore all warnings
  weightedEL = calcEdgeList(rscans, masterEL)
  hrs = meta_data$hrs.focalfollowed[match(weightedEL$alter, meta_data$id)]
  weightedEL$weight <- round(weightedEL$count / hrs, 5) #add weight information by dividing by the #hrs spent observing
  weightedEL$count <- NULL;weightedEL$conc <- NULL #delete those calumn variables
  
  #Need to upload an adjacency matrix, rather than socprog style data...
  adjMat = dils::AdjacencyFromEdgelist(weightedEL)
  data = adjMat[["adjacency"]]; rownames(data) = adjMat[["nodelist"]]; colnames(data) = adjMat[["nodelist"]]
  
  #read adjacency matrix
  m=as.matrix(data) # coerces the data set as a matrix
  am.g=graph.adjacency(m,mode="undirected",weighted=T) # this will create an directed 'igraph object'. Change qualifiers to make "undirected" or unweighted (null)
  graph = am.g 
  
  #Get the network measures
  NetworkMetrics = data.frame(matrix(NA, nrow = length(V(graph)), ncol = 7)); names(NetworkMetrics)=c("id","deg","INdeg","OUTdeg","between","eig.cent", "clusterCoeff")
  NetworkMetrics$id = as_ids(V(graph))
  #Unweighted degree (undirected)
  NetworkMetrics$deg<-igraph::degree(graph)
  mean_numP = mean(NetworkMetrics$deg, na.rm=T)
  NetworkMetrics$std.deg = NetworkMetrics$deg/mean_numP #standardize partner number
  #Weighted degree (undirected)
  NetworkMetrics$strength<-igraph::strength(graph)
  mean_strength = mean(NetworkMetrics$strength, na.rm=T)
  NetworkMetrics$std.strength = NetworkMetrics$strength/mean_strength #standardize strength of proximity relationship
  #Unweighted betweenness
  NetworkMetrics$between<-igraph::betweenness(graph, v=V(graph), directed=F, normalized=T)
  #Unweighted eigenvector centrality
  A <-igraph::eigen_centrality(graph, directed=F, scale=T)
  eig.cen = as.data.frame(A["vector"])
  NetworkMetrics$eig.cent = eig.cen$vector
  # #Weighted clustering coeff
  # NetworkMetrics$clusterCoeff = igraph::transitivity(graph, type = "localundirected")
  # NetworkMetrics$clusterCoeff[is.nan(NetworkMetrics$clusterCoeff)]=0
  
  SocialCapitalData[,c("numPartnersProx","std.numPartnersProx","DSIprox","std.DSIprox","between.prox","eig.cent.prox")] = 
    NetworkMetrics[match(meta_data$id, NetworkMetrics$id), c("deg","std.deg","strength","std.strength","between","eig.cent")]
  
  #####################################################################
  ## Combining proximity and grooming data
  #(taking the average of the grooming and proximity data)
  
  SocialCapitalData$numPartners = (SocialCapitalData$numPartnersProx + SocialCapitalData$numPartnersGroom)/2
  SocialCapitalData$std.numPartners = (SocialCapitalData$std.numPartnersProx + SocialCapitalData$std.numPartnersGroom)/2
  SocialCapitalData$DSI = (SocialCapitalData$DSIprox + SocialCapitalData$DSIgroom)/2
  SocialCapitalData$std.DSI = (SocialCapitalData$std.DSIprox + SocialCapitalData$std.DSIgroom)/2
  
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
  SocialCapitalData$std.AggIN=SocialCapitalData$AggIN/mean(SocialCapitalData$AggIN)
  SocialCapitalData$std.AggOUT=SocialCapitalData$AggOUT/mean(SocialCapitalData$AggOUT)
  
  #aggression reciprocity index
  IN=SocialCapitalData$AggIN; OUT=SocialCapitalData$AggOUT
  SocialCapitalData$Agg.Recip = 1-(OUT-IN)/(OUT+IN)
  SocialCapitalData$Agg.Recip[SocialCapitalData$Groom.Recip=="NaN"]=1
  #Output index is between 2 (All aggression IN, no aggression OUT) and 0 (all aggression OUT, no IN). 
  # 1-> Groom IN=Groom OUT
  #NaN is when there is not aggression recieved or given. I decided to set them to 1 for now (i.e. equal amount given & received). 
  
  #Dyadic Aggression Index
  TotalAgg = SocialCapitalData$AggIN + SocialCapitalData$AggOUT
  meanAggRate = mean(TotalAgg)
  SocialCapitalData$DSIAgg = TotalAgg
  SocialCapitalData$std.DSIAgg = TotalAgg/meanAggRate
  
  #####################################################################
  ## Relationship tenor
  Agg=SocialCapitalData$std.DSIAgg
  Aff=SocialCapitalData$std.DSI
  
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
  
  SocialCapitalData$vig.ra = vig.sdb$vig.ra[match(meta_data$id,vig.sdb$id)]
  SocialCapitalData$sdb.ra = vig.sdb$sdb.ra[match(meta_data$id,vig.sdb$id)]
  SocialCapitalData$std.vig.ra = vig.sdb$vig.I[match(meta_data$id,vig.sdb$id)]
  SocialCapitalData$std.sdb.ra = vig.sdb$sdb.I[match(meta_data$id,vig.sdb$id)]
  
  ###################################################################
  # Merge and save data
  SocialCapital.ALL = rbind(SocialCapital.ALL, SocialCapitalData)
  save(SocialCapital.ALL,file ="C:/Users/Camille Testard/Documents/GitHub/Cayo-Maria/R.Data/SocialCapital.RData")
}

#Combine biobanking and social metric data
SocialCapital.ALL[,names(biobanking_data)[-1]] = biobanking_data[match(SocialCapital.ALL$id,biobanking_data$id),c(2:4)]

# #Delete IDs with no rank (there should be very few)
# SocialCapital.ALL = SocialCapital.ALL[!is.na(SocialCapital.ALL$percentrank),]

# Create correlation matrix:
correl_socialCapital = SocialCapital.ALL[,c("age","percentrank","std.DSIgroom","Groom.Recip","std.numPartnersGroom","between.groom","eig.cent.groom",
                                            "std.DSIprox","std.AggIN","Agg.Recip","tenor","std.vig.ra","std.sdb.ra","brain_weight_grams")]
names(correl_socialCapital)=c("age","rank (%)","DSI Groom","Reciprocity","# Grooming partners","betweenness (groom)","eig. centrality (groom)",
                              "DSI Prox","Aggression","Agg. recip.","Tenor","vigilance","self-dir. behav.","brain_weight")

hist(correl_socialCapital[,"# Grooming partners"], breaks=20)

corel_matrix=as.data.frame(cor(correl_socialCapital))
#Get p-value for correlation
p_matrix=as.data.frame(cor_pmat(correl_socialCapital))

#Mask correlation matrix by p-value:
corr = corel_matrix
corr[p_matrix > 0.05]=0
corr[corr == 1]=0

#Visualize
ggcorrplot(corr, method='circle', lab=T,legend.title=c("corr. coeff."))
