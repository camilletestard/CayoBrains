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
library(ggplot2)
library(ggridges)
library(reshape)
library(reshape2)
# library(writexl)
library(CePa)
library(doBy)

#Load scan data and population info
setwd("C:/Users/Camille Testard/Documents/GitHub/CayoBrains") 
source("Functions/functions_GlobalNetworkMetrics.R")
source("Functions/KinshipPedigree.R")

setwd("C:/Users/Camille Testard/Desktop/Desktop-Cayo-Maria/") 
bigped <- read.delim("Behavioral_Data/SubjectInfo_2010-2017/PEDIGREE.txt", sep="\t")

setwd("C:/Users/Camille Testard/Desktop/Desktop_CayoBrains/Data/Biobanking_info/")
biobanking_data = read.csv("Cayo Biobank Tissue Catalog_ID_Tissues.csv");names(biobanking_data)[1]="id"

biobanking_data$id = as.character(biobanking_data$id)
biobanking_data$MOM = bigped$BehaviorMom[match(biobanking_data$id, bigped$ID)]

IDs = as.data.frame(biobanking_data$id); names(IDs)="id"

group = c("HH","KK","S")
years = c(2016,2017,2019)
groupyears = c("HH2016","KK2017","S2019")
SocialCapital.ALL = data.frame()

#####################################################################
# Compute Social Capital Metrics, per individual, per year
#####################################################################

gy=3
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
  SocialCapitalData$percentrank = SocialCapitalData$percentrank/100 #Make sure rank is in decimals
  
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
  x = as.character(groom_data$observation_name)
  groom_data$focalID = toupper(as.character(substr(x,12, 14))) #find focal ID from observation session name
  groom.ID = as.character(meta_data$id)
  groom.give = data.frame(); groom.receive = data.frame(); id=1 #Initialize
  for (id in 1:length(groom.ID)){ #For all IDs
    groom.give[id,"id"] = groom.ID[id]; groom.receive[id,"id"] = groom.ID[id]; #Initialize groom give and groom receive df for ID "id"
    groom.give[id,"duration"] = sum(groom_data$constrained_duration[groom_data$groom_giver == groom.ID[id]], na.rm =T) # & is a groom giver
    groom.receive[id,"duration"] = sum(groom_data$constrained_duration[groom_data$groom_reciever == groom.ID[id]], na.rm =T)# & is a groom receiver
  }#Outputs 2 structures: groom.give and groom receive - all IDs and duration engaged in both states.
  #Note: I also include gooming done outside of focal sampling. Otherwise there is a mismatch between time spent grooming and other social network measures
  #Note 2: There is grooming with unidentified partners or partners outside of the group. This leads to some 
  #discrepancies with social network measures.
  
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
  SocialCapitalData$DSIgroom = SocialCapitalData$GroomIN + SocialCapitalData$GroomOUT
  meanGroomRate = mean(SocialCapitalData$DSIgroom) #for standardization. If we wish to!
  SocialCapitalData$std.DSIgroom = SocialCapitalData$DSIgroom/meanGroomRate
  
  SocialCapitalData$isIsolated=0; SocialCapitalData$isIsolated[SocialCapitalData$std.DSIgroom<0.3]=1

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
  
  # #Compute # weak connections and strong connections
  # EdgeList = masterEL[which(masterEL$weight != 0),c("givingID","receivingID","stdWeight")]
  # partners=data.frame(); id=1
  # for (id in 1:length(groom.ID)){ #For all IDs
  #   idx=c(which(EdgeList$givingID==groom.ID[id]), which(EdgeList$receivingID==groom.ID[id]))
  #   partners[id,"id"]=groom.ID[id]; strong_partners[id,"id"]=groom.ID[id];
  #   partners[id,"weak.P"] = length(which(EdgeList$stdWeight[idx]<0.5))
  #   partners[id,"strong.P"] = length(which(EdgeList$stdWeight[idx]>1.5))
  # }
  
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
  #Unweighted degree (undirected)
  NetworkMetrics$deg<-igraph::degree(graph)
  hrs.followed = meta_data$hrs.focalfollowed[match(NetworkMetrics$id, meta_data$id)]
  NetworkMetrics$deg/hrs.followed #account for # hrs followed!
  mean_numP = mean(NetworkMetrics$deg, na.rm=T)
  NetworkMetrics$std.deg = NetworkMetrics$deg/mean_numP #standardize partner number by group average
  #Weighted betweenness
  NetworkMetrics$between<-igraph::betweenness(graph, v=V(graph), directed=T,normalized=T)
  #Weighted eigenvector centrality
  A <-igraph::eigen_centrality(graph, directed=T, scale=T)
  eig.cen = as.data.frame(A["vector"])
  NetworkMetrics$eig.cent = eig.cen$vector
  #Weighted clustering coeff
  NetworkMetrics$clusterCoeff = transitivity(graph, type = "weighted", isolates='zero')
  #Closeness.
  NetworkMetrics$closeness<-igraph::closeness(graph, v=V(graph), normalized=T)
  
  SocialCapitalData[,c("numPartnersGroom","std.numPartnersGroom","between.groom","eig.cent.groom","clusterCoeff.groom","closeness.groom")] = 
    NetworkMetrics[match(meta_data$id, NetworkMetrics$id), c("deg","std.deg","between","eig.cent","clusterCoeff","closeness")]

  ##################################################################
  #Visualize network
  #increase space between nodes if overlapping. Choose graph layout.
  l <- layout_nicely(am.g,niter=500,area=vcount(am.g)^6,repulserad=vcount(am.g)^6) 
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
  V(am.g)$DSI=SocialCapitalData$std.DSIgroom[match(V(am.g)$name,SocialCapitalData$id)]; V(am.g)$DSI=V(am.g)$DSI/mean(V(am.g)$DSI)
  V(am.g)$eigcent= NetworkMetrics$eig.cent; V(am.g)$eigcent[V(am.g)$eigcent>0.9] = 0.4; V(am.g)$eigcent=V(am.g)$eigcent/mean(V(am.g)$eigcent)
  V(am.g)$clustercoeff= NetworkMetrics$clusterCoeff; V(am.g)$clustercoeff=V(am.g)$clustercoeff/mean(V(am.g)$clustercoeff)
  V(am.g)$closeness= NetworkMetrics$closeness; V(am.g)$closeness=V(am.g)$closeness/mean(V(am.g)$closeness)
  V(am.g)$between= NetworkMetrics$between; V(am.g)$between=V(am.g)$between/mean(V(am.g)$between)
  
  #set path for saving and plot graph
  setwd("C:/Users/Camille Testard/Desktop/Desktop_CayoBrains/Results") 
  tiff("NetworkGraph_DSI.tiff", units="in", width=10, height=8, res=300, compression = 'lzw')
  plot.igraph(am.g,layout=l, vertex.label=V(am.g)$name, vertex.color=V(am.g)$color, vertex.size=5*V(am.g)$DSI,edge.color="grey20", 
              edge.width=E(am.g)$weight*2,edge.arrow.size = 0.5, main = "Social Network HH2016")
  dev.off()

  ##################################################################
  #Kinship Ratio amongst grooming partners
  
  #Compute pedigree for all IDs
  pedigree = bigped[,c("ID","DAM","SIRE")]
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
  SocialCapitalData$Agg.Recip[SocialCapitalData$Groom.Recip=="NaN"]=1
  #Output index is between 2 (All aggression IN, no aggression OUT) and 0 (all aggression OUT, no IN). 
  # 1-> Groom IN=Groom OUT
  #NaN is when there is not aggression recieved or given. I decided to set them to 1 for now (i.e. equal amount given & received). 
  
  #Dyadic Aggression Index
  SocialCapitalData$DSIAgg = SocialCapitalData$AggIN + SocialCapitalData$AggOUT
  meanAggRate = mean(SocialCapitalData$DSIAgg)
  SocialCapitalData$std.DSIAgg = SocialCapitalData$DSIAgg/meanAggRate
  
  #####################################################################
  ## Relationship tenor
  Agg=SocialCapitalData$DSIAgg
  Aff=SocialCapitalData$DSIgroom
  
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
  social_interest=Loneliness$social.interest/mean(Loneliness$social.interest)
  social_attainment=Loneliness$social.attainment/mean(Loneliness$social.attainment)
  
  #Get loneliness measure
  SocialCapitalData$loneliness = 1-((social_attainment-social_interest)/(social_interest+social_attainment))
  SocialCapitalData$social.interest = social_interest
  SocialCapitalData$social.attainment = social_attainment
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
  # df_DSI=data.frame(matrix(nrow=length(unqIDs), ncol=length(unqIDs))); colnames(df_DSI)=unqIDs; rownames(df_DSI)=unqIDs; 
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
  # data.all = data.frame(matrix(ncol=5,nrow=70*70)); colnames(data.all)=c("dyad","diversity","DSI","freq","tenor")
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
  # df_DSI = ((df_grooming/meta_data$hrs.focalfollowed)/groom_mean +
  #   (df_proximity/meta_data$hrs.focalfollowed)/prox_mean)/2
  # DSI = dils::EdgelistFromAdjacency(df_DSI)
  # data.all$DSI=DSI$weight
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
}

#For now remove data from individuals with no ranking info 
#SocialCapital.ALL=SocialCapital.ALL[SocialCapital.ALL$age>6,]
SocialCapital.ALL=SocialCapital.ALL[!is.na(SocialCapital.ALL$ordrank),]

#Z-score regressors.
SocialCapital.ALL[,c(3,5,8:ncol(SocialCapital.ALL))]=scale(SocialCapital.ALL[,c(3,5,8:ncol(SocialCapital.ALL))])


###################################################################
# Select individuals

#Depending on who we have biological data for, we might want to select a subset of individuals.
SocialCapital.ALL=merge(IDs,SocialCapital.ALL, by="id")

###################################################################
# Select regressors

# #ALL COMPUTED PARAMETERS
selected_regressors = SocialCapital.ALL[,c("id","sex","age","percentrank","GroomOUT","GroomIN",
                                           "numPartnersGroom","between.groom","eig.cent.groom","clusterCoeff.groom",
                                           "closeness.groom","numKin","AggOUT","AggIN","tenor","vig.ra","sdb.ra","loneliness")]
names(selected_regressors)=c("id","sex","age","rank","groom OUT","groom IN","#partners","between","EVC","clusterCoeff",
                             "closeness","#kin","aggression give","aggression rec.","tenor","vigilance","SDB","loneliness")
setwd('C:/Users/Camille Testard/Desktop/Desktop_CayoBrains/')
write.csv(selected_regressors,'social_metrics_allAnimals.csv')

#Format for full glm design matrices:
output = selected_regressors
#Format sex column
output$sex=as.character(output$sex)
output$sex[output$sex=="F"]=1; output$sex[output$sex=="M"]=0; output$sex=as.numeric(output$sex)
#Scale
output[,c("sex","rank")]=scale(output[,c("sex","rank")])


##################################################################################################
# #Visualize values for each regressor, for each individual. Are there obvious patterns across IDs?:

# 1. Heatmap
matrix_plot = as.matrix(output[,-1])
row.names(matrix_plot)=SocialCapital.ALL$id
# heatmap(matrix_plot,Colv = NA, Rowv = NA,)

long.form.plot = melt(matrix_plot)
ggplot(long.form.plot,aes(x=Var1,y=Var2,fill=value))+
  geom_tile()+
  scale_fill_gradient(low="yellow", high="blue")

#Ridgline plot
mdata <- melt(output[,-1]); #mdata[is.na(mdata$value),2]=0
ggplot(mdata, aes(x=value, y=variable, fill=variable))+
  geom_density_ridges(jittered_points = TRUE, position = "raincloud",
                      alpha = 0.7, scale = 0.9)+
  theme_classic(base_size=25)+
  theme(legend.position="none")+
  # ggtitle('Social Brain Model')+
  ylab('')+ xlab('')#+xlim(c(-3, 3))
ggsave('Social Metrics All Culled Animals 2016-2019.tiff')

#Compute and visualize correlation matrix
corel_matrix=as.data.frame(cor(output[,-1]))
#Get p-value for correlation
p_matrix=as.data.frame(cor_pmat(output[,-1]))

#Mask correlation matrix by p-value:
corr = corel_matrix
# corr[p_matrix > 0.01]=0
# corr[corr == 1]=0

#Visualize
ggcorrplot(corr, method='circle', lab=T,legend.title=c("corr. coeff."))
ggsave('AllCulledAnimals_CorrPlot.tiff')

#Separated by sex
data_male = output[output$sex<0,]; 
corel_matrix=as.data.frame(cor(data_male[,-c(1,2)]))
ggcorrplot(corel_matrix, method='circle', lab=T,legend.title=c("corr. coeff."))
ggsave('MalesCulled_CorrPlot.tiff')

data_female = output[output$sex>0,];
corel_matrix=as.data.frame(cor(data_female[,-c(1,2)]))
ggcorrplot(corel_matrix, method='circle', lab=T,legend.title=c("corr. coeff."))
ggsave('FemalesCulled_CorrPlot.tiff')
