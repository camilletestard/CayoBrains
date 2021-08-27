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
library(kinship2)
# library(writexl)
library(CePa)
library(doBy)
library(forcats)
library(dplyr)

#Load scan data and population info
setwd("/Users/camilletestard/Documents/GitHub/CayoBrains") 
source("Functions/functions_GlobalNetworkMetrics.R")
source("Functions/KinshipPedigree.R")

setwd("/Users/camilletestard/Dropbox/Cleaned Cayo Data/Raw data") 
bigped <- read.delim("PEDIGREE_2021.txt", sep="\t")

setwd("/Users/camilletestard/Desktop/CayoBrains/Biobanking_info/")
brainweights <- read.csv("CBRU_Brain_Weights.csv")
biobanking_data = read.csv("Cayo Biobank Tissue Catalog_ID_Tissues.csv");names(biobanking_data)[1]="id"
biobanking_data$id = as.character(biobanking_data$id)
biobanking_data$MOM = bigped$BehaviorMom[match(biobanking_data$id, bigped$ID)]

#Add HH dominance for juveniles. Unfortunately these were not calculated for KK or S. 
hh.dominance <- read.csv("HH_Dominance.csv");names(hh.dominance)[1]="id"
hh.dominance$id = as.character(hh.dominance$id)
hh.dominance$id[hh.dominance$id=="2.00E+09"]="2E9"
hh.dominance$id[hh.dominance$id=="2.00E+08"]="2E8"
hh.dominance$id[hh.dominance$id=="7.00E+00"]="7E0"
hh.dominance$id[hh.dominance$id=="7.00E+03"]="7E3"
hh.dominance$id[hh.dominance$id=="8.00E+02"]="8E2"

#Get IDs wanted:
setwd("/Users/camilletestard/Desktop/CayoBrains/")
marina_IDs <- read.csv("hurr_meta_group_year.csv"); unique(marina_IDs$group_year)
IDs = as.data.frame(marina_IDs$animalID); names(IDs)="id"

group = c("KK")
years = c(2018)
groupyears = c("KK2018") #groupyears = c("HH2016","KK2017","S2019")
SocialCapital.ALL = data.frame()

#####################################################################
# Compute Social Capital Metrics, per individual, per year
#####################################################################

gy=1
#for (gy in 1:length(groupyears)){
  
  print(paste("%%%%%%%%%%%%%%%%%% ",groupyears[gy], "%%%%%%%%%%%%%%%%%%"))
  
  #Load data
  setwd("/Users/camilletestard/Desktop/Desktop-Cayo-Maria/Behavioral_Data/Data All Cleaned/")
  allScans <- read.csv("allScans.txt"); allScans<-allScans[allScans$group=="KK" & allScans$year ==2018,]
  meta_data = read.csv(paste("Group",groupyears[gy],"_GroupByYear.txt", sep = ""))
  
  #Create Social Capital Data frame & add Sex, Age, Rank, Group and Year
  SocialCapitalData= meta_data[,c("id","sex","age","ordinal.rank","percofsex.dominanted")]
  names(SocialCapitalData)=c("id","sex","age","ordrank","percentrank")
  SocialCapitalData$group = group[gy]
  SocialCapitalData$year = years[gy]
  # SocialCapitalData$percentrank = SocialCapitalData$percentrank/100 #Make sure rank is in decimals
  
  #####################################################################
  ## For GROOMING DATA
  #####################################################################
  
  ##############################
  #Using a non-network approach: get groom strength IN, OUT and TOTAL & number of grooming partners.
  
  activity.partners = str_split(allScans$partner.ID, c(", "), simplify = TRUE) 
  
  # 1. Output weighted edgelist from the groom data.
  #x = as.character(groom_data$observation_name)
  #groom_data$focalID = toupper(as.character(substr(x,12, 14))) #find focal ID from observation session name
  groom.ID = as.character(meta_data$id)
  groom.give = data.frame(); groom.receive = data.frame(); groom.partner = data.frame(); id=1 #Initialize
  for (id in 1:length(groom.ID)){ #For all IDs
    groom.give[id,"id"] = groom.ID[id]; groom.receive[id,"id"] = groom.ID[id];groom.partner[id, "id"] = groom.ID[id]#Initialize groom give and groom receive df for ID "id"
    groom.give[id,"freq"] = sum(allScans[which(allScans$focalID==groom.ID[id]),"isSocialGive"])
    groom.receive[id,"freq"] = sum(allScans[which(allScans$focalID==groom.ID[id]),"isSocialGet"])# & is a groom receiver
    unq.partner <- unique(as.vector(activity.partners[which(allScans$focalID==groom.ID[id] & allScans$isSocial==1),]))
    unq.partner = unq.partner[unq.partner!=""]; unq.partner = unq.partner[!is.na(unq.partner)]
    groom.partner[id, "numPartners"] = length(unq.partner)
  }#Outputs 3 structures: groom.give and groom receive - all IDs and duration engaged in both states. And groom.partner which gets the number of unique partners per individual.
  #(Note that I may be underestiamting the number of partners since unknowns "juvenile females" are considered all 1 partner)
  
  #IMPORTANT NOTE: I also include gooming done outside of focal sampling. Otherwise there is a mismatch between time spent grooming and other social network measures
  
  #IMPORTANT NOTE 2: There is grooming with unidentified partners or partners outside of the group. This leads to some 
  #discrepancies with social network measures. To get a more accurate picture of number of partners (i.e. including unidentified individuals), I compute it manually instead 
  # of through igraph. 
  
  #GROOM GIVE
  numObs = meta_data$numObs #find the number of hours followed for each groom giver ID
  groom.give$weight <- round(groom.give$freq / numObs, 5) #add weight information by dividing by the #hrs spent observing --> this yields rate
  
  #GROOM RECEIVE
  groom.receive$weight <- round(groom.receive$freq / numObs, 5) #add weight information by dividing by the #hrs spent observing
  
  #NUMBER UNIQUE PARTNERS
  groom.partner$adjusted.numPartners <- round(groom.partner$numPartners/ meta_data$numObs, 5) #Number of grooming partners adjusted for hours observed
  
  # 2. Add Groom weighted in-degree and out-degree (weighted)
  SocialCapitalData$GroomIN = groom.receive$weight #Grooming rate received
  SocialCapitalData$GroomOUT = groom.give$weight #Groming rate given
  SocialCapitalData$GroomIN[is.na(SocialCapitalData$GroomIN)]=0; SocialCapitalData$GroomOUT[is.na(SocialCapitalData$GroomOUT)]=0
  SocialCapitalData$std.GroomIN = SocialCapitalData$GroomIN/mean(SocialCapitalData$GroomIN) #Grooming rate RELATIVE to the group
  SocialCapitalData$std.GroomOUT = SocialCapitalData$GroomOUT/mean(SocialCapitalData$GroomOUT) #Grooming rate RELATIVE to the group
  
  SocialCapitalData$numPartnersGroom = groom.partner$numPartners
  mean_numP = mean(groom.partner$adjusted.numPartners)
  SocialCapitalData$std.numPartnersGroom = groom.partner$adjusted.numPartners/mean_numP
  
  #grooming reciprocity index
  IN=SocialCapitalData$GroomIN; OUT=SocialCapitalData$GroomOUT
  SocialCapitalData$Groom.Recip = 1-(OUT-IN)/(OUT+IN)
  SocialCapitalData$Groom.Recip[SocialCapitalData$Groom.Recip=="NaN"]=1
  #Output index is between 2 (All groom IN, no groom OUT) and 0 (all groom OUT, no IN). 
  # 1-> Groom IN=Groom OUT
  #NaN is when there is no grooming recieved or given.
  
  #Grooming index or CSI
  SocialCapitalData$CSIgroom = SocialCapitalData$GroomIN + SocialCapitalData$GroomOUT
  meanGroomRate = mean(SocialCapitalData$CSIgroom) #for standardization. If we wish to!
  SocialCapitalData$std.CSIgroom = SocialCapitalData$CSIgroom/meanGroomRate
  
  #Measure of social isolation
  #If only have one partner or less
  SocialCapitalData$isIsolated=0; SocialCapitalData$isIsolated[SocialCapitalData$numPartnersGroom<=5]=1
  length(which(SocialCapitalData$isIsolated==1))
  
  ####################################  
  # 4. Create grooming network metrics (i.e. measure of indirect connectedness). This uses a network approach.
  
  #Find all unique IDs
  unqIDs = groom.ID
  
  # Output the Master Edgelist of all possible pairs given the unique IDs.
  masterEL = calcMasterEL_groom(unqIDs)
  el <- calcEdgeList_groom(allScans, masterEL)
  
  el$numObs = rowSums(cbind(meta_data$numObs[match(el$givingID, meta_data$id)], #hours followed is the mean num hrs 
                                        meta_data$numObs[match(el$receivingID, meta_data$id)]))/2 #followed between the pair of each edge
  el$weight =  el$count/el$numObs #groom rate
  #Get a standard weight
  mean_weight = mean(el$weight[which(el$weight != 0)])
  el$stdWeight = el$weight/mean_weight;
  # NAidx = which(is.na(masterEL$stdWeight)); if(length(NAidx)!=0){break} #Check to make sure we don't have NAs
  
  #Compute # weak connections and strong connections
  #Using arbitrary thresholds of half the group mean and 1.5 times the group mean
  EdgeList = el[which(el$weight != 0),c("givingID","receivingID","stdWeight")]
  partners=data.frame(); id=1
  for (id in 1:length(groom.ID)){ #For all IDs
    idx=c(which(EdgeList$givingID==groom.ID[id]), which(EdgeList$receivingID==groom.ID[id]))
    partners[id,"id"]=groom.ID[id]; 
    partners[id,"weak.P"] = length(which(EdgeList$stdWeight[idx]<0.5))
    partners[id,"strong.P"] = length(which(EdgeList$stdWeight[idx]>1.5))
  }
  SocialCapitalData$numWeakPartners = partners$weak.P
  SocialCapitalData$numStrongPartners = partners$strong.P
  
  #Create an adjacency matrix to generate igraph object for social network measures
  weightedEL =el[,c("givingID","receivingID","stdWeight")] #only keep columns of interest
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
  # SocialCapitalData[,c("numPartnersGroom","std.numPartnersGroom","between.groom","eig.cent.groom","clusterCoeff.groom","closeness.groom")] = 
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
  nonzero.el=weightedEL[weightedEL$stdWeight!=0,]
  KC      <- NULL; for(i in 1:length(nonzero.el[,1])){ 
    KC[i] <-  ped[which(rownames(ped)==as.character(nonzero.el$givingID[i])) , which(colnames(ped)==as.character(nonzero.el$receivingID[i]))]
  }
  nonzero.el$KC   <- round(KC, 4)
  nonzero.el$KinPairClass <- "unrelated"
  nonzero.el$KinPairClass[which(nonzero.el$KC >= .125)] <- "rel"
  
  #Find number of kin amonsgt all partners
  kin=data.frame()
  for (id in 1:length(groom.ID)){ #For all IDs
    kin[id,"id"]=groom.ID[id]
    kin[id,"numKin"]=length(which(nonzero.el$givingID==groom.ID[id]|nonzero.el$receivingID==groom.ID[id] & nonzero.el$KinPairClass=="rel"))
    kin[id,"numTotal"]=length(which(nonzero.el$givingID==groom.ID[id]|nonzero.el$receivingID==groom.ID[id]))
    kin[id,"strengthKin"]=sum(nonzero.el$stdWeight[which(nonzero.el$givingID==groom.ID[id]|nonzero.el$receivingID==groom.ID[id] & nonzero.el$KinPairClass=="rel")])
    kin[id,"strengthTotal"]=sum(nonzero.el$stdWeight[which(nonzero.el$givingID==groom.ID[id]|nonzero.el$receivingID==groom.ID[id])])
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
  
  prox_partners = str_split(allScans$in.proximity, c(","), simplify = TRUE) 
  prox_partners = as.data.frame(cbind(as.character(allScans$focalID), prox_partners))
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
  agg.ID <- groom.ID
  agg.give = data.frame(); agg.receive = data.frame(); agg.partner = data.frame(); id=1 #Initialize
  for (id in 1:length(agg.ID)){ #For all IDs
    agg.give[id,"id"] = agg.ID[id]; agg.receive[id,"id"] = agg.ID[id];agg.partner[id, "id"] = agg.ID[id]#Initialize agg give and agg receive df for ID "id"
    agg.give[id,"freq"] = length(which(allScans$focalID==agg.ID[id] & (allScans$focal.activity.isPost=="AG" | allScans$focal.activity.isPost=="SR")))
    agg.receive[id,"freq"] = length(which(allScans$focalID==agg.ID[id] & (allScans$focal.activity.isPost=="AR" | allScans$focal.activity.isPost=="SG")))# & is a agg receiver
    unq.partner <- unique(as.vector(activity.partners[which(allScans$focalID==agg.ID[id] & allScans$focal.activity=="aggression"),]))
    unq.partner = unq.partner[unq.partner!=""]; unq.partner = unq.partner[!is.na(unq.partner)]
    agg.partner[id, "numPartners"] = length(unq.partner)
  }
  #Add weights
  agg.give$weight = agg.give$freq/numObs
  agg.receive$weight = agg.receive$freq/numObs
  
  # 2. Add Aggression weighted in-degree and out-degree (weighted)
  SocialCapitalData$AggOUT = agg.give$weight
  SocialCapitalData$AggIN = agg.receive$weight
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
  #NaN is when there is no aggression recieved or given. I decided to set them to 1 for now (i.e. equal amount given & received). 
  
  #Composite Aggression Index
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
  
  unqIDs = as.character(meta_data$id); sdb=data.frame(matrix(NA,length(unqIDs), 3)); colnames(sdb)=c("id","sdb.freq")
  for (id in 1:length(unqIDs)){
    sdb$id[id] = unqIDs[id]
    sdb$sdb.freq[id] = length(which(allScans$focalID == unqIDs[id] & allScans$focal.activity.isPost == "SD"))
  }
  sdb$sdb.ra = sdb$sdb.freq/numObs
  sdb$sdb.I = sdb$sdb.ra/mean(sdb$sdb.ra)
  
  #Keep the standerdized version
  SocialCapitalData$sdb.ra = sdb$sdb.ra
  SocialCapitalData$std.sdb.ra = sdb$sdb.I
  
  
  ###################################################################
  # Merge and save data
  SocialCapital.ALL = rbind(SocialCapital.ALL, SocialCapitalData)
#}

SocialCapital.ALL$idyear <- paste(SocialCapital.ALL$id, SocialCapital.ALL$year)
#Remove duplicates
SocialCapital.ALL<-SocialCapital.ALL[!duplicated(SocialCapital.ALL$idyear),]

#For now remove data from individuals with no ranking info 
#SocialCapital.ALL=SocialCapital.ALL[SocialCapital.ALL$age>6,]
#SocialCapital.ALL=SocialCapital.ALL[!is.na(SocialCapital.ALL$ordrank),]

#Z-score regressors.
#SocialCapital.ALL[,c(3,5,8:ncol(SocialCapital.ALL))]=scale(SocialCapital.ALL[,c(3,5,8:ncol(SocialCapital.ALL))])


###################################################################
# Select individuals

#Depending on who we have biological data for, we might want to select a subset of individuals.
SocialCapital.ALL=merge(IDs,SocialCapital.ALL, by="id")
setwd('/Users/camilletestard/Desktop/CayoBrains/')
write.csv(SocialCapital.ALL,'ScanBased_SocialMetrics_KK2018.csv', row.names = F)
