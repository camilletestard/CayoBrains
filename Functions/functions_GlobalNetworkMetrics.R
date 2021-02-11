###############################################
#Functions to generate networks from unique IDs
###############################################
#functions_GlobalNetworkMetrics: functions to generate social networks from unique IDs + compute global network
#metrics & partner preference indices (sex, kin, rank).
#- calcMasterEL: output non-directional EL (for proximity)
#- calcEdgeList: Counts edge weights in non-directional EL (just count, division by #scans 
#happens in "Generate_GlobalNetworkMetrics"
#- calcMasterEl_groom: output directional EL (for grooming)
#- calcEdgeList_groom: counts edge weights in directional EL (just counts as well)
#- createIG / createIG_groom (nondirection and directional): create igraph and tnet objects
#- calcGenStats: calculate generic network values: density, community size, clustering coefficient for whole network,
#female only and male only networks
#- calcSexProps: calculate ratio of observed MM/MF/FF pairs over expected given network (weighted)
#- calcKinProps: calculate ratio of observed ck/dk/unrel pairs over expected given network (weighted)
#- calcRankProps: calculate ratio of observed HH/LL/HL pairs over expected given network (weighted)
#The last 3 functions allow to estimate how much more than chance are each pair category occurring, comparing pre- 
#and post-hurricane.
#Note: expected is computed based on the distribution of sex/kinship/rank attributes in the network IDs.


#1. Output the NON-DIRECTIONAL Master Edgelist of all possible pairs given the unique IDs
calcMasterEL    <- function(unqIDs){ #unq IDs only include focal individuals. I.e. We will only take into account focal individuals in our social networks.
  ego <- NULL; alter <- NULL #Initialize ego and alter
  for(i in 1:length(unqIDs)){ #for all unique IDs
    #Create a list of pairsfor each individual and all other IDs, directionality does not matter
    ego <- append(ego, rep(unqIDs[i], length(unqIDs) - i)) #append: add to the variable "ego"
    alter   <- append(alter  , unqIDs[(i+1):length(unqIDs)])
  }
  alter <- alter[1:length(ego)]#Make sure ego and alter are the same length
  
  masterEL <- data.frame(ego, alter) #combine ego and alter
  masterEL$conc <- paste(masterEL[,1],masterEL[,2],sep=".") #create "pair" or "edge" column 
  masterEL$count <- 0 #initialize count for each pair
  
  return(masterEL)
}

#2. Calculate NON-DIRECTIONAL Edgelist (for proximity data)
calcEdgeList    <- function(rscans, masterEL){
  
  partners = str_split(rscans$in.proximity, c(","), simplify = TRUE) #split proximity partner by ","
  focalID = as.character(rscans$focalID)
  a = cbind(focalID,partners) #bind focal ID and partner together
  
  PP <- NULL  
  for(ii in 1:nrow(a)){ #for all observations
    for(p in 2:ncol(a)){ #for all proximity partners (not counting the first column, which is the focal ID column)
      if(!is.na(a[ii,p])) { #if not NA
        if(a[ii,p] != "" ) { #if not empty
          S1 <- data.frame(ego = as.character(a[ii, 1]), #ego is the focal ID
                           alter = as.character(a[ii, p])) #alter for proximity partners. There will be a separate row for each partner. 
          PP <- dplyr::bind_rows(PP, S1) #bind rows, adds values to PP
        }
      }
    }
  }
  #create an edge in each direction, since the relationship is NOT directional. Both "directions" should count.
  PP$conc1 <- paste(str_trim(PP$alter), str_trim(PP$ego), sep=".") 
  PP$conc2 <- paste(str_trim(PP$ego), str_trim(PP$alter), sep=".")
  head(PP)
  head(masterEL)  
  
  for(i in 1:nrow(PP)){ #for all pairs
    if(PP$conc1[i] %in% masterEL$conc){ #If this edge exists in the master list
      masterEL$count[which(masterEL$conc == PP$conc1[i])] <- masterEL$count[which(masterEL$conc == PP$conc1[i])] +1 #Find the index and add counts
    }
    
    if(PP$conc2[i] %in% masterEL$conc){
      masterEL$count[which(masterEL$conc == PP$conc2[i])] <- masterEL$count[which(masterEL$conc == PP$conc2[i])] +1
    }
  }
  
  return(masterEL)
}

###############################################
#For directed  graph (grooming data)

#1. Output the Master Edgelist of all possible pairs given the unique IDs for grooming data (directional)
calcMasterEL_groom    <- function(unqIDs){ #unq IDs only include focal individuals. I.e. We will take into account focal individuals in our social networks.
  givingID <- NULL; receivingID <- NULL
  for(i in 1:length(unqIDs)){
    givingID <- append(givingID, rep(unqIDs[i], length(unqIDs) - 1)) #append: add to the variable "givingID"
    if (i==1) {receivingID   <- append(receivingID  , unqIDs[2:length(unqIDs)]) }
    else {receivingID   <- append(receivingID  , unqIDs[c(1:i-1,(i+1):length(unqIDs))])  }
  }
  receivingID <- receivingID[1:length(givingID)]
  
  masterEL <- data.frame(givingID,receivingID)
  masterEL$conc <- paste(masterEL[,1],masterEL[,2],sep=".")
  masterEL$count <- 0
  
  return(masterEL)
}

#2. Calculate DIRECTIONAL Edgelist
calcEdgeList_groom    <- function(rscans, masterEL){
  
  partners = str_split(rscans$partner.ID, c(","), simplify = TRUE)
  focalID = as.character(rscans$focalID)
  a = cbind(focalID,partners,rscans$isSocialGive, rscans$isSocialGet,rscans$isSocial)
  
  PP <- NULL;count1=0; count2=0; count3=0  
  for(ii in 1:nrow(a)){ #for all observations
    if(!is.na(a[ii,2])) { #if not NA
      if(a[ii,ncol(a)] == 1 ) { #if isSocial
        if(a[ii,ncol(partners)+2] == 1){count1=count1+1; S1 <- data.frame(givingID = as.character(a[ii, 1]), receivingID = as.character(a[ii, 2]))}
        else if (a[ii,ncol(partners)+3] == 1) {count2=count2+1; S1 <- data.frame(givingID = as.character(a[ii, 2]), receivingID = as.character(a[ii, 1]))}
        else {count3=count3+1}#; S1 <- data.frame(givingID = as.character(a[ii, 2]), receivingID = as.character(a[ii, 1]))}
        PP <- dplyr::bind_rows(PP, S1) #bind rows, adds values to PP
      }
    }
  }
  PP$conc <- paste(str_trim(PP$givingID), str_trim(PP$receivingID), sep=".")#creates directed edges
  head(PP)
  head(masterEL)  
  
  for(i in 1:nrow(PP)){
    if(PP$conc[i] %in% masterEL$conc){ #If this edge exists in the master list
      masterEL$count[which(masterEL$conc == PP$conc[i])] <- masterEL$count[which(masterEL$conc == PP$conc[i])] +1 #Find the index and add counts
    }
  }
  
  return(masterEL)
}

##################################################

#3. Create UNDIRECTED igraph and tnet objects (for proximity networks)
createIG        <- function(edgelist, unqIDs, sexage){
  
  # create igraph object
  ig <- simplify(graph.data.frame(d=edgelist, directed = F), #create UNDIRECTED igraph object 
                 remove.loops = T) #simple graph with no loop
  ig <- delete_edges(ig, which(E(ig)$weight == 0)) #delete edges that have a weight of 0 [delete.edges]
  
  # Add in isolated individuals
  #add.vertices
  ig <- add_vertices(ig, length(unqIDs[which(!unqIDs %in% V(ig)$name)]), #find unique IDs that are not in ig vertices (i.e. isolated individuals whose edge weight are all 0)
                     attr = list(name = as.character(unqIDs[which(!unqIDs %in% V(ig)$name)]))) #add names of those isolated individuals
  
  # Add Sex to IG
  ig <-set_vertex_attr(ig, name = "isFemale", value = sexage$isFemale[match(V(ig)$name, as.character(sexage$id))])#set.vertex.attribute
  
  # Add Age to IG
  ig <- set_vertex_attr(ig, name = "age", value = sexage$age[match(V(ig)$name, as.character(sexage$id))])#set.vertex.attribute
  
  # Create tnet object
  tnet <- cbind(get.edgelist(ig, names=F), #get the list of edges from object ig
                E(ig)$weight) #get the weights of edges
  if(!is_directed(ig)){ #if ig is not directed, which is the case in proximity networks
    tnet <- symmetrise_w(tnet)} #make sure the network is symmetrical (i.e. undirected)
  tnet  <- as.tnet(tnet, type="weighted one-mode tnet") 
  
  # create female only networks
  igFem <-  igraph::delete_vertices(ig, which(V(ig)$isFemale==0))    
  igFem <-  delete_vertex_attr(igFem, name="isFemale")
  
  tnetFem <- cbind(get.edgelist(igFem, names=FALSE), E(igFem)$weight)
  if(!is_directed(igFem)){tnetFem <- symmetrise_w(tnetFem)}
  tnetFem  <- as.tnet(tnetFem, type="weighted one-mode tnet") 
  
  # create male only networks: 
  igMal <-  igraph::delete.vertices(ig, which(V(ig)$isFemale==1))    
  igMal <-  delete_vertex_attr(igMal, name="isFemale")
  
  tnetMal <- cbind(get.edgelist(igMal, names=FALSE), E(igMal)$weight)
  if (length(tnetMal)>0){
    if(!is_directed(igMal)){tnetMal <- symmetrise_w(tnetMal)} 
    tnetMal  <- as.tnet(tnetMal, type="weighted one-mode tnet") 
    netList <- list(ig, tnet, igFem, tnetFem, igMal, tnetMal)}
  else {
    netList <- list(ig, tnet, igFem, tnetFem)
  }
  
  return(netList)
}

#3. Create igraph and tnet objects (for grooming networks)
createIG_groom        <- function(edgelist, unqIDs, sexage){
  
  # create igraph object
  ig <- simplify(graph.data.frame(d=edgelist, directed =T), #create DIRECTED igraph object 
                 remove.loops = T) #simple graph with no loop
  ig <- delete_edges(ig, which(E(ig)$weight == 0)) #delete edges that have a weight of 0 [delete.edges]
  
  # Add in isolated individuals
  #add.vertices
  ig <- add_vertices(ig, length(unqIDs[which(!unqIDs %in% V(ig)$name)]), #find unique IDs that are not in ig vertices (i.e. isolated individuals whose edge weight are all 0)
                     attr = list(name = as.character(unqIDs[which(!unqIDs %in% V(ig)$name)]))) #add names of those isolated individuals
  
  # Add Sex to IG
  ig <-set_vertex_attr(ig, name = "isFemale", value = sexage$isFemale[match(V(ig)$name, as.character(sexage$id))])#set.vertex.attribute
  
  # Add Age to IG
  ig <- set_vertex_attr(ig, name = "age", value = sexage$age[match(V(ig)$name, as.character(sexage$id))])#set.vertex.attribute
  
  # Create tnet object
  tnet <- cbind(get.edgelist(ig, names=F), #get the list of edges from object ig
                E(ig)$weight) #get the weights of edges
  if(!is_directed(ig)){ #if ig is not directed, which is the case in proximity networks
    tnet <- symmetrise_w(tnet)} #make sure the network is symmetrical (i.e. undirected)
  tnet  <- as.tnet(tnet, type="weighted one-mode tnet") 
  
  # create female only networks
  igFem <-  igraph::delete_vertices(ig, which(V(ig)$isFemale==0))    
  igFem <-  delete_vertex_attr(igFem, name="isFemale")
  
  tnetFem <- cbind(get.edgelist(igFem, names=FALSE), E(igFem)$weight)
  if(!is_directed(igFem)){tnetFem <- symmetrise_w(tnetFem)}
  tnetFem  <- as.tnet(tnetFem, type="weighted one-mode tnet") 
  
  # create male only networks: 
  igMal <-  igraph::delete.vertices(ig, which(V(ig)$isFemale==1))    
  igMal <-  delete_vertex_attr(igMal, name="isFemale")
  
  tnetMal <- cbind(get.edgelist(igMal, names=FALSE), E(igMal)$weight)
  if (length(tnetMal)>0){
    if(!is_directed(igMal)){tnetMal <- symmetrise_w(tnetMal)} 
    tnetMal  <- as.tnet(tnetMal, type="weighted one-mode tnet") 
    netList <- list(ig, tnet, igFem, tnetFem, igMal, tnetMal)}
  else {netList <- list(ig, tnet, igFem, tnetFem)}
  
  return(netList)
}

###############################################
#Caculate Global network proporties
###############################################

# Calculate generic network values
calcGenStats <- function(netList){
  
  # total network
  dens    <- length(E(netList[[1]])) / (length(V(netList[[1]]))^2 - length(V(netList[[1]]))) # Density = num edges / all possible edges
  dens.w  <- sum(E(netList[[1]])$weight) * dens #network density * sum of all weights
  kcomm   <- length(fastgreedy.community( #Find dense subgraph, also called communities in graphs
    delete.vertices(as.undirected(netList[[1]]), degree(netList[[1]]) == 0))[]) #delete vertices that have no edges
  clust.w <- as.numeric(tnet::clustering_w(netList[[2]])) #clustering coefficient
  
  # female network
  dens.f    <- length(E(netList[[3]])) / (length(V(netList[[3]]))^2 - length(V(netList[[3]])))
  dens.f.w  <- sum(E(netList[[3]])$weight) * dens.f
  kcomm.f   <- length(fastgreedy.community(delete.vertices(as.undirected(netList[[3]]), degree(netList[[3]]) == 0))[])
  clust.f.w <- as.numeric(tnet::clustering_w(netList[[4]]))
  
  if (length(netList)>4){
    # male network
    dens.m    <- length(E(netList[[5]])) / (length(V(netList[[5]]))^2 - length(V(netList[[5]])))
    dens.m.w  <- sum(E(netList[[5]])$weight) * (length(E(netList[[5]])) / (length(V(netList[[5]]))^2 - length(V(netList[[5]]))))
    kcomm.m   <- length(fastgreedy.community(delete.vertices(as.undirected(netList[[5]]), degree(netList[[5]]) == 0))[]) 
    clust.m.w <- as.numeric(tnet::clustering_w(netList[[6]]))
  }
  else {
    dens.m = NA; dens.m.w = NA ; gini.m = NA ; gini.m.w = NA ; kcomm.m = NA ; clust.m.w = NA
  }
  df <- data.frame(dens, dens.w, kcomm, clust.w,
                   dens.f, dens.f.w, kcomm.f, clust.f.w,
                   dens.m, dens.m.w, kcomm.m, clust.m.w) 
  
  return(df)
}

###############################################
#Caculate partner preference
###############################################

# Calculate sex-based joint-counts
calcSexProps <- function(netList){
  # Classifying pairs as FF, MM or Opposite
  el            <- data.frame(get.edgelist(netList[[1]]), E(netList[[1]])$weight); colnames(el) <- c("ego", "alter", "weight") #tranform tnet object to dataframe
  el$isFemEgo   <- sexage$isFemale[match(as.character(el$ego),   sexage$id)]
  el$isFemAlter <- sexage$isFemale[match(as.character(el$alter), sexage$id)]
  el$pairClass  <- "opp"; el$pairClass[which(el$isFemEgo == 1 & el$isFemAlter == 1)] <- "bothFem"; el$pairClass[which(el$isFemEgo == 0 & el$isFemAlter == 0)] <- "bothMal"
  
  # Calculating Sex Pair weights
  weight.FF     <- sum(el$weight[el$pairClass=="bothFem"]) #sum of all the weights for class FF
  weight.MM     <- sum(el$weight[el$pairClass=="bothMal"]) #sum of all the weights for class MM                       
  weight.cross  <- sum(el$weight[el$pairClass=="opp"]) #sum of all the weights for class MF  
  
  #Compute expected proportions of pair classes if the network was completely random
  possFemPairs        <-    (length(V(netList[[3]]))^2)      - length(V(netList[[3]])) #all possible FF pairs
  possMalPairs        <-    (length(V(netList[[5]]))^2)      - length(V(netList[[5]])) #all possible MM pairs
  allPossiblePairs    <-    (length(V(netList[[1]]))^2)      - length(V(netList[[1]])) #all possible pairs
  possCrossPairs      <-    allPossiblePairs - (possFemPairs + possMalPairs) #all possible MF pairs
  
  #Weighted expected proportions
  exp.FF      <- (possFemPairs/allPossiblePairs)    * sum(el$weight)
  exp.MM      <- (possMalPairs/allPossiblePairs)    * sum(el$weight)
  exp.cross   <- (possCrossPairs/allPossiblePairs)  * sum(el$weight)
  
  #Compute ratio actual/expected. If ratio = 1 then proportion of FF pairs is exactly as expected if network was random
  eo.FF     <- weight.FF/exp.FF
  eo.MM     <- weight.MM/exp.MM
  eo.cross  <- weight.cross/exp.cross
  
  sexPairStats <- data.frame(weight.FF, weight.MM, weight.cross, eo.FF, eo.MM, eo.cross)

  return(sexPairStats)        
}

# Calculate kin-based joint-counts
calcKinProps <- function(netList, ped){
  
  el <- data.frame(get.edgelist(netList[[1]]), E(netList[[1]])$weight); colnames(el) <- c("ego", "alter", "weight") #tranform tnet object to dataframe
  
  KC      <- NULL; for(i in 1:length(el[,1])){ 
    KC[i] <-  ped[which(rownames(ped)==as.character(el$ego[i])) , which(colnames(ped)==as.character(el$alter[i]))]
  }
  el$KC   <- round(KC, 4)
  el$pairClass <- "unrelated"
  el$pairClass[which(el$KC >= .125 & el$KC < .25)] <- "dRel"
  el$pairClass[which(el$KC >= .25)] <- "rel"
  
  #Compute observed weight of related pair classes
  obs.ck     <- sum( el$weight[el$pairClass =="rel"])
  obs.dk     <- sum( el$weight[el$pairClass =="dRel"])
  obs.u      <- sum( el$weight[el$pairClass =="unrelated"])
  
  #Computer expected weight of related pair classes based on the distribution of kinship relationship in the whole population
  ckPairs    <- length(which(ped  >= .25))               ; exp.ck      <- (ckPairs   / length(ped))   * sum(el$weight)
  dkPairs    <- length(which(ped  >= .125 & ped <.25))   ; exp.dk      <- (dkPairs   / length(ped))   * sum(el$weight)
  uPairs     <- length(which(ped  <= .125))              ; exp.u       <- (uPairs    / length(ped))   * sum(el$weight)
  
  #Compute observed/expected ratio
  eo.ck      <- obs.ck    /   exp.ck
  eo.dk      <- obs.dk    /   exp.dk
  eo.u       <- obs.u     /   exp.u
  
  el$weightKC    <- el$weight * el$KC
  kinDegree      <- sum(el$weightKC) / sum(el$weight)
  
  kinPairStats <- data.frame(obs.ck, obs.dk, obs.u, eo.ck, eo.dk, eo.u, kinDegree)
  
  return(kinPairStats)
}

# Calculate Rank-based joint-counts
calcRankProps <- function(netList, dominance_info, year){
  
  #Get edgelist
  el <- data.frame(get.edgelist(netList[[1]]), E(netList[[1]])$weight); colnames(el) <- c("ego", "alter", "weight") #tranform tnet object to dataframe
  
  #Set high rank vs low rank
  dominance_info$ORD_RANK2 = "L"
  dominance_info$ORD_RANK2[which(dominance_info$X.DOMINATED>=70)]="H"
  # Classifying pairs as HH, LL or Opposite
  el$RankEgo   <- dominance_info$ORD_RANK2[match(paste(as.character(el$ego),year,sep=""), as.character(dominance_info$IDyear))]
  el$RankAlter <- dominance_info$ORD_RANK2[match(paste(as.character(el$alter),year,sep=""), as.character(dominance_info$IDyear))]
  el$pairClass  <- "opp"; el$pairClass[which(el$RankEgo == "H" & el$RankAlter == "H")] <- "bothH"; el$pairClass[which(el$RankEgo == "L" & el$RankAlter == "L")] <- "bothL"
  
  # Calculating Sex Pair weights
  weight.HH     <- sum(el$weight[el$pairClass=="bothH"]) #sum of all the weights for class FF
  weight.LL     <- sum(el$weight[el$pairClass=="bothL"]) #sum of all the weights for class MM                       
  weight.crossR  <- sum(el$weight[el$pairClass=="opp"]) #sum of all the weights for class MF  
  
  #Compute expected proportions of pair classes if the network was completely random
  nodes = data.frame(matrix(0, nrow = length(V(netList[[1]])), ncol = 2)); names(nodes)=c("id","rank")
  nodes$id = as_ids(V(netList[[1]]))
  nodes$rank = dominance_info$ORD_RANK2[match(paste(nodes$id,year,sep=""), as.character(dominance_info$IDyear))]
  
  possHHPairs        <-    (length(which(nodes$rank=="H"))^2) - length(which(nodes$rank=="H")) #all possible HH pairs
  possLLPairs        <-    (length(which(nodes$rank=="L"))^2) - length(which(nodes$rank=="L"))  #all possible LL pairs
  allPossiblePairs   <-    (length(V(netList[[1]]))^2)      - length(V(netList[[1]])) #all possible pairs
  possCrossPairs     <-    allPossiblePairs - (possHHPairs + possLLPairs) #all possible HL pairs
  
  #Weighted expected proportions
  exp.HH      <- (possHHPairs/allPossiblePairs)    * sum(el$weight)
  exp.LL      <- (possLLPairs/allPossiblePairs)    * sum(el$weight)
  exp.crossR  <- (possCrossPairs/allPossiblePairs) * sum(el$weight)
  
  #Compute ratio actual/expected. If ratio = 1 then proportion of FF pairs is exactly as expected if network was random
  eo.HH     <- weight.HH/exp.HH
  eo.LL     <- weight.LL/exp.LL
  eo.crossR <- weight.crossR/exp.crossR 
  
  RankPairStats <- data.frame(weight.HH, weight.LL, weight.crossR, eo.HH, eo.LL, eo.crossR)
  
  return(RankPairStats)
}
