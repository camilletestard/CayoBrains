#Compute rank using perc package

#Load library
library(Perc)
library(dplyr)
library(tidyr)
library(stringr)

#Load data
setwd("~/Dropbox (Penn)/CayoBrains/BehaviorData_HH2016/") 
agg_data = read.csv(paste("GroupHH2016_AgonisticActions.txt", sep = ""))
meta_data = read.csv(paste("GroupHH2016_GroupByYear.txt", sep = ""))

# Get weighted edgelist
head(sampleWeightedEdgelist, 5)

edgelist = paste(agg_data$agonism_winner, agg_data$agonism_loser, sep = ".")
edgelist.df = as.data.frame(table(edgelist))
edgelist.df.final = edgelist.df[nchar(as.character(edgelist.df$edgelist))==7,]
edgelist.df.final = edgelist.df.final %>% separate(edgelist, c('winner', 'loser'))

# Convert to win-loss matrix
confmatrix <- as.conflictmat(edgelist.df.final, weighted = TRUE)
conftrans <- transitivity(confmatrix)

# Finding rank order using perc simulations
DominanceProbability <- conductance(confmatrix, maxLength = 4)
# substracting the original conflict matrix from imputed conflict matrix.
informationGain <- DominanceProbability$imputed.conf - confmatrix
# generating a heatmap representing information gained by using informatio from indirect pathways.
plotConfmat(informationGain, ordering = NA, labels = TRUE)
# displaying the converted long format win-loss probability.
dyadicLongConverter(DominanceProbability$p.hat)
# find simRankOrder
s.rank <- simRankOrder(DominanceProbability$p.hat, num = 10, kmax = 10)
s.rank$BestSimulatedRankOrder # displaying the best simulated rank order
s.rank$AllSimulatedRankOrder # all simulated rank order (for 10 tries)
plotConfmat(DominanceProbability$p.hat, ordering = s.rank[[1]]$ID, labels = TRUE) # plot confidence. Brown shows low confidence.

# Find the correlation between perc-generated rank and manual ranking.
perc_rank = s.rank$BestSimulatedRankOrder
rank = meta_data[,c("id","rank")]

for (id in 1:nrow(meta_data)){
  rank$perc_rank[id] = perc_rank$rank[which(!is.na(match(perc_rank$ID,rank$id[id])))]
}
cor.test(rank$rank, rank$perc_rank)

#Save to meta data file
meta_data$perc_rank = rank$perc_rank
meta_data$std.perc_rank = scale(meta_data$perc_rank)
write.csv(meta_data, "GroupHH2016_GroupByYear.csv")
