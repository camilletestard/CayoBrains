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
bigped <- read.delim("C:/Users/Camille Testard/Desktop/Desktop-Cayo-Maria/Behavioral_Data/Data All Raw/PEDIGREE_2021.txt", sep="\t")

#Get and preprocess biobanking data
setwd("C:/Users/Camille Testard/Desktop/Desktop_CayoBrains/Data/Biobanking_info/")
biobanking_data = read.csv("2016_Biobanking_metadata.csv");names(biobanking_data)[1]="id"

biobanking_data$id = as.character(biobanking_data$id)
biobanking_data$id[biobanking_data$id=="2.00E+09"]="2E9"
biobanking_data$id[biobanking_data$id=="7.00E+00"]="7E0"
biobanking_data$id[biobanking_data$id=="7.00E+03"]="7E3"
biobanking_data$id[biobanking_data$id=="8.00E+02"]="8E2"
biobanking_data = biobanking_data[, c("id","age","days_between","weight_kg", "brain_weight_grams", "trapping_date")]

biobanking_data$id = as.character(biobanking_data$id)
biobanking_data$MOM = bigped$BehaviorMom[match(biobanking_data$id, bigped$AnimalId)]

group = c("F","F","F","F","F","F","F","F",
          "HH","HH","KK","KK","KK","R","R","S","S",
          "V","V","V","V")
years = c(2010, 2011, 2012, 2013, 2014, 2015,2016,2017,
          2014, 2016, 2013, 2015, 2017, 2015, 2016, 2011, 2019,
          2015,2016,2017,2019)

meta_data_pooled = data.frame()
for (gy in 1:length(group)){

setwd("C:/Users/Camille Testard/Desktop/Desktop-Cayo-Maria/Behavioral_Data/Data All Cleaned") 
meta_data = read.csv(paste("Group",group[gy], years[gy],"_GroupByYear.txt", sep = ""))

meta_data_pooled = rbind(meta_data_pooled, meta_data)
}

biobanking_data$behavior.id = 0; 
biobanking_data$behavior.id[which(!is.na(match(biobanking_data$id,meta_data_pooled$id)))]=1;

biobanking_data$behavior.mom = 0
biobanking_data$behavior.mom[which(!is.na(match(biobanking_data$MOM,meta_data_pooled$id)))]=1