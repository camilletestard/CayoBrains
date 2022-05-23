#Generate pre-hurricane social capital metrics
# This script generates pre-hurricane social capital metrics with in mind the following objectives:
#(1) 
#Functions called: functions_GlobalNetworkMetrics
#Input: GroomingEvents.txt, AgonisticActions.txt, FocalData.txt, GroupByYear.txt, 
#Output Social Capital Metrics: "SocialCapital.RData"

library(stringr)
library(igraph)
library(ggcorrplot)
library(ggplot2)
library(ggridges)
library(reshape)
library(reshape2)
library(kinship2)
library(doBy)
library(forcats)
library(dplyr)

#Get monkey list
setwd('~/Dropbox (Penn)/CayoBrains/Pre-post-hurr/DBM_pre-post_hurricane/monkeys')
setwd('~/Desktop/DBM_pre-post_hurricane/monkeys')
IDs = list.dirs()
IDs = data.frame(id=as.character(substring(IDs[2:length(IDs)],first=3, last=5)))

#Load scan data and population info
setwd("/Users/camilletestard/Documents/GitHub/CayoBrains") 
source("Functions/functions_GlobalNetworkMetrics.R")
source("Functions/KinshipPedigree.R")

#bigped <- read.delim("PEDIGREE_2021.txt", sep="\t")
#biobanking_data = read.csv("Cayo Biobank Tissue Catalog_ID_Tissues.csv");names(biobanking_data)[1]="id"
#biobanking_data$id = as.character(biobanking_data$id)
#biobanking_data$MOM = bigped$BehaviorMom[match(biobanking_data$id, bigped$ID)]

setwd("/Users/camilletestard/Dropbox (Penn)/CayoBrains/Pre-post-hurr/BiobankData")
biobank2016<- read.csv("2016_Biobanking_Master.csv"); biobank2016$group="HH"; biobank2016$year=2016; names(biobank2016)[1]="id"
biobank2018<- read.csv("2018_Biobanking_Master.csv"); biobank2018$group="KK"; biobank2018$year=2018; names(biobank2018)[1]="id"
biobank2019<- read.csv("2019_Biobanking_Master.csv"); biobank2019$group="S"; biobank2019$year=2019; names(biobank2019)[2]="id"
id.list=rbind(biobank2016[,c("id","group","year")],biobank2018[,c("id","group","year")],biobank2019[,c("id","group","year")])

brainweights <- read.csv("CBRU_Brain_Weights.csv")
brainweights<- brainweights[,c("animal_ID","age","sex","brain_wt")]; names(brainweights)[1]="id"

#Basic hurricane df
temp= merge(brainweights, id.list, by="id")
hurricane.df = merge(temp, IDs, by="id")
hurricane.df$hurricane.status=1;  if (years[gy]<2018){SocialCapitalData$hurricane.status=0;}


#Load population info
setwd("~/Dropbox (Penn)/CayoBehavior/Data/")
bigped <- read.delim("PEDIGREE_2021.txt", sep="\t")

#Add HH dominance for juveniles. Unfortunately these were not calculated for KK or S. 
hh.dominance <- read.csv("HH_Dominance.csv");names(hh.dominance)[1]="id"
hh.dominance$id = as.character(hh.dominance$id)
hh.dominance$id[hh.dominance$id=="2.00E+09"]="2E9"
hh.dominance$id[hh.dominance$id=="2.00E+08"]="2E8"
hh.dominance$id[hh.dominance$id=="7.00E+00"]="7E0"
hh.dominance$id[hh.dominance$id=="7.00E+03"]="7E3"
hh.dominance$id[hh.dominance$id=="8.00E+02"]="8E2"

group = c("HH","KK","KK","S")
years = c(2016,2017,2018,2019)
groupyears = c("HH2016","KK2018","KK2017","S2019")
SocialCapital.ALL = data.frame()

#####################################################################
# Compute Social Capital Metrics, per individual, per year
#####################################################################

gy=1
for (gy in 1:length(groupyears)){
  
  print(paste("%%%%%%%%%%%%%%%%%% ",groupyears[gy], "%%%%%%%%%%%%%%%%%%"))
  
  #Load data
  setwd("~/Dropbox (Penn)/CayoBehavior/Data/BehaviorAllGroups_2010-2019/BehaviouralData/")
  
  #groom_data = read.csv(paste("Group",groupyears[gy],"_GroomingEvents.txt", sep = ""))
  #agg_data = read.csv(paste("Group",groupyears[gy],"_AgonisticActions.txt", sep = ""))
  #focal_data = read.csv(paste("Group",groupyears[gy],"_FocalData.txt", sep = ""))
  meta_data = read.csv(paste("Group",groupyears[gy],"_GroupByYear.txt", sep = ""))
  #prox_data = read.csv(paste("Group",groupyears[gy],"_ProximityGroups.txt", sep = ""))
  

  if (group[gy] == "HH"){
    #Get rank for HH subadults
    meta_data$percofsex.dominanted <- hh.dominance$percofsex.domianted[match(meta_data$id, hh.dominance$id)]
    meta_data$ordinal.rank <- hh.dominance$ordinal.rank[match(meta_data$id, hh.dominance$id)]
  }
  
  #Create Social Capital Data frame & add Sex, Age, Rank, Group and Year
  SocialCapitalData= meta_data[,c("id","ordinal.rank","percofsex.dominanted","hrs.focalfollowed")]
  names(SocialCapitalData)=c("id","ordrank","percentrank","hrs.followed")
  #SocialCapitalData$group = group[gy]
  SocialCapitalData$year = years[gy]
  SocialCapitalData$id.group.year = paste(SocialCapitalData$id, group[gy], years[gy], sep=".")
    
  ###################################################################
  # Merge and save data
  SocialCapital.ALL = rbind(SocialCapital.ALL, SocialCapitalData)
}


###################################################################
# Select individuals

#Depending on who we have biological data for, we might want to select a subset of individuals.
hurricane.social.df=merge(hurricane.df,SocialCapital.ALL, by=c("id"))
hurricane.social.df$year = as.numeric(substr(hurricane.social.df$id.group.year, nchar(hurricane.social.df$id.group.year)-3, nchar(hurricane.social.df$id.group.year)))

setwd("~/Dropbox (Penn)/CayoBehavior/Data/BehaviorAllGroups_2010-2019/BehaviouralData/")
write.csv(hurricane.df,'HurricaneAdult_MetaData.csv', row.names = F)


hurricane.df[is.na(match(hurricane.df$id, hurricane.social.df$id)),]


###################################################################
# Select regressors

# #ALL COMPUTED PARAMETERS
#selected_regressors = SocialCapital.ALL
# selected_regressors = SocialCapital.ALL[,c("id","group","sex","age","percentrank","std.GroomOUT","std.GroomIN","std.CSIgroom",
#                                            "std.numPartnersGroom","Groom.Recip","between.groom","eig.cent.groom","clusterCoeff.groom",
#                                            "closeness.groom","isIsolated","std.Kin","std.CSIprox","std.CSI","std.AggOUT","std.AggIN",
#                                            "Agg.Recip","std.CSIAgg","tenor","vig.ra","sdb.ra","std.social.interest","std.social.attainment","loneliness")]
# names(selected_regressors)=c("id","group","sex","age","rank","groom OUT","groom IN","Total Groom","#partners","Groom reciprocity","betweenness","eigenvector centrality","clusterCoeff",
#                              "closeness","isIsolated","Kin Strength","Proximity rate","Dyadic Sociality Index (prox+groom) ","aggression give","aggression rec.","Agg Recip",
#                              "Total aggresion","tenor","vigilance","self-directed behavior","social interest","social attainment","loneliness")

selected_regressors = SocialCapital.ALL[,c("id","group","sex","age","std.CSI",
                                           "std.numPartnersGroom","isIsolated")]
names(selected_regressors)=c("id","group","sex","age","Dyadic Sociality Index","#Grooming Partners","isIsolated")


setwd('/Users/camilletestard/Desktop/CayoBrains/')
write.csv(selected_regressors,'social_metrics_allAnimals.csv', row.names = F)

#Format for full glm design matrices:
output = selected_regressors[, -c(1,2,3,4)]
# #Format sex column
# output$sex=as.character(output$sex)
# output$sex[output$sex=="F"]=1; output$sex[output$sex=="M"]=0; output$sex=as.numeric(output$sex)
# output[,-1]=scale(output[,-1])

##################################################################################################
# #Visualize values for each regressor, for each individual. Are there obvious patterns across IDs?:

# 1. Heatmap
matrix_plot = as.matrix(output)
row.names(matrix_plot)=SocialCapital.ALL$id
order_id_partners <- as.character(SocialCapital.ALL$id[order(-matrix_plot[,2])])

long.form.plot = melt(matrix_plot)
long.form.plot$Var1<- factor(long.form.plot$Var1, levels = order_id_partners)
long.form.plot$Var2<- factor(long.form.plot$Var2, levels = c("Aggression recieved","Dyadic Sociality Index","isIsolated","#Grooming Partners"))

ggplot(long.form.plot, aes(x=factor(Var1),y=Var2,fill=value))+
  geom_tile()+
  scale_fill_gradient(low="yellow", high="blue")+
  #scale_y_discrete(limits = rev(levels(Var2)))+
  ylab("Sociality Metrics")+
  theme_classic(base_size = 16)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave('Social Metrics HeatMap All Groups 2016-2019.tiff')

#Ridgline plot
mdata <- melt(output); #mdata[is.na(mdata$value),2]=0
mdata$variable<- factor(mdata$variable, levels = c("Aggression recieved","Composite Sociality Index","isIsolated","#Grooming Partners"))

ggplot(mdata, aes(x=value, y=variable, fill=variable))+
  geom_density_ridges(jittered_points = TRUE, position = "raincloud",
                      alpha = 0.5, scale = 0.9)+
  theme_classic(base_size=20)+
  theme(legend.position="none")+
  # ggtitle('Social Brain Model')+
  ylab('')+ xlab('')#+xlim(c(-3, 3))
ggsave('Social Metrics Ridgeline All Groups 2016-2019.tiff')

#Compute and visualize correlation matrix
corel_matrix=as.data.frame(cor(output))
#Get p-value for correlation
p_matrix=as.data.frame(cor_pmat(output))

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
