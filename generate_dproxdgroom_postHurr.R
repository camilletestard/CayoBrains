#generate_dpGroomProx: Generate change in probability of grooming and proximity over 500 iterations for all individuals or only including individuals with at least X observations.
#p(grooming) -> #scans grooming/total #scans
#  Functions called: CalcSubsampledScans
#Input: allScans.txt and GroupByYear.txt
#Output: ChangeP_minX.RData


#Load functions and data
setwd("~/Documents/Github/Cayo-Maria/")
source("cleaned_code/Functions/CalcSubsampledScans.R")
allScans = read.csv("Data All Cleaned/allScans.txt")

#Set parameters:
num_iter = 500; iter =1
group = c("KK")
years = c(2017)
groupyears = c("KK2017")

#Get IDs wanted:
setwd("/Users/camilletestard/Desktop/CayoBrains/")
marina_IDs <- read.csv("hurr_meta_group_year.csv"); unique(marina_IDs$group_year)
marina_IDs<- marina_IDs[marina_IDs$group_year== "KK_2018",]
IDs = as.data.frame(marina_IDs$animalID); names(IDs)="id"


dprob.ALL = data.frame();
for (iter in 1:num_iter){
  
  print(paste("%%%%%%%%%%%%%%%%%% iter",iter, "%%%%%%%%%%%%%%%%%%"))
  #####################################################################
  # 1. Compute change in p(Acc) and p(Social), per individual, per year
  #####################################################################
  
  #Calculate random subsamples
  randomScans = calcRandomScans(allScans)
  gy=1
  for (gy in 1:length(groupyears)){
    
    rscans = randomScans[which(randomScans$year == years[gy] & randomScans$group == group[gy]),]
    #Load data
    setwd("~/Documents/Github/Cayo-Maria/Data All Cleaned/BehavioralDataFiles/")
    meta_data = read.csv(paste("Group",groupyears[gy],"_GroupByYear.txt", sep = ""))
    
    unqIDs = as.character(meta_data$id)
    dprob=data.frame(matrix(NA, nrow=length(unqIDs),ncol=8)); colnames(dprob)=c(c("id","dpAcc","dpSocial","pACC.pre", "pACC.post", "pSocial.pre", "pSocial.post", "num_obs"))
    for (id in 1:length(unqIDs)){ #For all individuals
      isProx.pre = rscans$isProx[which(as.character(rscans$focalID) == unqIDs[id] & rscans$isPost == 0)] #get all pre-hurricane data for that individuals
      isProx.post = rscans$isProx[which(as.character(rscans$focalID) == unqIDs[id] & rscans$isPost == 1)]#get all post-re-hurricane data for that individuals
      isSocial.pre = rscans$isSocial[which(as.character(rscans$focalID) == unqIDs[id] & rscans$isPost == 0)] #get all pre-hurricane data for that individuals
      isSocial.post = rscans$isSocial[which(as.character(rscans$focalID) == unqIDs[id] & rscans$isPost == 1)]#get all post-re-hurricane data for that individuals
      dpAcc=NA; dpSocial=NA; num_obs = length(isProx.pre) *2
      if (length(isProx.pre)>=10) { #If there are more than 20 observations for that individual pre hurricane
      # This if clause was added as a response to reviewer who worried our results were driven by individuals
      # with few observations. If we only include individuals with at least 40 observations(20 pre-hurricane, 
      # 20 post-hurricane), our results hold. 
      pACC.pre = sum(isProx.pre)/length(isProx.pre)
      pACC.post = sum(isProx.post)/length(isProx.post)
      dpAcc = pACC.post - pACC.pre
      pSocial.pre = sum(isSocial.pre)/length(isSocial.pre)
      pSocial.post = sum(isSocial.post)/length(isSocial.post)
      dpSocial = pSocial.post - pSocial.pre
      } #end of min obs clause
      dprob[id,]=c(unqIDs[id],dpAcc,dpSocial,pACC.pre, pACC.post, pSocial.pre, pSocial.post, num_obs)
    } #end of id for loop
    dprob$group = group[gy]; dprob$year = years[gy]; dprob$iter=iter
    dprob.ALL = rbind(dprob.ALL, dprob)
  } #end of groupyear for loop
}

if (length(which(is.na(dprob.ALL$dpAcc)))!=0) {dprob.ALL = dprob.ALL[-which(is.na(dprob.ALL$dpAcc)),]} #remove NA
unqIDs_minObs<- unique(dprob.ALL$id)

intersectingIDs<-intersect(unqIDs_minObs, as.character(IDs$id))
output.all = data.frame();
for (id in 1:length(intersectingIDs)){
  output <- dprob.ALL[dprob.ALL$id==intersectingIDs[id],]
  output.all<- rbind(output.all, output)
}

# Select individuals
setwd('/Users/camilletestard/Desktop/CayoBrains/')
write.csv(output.all,'ChangeSociality_KK_2017to2018.csv', row.names = F)

