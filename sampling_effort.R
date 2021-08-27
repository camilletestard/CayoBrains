library(ggplot2)
library(scales)

setwd("/Users/camilletestard/Desktop/CayoBrains/Biobanking_info/")
brainweights <- read.csv("CBRU_Brain_Weights.csv")
biobanking_data = read.csv("Cayo Biobank Tissue Catalog_ID_Tissues.csv");names(biobanking_data)[1]="id"
biobanking_data$id = as.character(biobanking_data$id)

unqIDs = as.character(biobanking_data$id);

#PRE-HURRICANE NUMBER OF GROOMING BOUTS, FOCAL HOURS AND NUMBER OF SCANS PER IDYEAR IN OUR SAMPLE
#Load proximity scans from groups and years of interest in the focal format: 
setwd("/Users/camilletestard/Dropbox/Cleaned Cayo Data/Output")  
group = c("HH","KK","S")
years = c(2016,2017,2019)
groupyears = c("HH2016","KK2017","S2019")

gy=1; grooming_bouts_per_group=data.frame(matrix(data=NA, nrow=length(groupyears), ncol=0)); 
all_grooming_bouts = data.frame(matrix(data=NA, nrow=0, ncol=7)); 
names(all_grooming_bouts) = c("id","groupyear","num.groom.bouts","groom.time.in.s","unq.partners","hrs.followed","num.scans")
mean_hrs_followed = vector(); sd_hrs_followed = vector(); min_hrs_followed = vector();
mean_numScans = vector(); sd_numScans = vector(); min_numScans = vector();
mean_groomingBouts = vector(); sd_groomingBouts = vector()
total_numScans = vector(); total_hrs_followed = vector()
for (gy in 1:length(groupyears)){ #for all groups & years
  
  meta_data = read.csv(paste("Group",groupyears[gy],"_GroupByYear.txt", sep = ""))
  # meta_data = meta_data[!is.na(meta_data$rank),]
  groom_data = read.csv(paste("Group",groupyears[gy],"_GroomingEvents.txt", sep = "")) #load prox data from groupyear gy
  prox_data <- read.csv(paste("Group",groupyears[gy],"_ProximityGroups.txt", sep = ""))
  
  grooming_bouts_per_group$group[gy]=groupyears[gy]
  grooming_bouts_per_group$total_bouts[gy]= nrow(groom_data)
  grooming_bouts_per_group$num_focals[gy]=length(which(!is.na(match(meta_data$id,unqIDs))))
  
  id_list = as.character(meta_data$id[match(unqIDs,meta_data$id)])
  IDs = id_list[!is.na(id_list)]; grooming_bouts = data.frame(matrix(data=NA, nrow=length(IDs), ncol=7)); 
  names(grooming_bouts) = c("id","groupyear","num.groom.bouts","groom.time.in.s","unq.partners","hrs.followed","num.scans");id=1
  for (id in 1:length(IDs)){
    grooming_bouts$id[id] = IDs[id]; grooming_bouts$groupyear = groupyears[gy]
    grooming_bouts$sex[id] <- as.character(meta_data$sex[meta_data$id==IDs[id]])
    grooming_bouts$age[id] <- meta_data$age[meta_data$id==IDs[id]]
    idx=which(groom_data$groom_giver== IDs[id] | groom_data$groom_reciever== IDs[id])
    grooming_bouts$num.groom.bouts[id] = length(idx)
    grooming_bouts$groom.time.in.s[id] = sum(groom_data$constrained_duration[idx])
    grooming_bouts$unq.partners[id] = length(unique(c(as.character(groom_data$groom_giver[idx]), as.character(groom_data$groom_reciever[idx]))))
    grooming_bouts$hrs.followed[id] = meta_data$hrs.focalfollowed[meta_data$id==IDs[id]]
    #grooming_bouts$num.scans[id] = length(which(SubScans$focalID == IDs[id] & SubScans$groupyear == groupyears[gy]))
  }
  
  all_grooming_bouts = rbind(all_grooming_bouts, grooming_bouts)
  
  total_hrs_followed[gy]=sum(grooming_bouts$hrs.followed)
  mean_hrs_followed[gy]=mean(grooming_bouts$hrs.followed) #mean(meta_data$hrs.focalfollowed)
  sd_hrs_followed[gy]=sd(grooming_bouts$hrs.followed)
  min_hrs_followed[gy]=min(grooming_bouts$hrs.followed)
  
  # total_numScans[gy]=sum(grooming_bouts$num.scans)
  # mean_numScans[gy]=mean(grooming_bouts$num.scans) #mean(meta_data$hrs.focalfollowed)
  # sd_numScans[gy]=sd(grooming_bouts$num.scans)
  # min_numScans[gy]=min(grooming_bouts$num.scans)
  
  mean_groomingBouts[gy]=mean(grooming_bouts$num.groom.bouts) #mean(meta_data$hrs.focalfollowed)
  sd_groomingBouts[gy]=sd(grooming_bouts$num.groom.bouts)
}

length(which(all_grooming_bouts$sex=="F"))

mean(all_grooming_bouts$num.groom.bouts); sd(all_grooming_bouts$num.groom.bouts)
mean(all_grooming_bouts$hrs.followed); sd(all_grooming_bouts$hrs.followed) ; range(all_grooming_bouts$hrs.followed); sum(all_grooming_bouts$hrs.followed)
mean(all_grooming_bouts$num.scans); sd(all_grooming_bouts$num.scans); range(all_grooming_bouts$num.scans)

#If exclude some individuals: 
test<-all_grooming_bouts[order(all_grooming_bouts$hrs.followed),]
new_data<-test[27:nrow(test),]
mean(new_data$hrs.followed); sd(new_data$hrs.followed) ; range(new_data$hrs.followed); sum(new_data$hrs.followed)
length(which(new_data$sex=="M"))
range(new_data$age)

#Plot distribution of samples across individuals
ggplot(all_grooming_bouts, aes(x=num.scans, fill=groupyear))+
  geom_histogram(color="#e9ecef")+facet_grid(~groupyear)+ theme_classic(base_size = 20)+ xlab('# scans')+
  theme(axis.text.x=element_text(color = "black", size=15, angle=40, vjust=.8, hjust=0.8)) 

ggplot(all_grooming_bouts, aes(x=hrs.followed, fill=groupyear))+
  geom_histogram(color="#e9ecef")+facet_grid(~groupyear)+ theme_classic(base_size = 20)+ xlab('# hours followed')+
  theme(axis.text.x=element_text(color = "black", size=15, angle=40, vjust=.8, hjust=0.8)) 