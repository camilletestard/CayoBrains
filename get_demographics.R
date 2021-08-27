#Get Demographics
setwd("/Users/camilletestard/Desktop/CayoBrains/Models/")
adults<-read.csv("integration_model.csv",sep = ",")
teens<-read.csv("teenHH_motherBehav_noScaled.csv",sep = ",")
babies<-read.csv("babyHH_motherBehav_notScaled.csv",sep = ",")

adults$age.cat<-NA
adults$age.cat[adults$age<9]<-"young adults"
adults$age.cat[adults$age>8 & adults$age<15]<-"mature adults"
adults$age.cat[adults$age>14]<-"old adults"

adults$age.sex.cat<- paste(adults$age.cat, adults$sex, sep = ".")

output<-as.data.frame(table(adults$age.sex.cat))

table(teens$sex)
table(babies$sex)

all_ID_age<- c(adults$age, babies$age, teens$age)
hist(all_ID_age)

#Check ID overlap between initial dataset and final dataset.
all_IDs<- c(as.character(adults$id), as.character(babies$id), as.character(teens$id))
setwd("/Users/camilletestard/Desktop/CayoBrains/Biobanking_info/")
biobanking_data<-read.csv("2016_Biobanking_metadata.csv")
biobanking_data$animal_id<-as.character(biobanking_data$animal_id)
biobanking_data$animal_id[biobanking_data$animal_id=="2.00E+09"]="2E9"
biobanking_data$animal_id[biobanking_data$animal_id=="7.00E+00"]="7E0"
biobanking_data$animal_id[biobanking_data$animal_id=="7.00E+03"]="7E3"
biobanking_data$animal_id[biobanking_data$animal_id=="8.00E+02"]="8E2"

all_IDs_biobanking<-as.character(biobanking_data$animal_id)

all_IDs_biobanking[!is.element(all_IDs_biobanking, all_IDs)]

#Get number of hours observed per individual:
setwd("/Users/camilletestard/Dropbox/Cleaned Cayo Data/Output") 
meta_data = read.csv(paste("GroupHH2016_GroupByYear.txt", sep = ""))
mean(meta_data$hrs.focalfollowed)
sd(meta_data$hrs.focalfollowed)

#Get dates of data points for HH
focal_data<-read.csv(paste("GroupHH2016_FocalData.txt", sep = ""))
dates<-unique(strtrim(as.character(focal_data$date), 7))
