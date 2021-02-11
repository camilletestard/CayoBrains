# Set parameters
num_subjects = 67
num_voxels = 187803 #size of mask
effect_size = 0 #i.e. no difference
thr = 0.001 #uncorrected p
cluster_size = 200 # number of contiguous voxels

#Create random values from normal distribution
group1 = rnorm(num_subjects*num_voxels, mean=0, sd=2)
group1=matrix(group1, nrow=num_subjects, byrow = T)

group2 = rnorm(num_subjects*num_voxels, mean=0+effect_size, sd=2)
group2 = matrix(group2, nrow=num_subjects, byrow = T)

#Check both matrices are not equal
sum(abs(group1-group2))

#Run ttest for each voxel
p=vector();vox=1
for (vox in 1:num_voxels){

  result = t.test(group1[,vox],group2[,vox])
  p[vox] = result[['p.value']]
  
}

#How many significant differences when there is none?
significant_vox = which(p<=thr); print(paste('# voxels with spurious differences',length(significant_vox)))

#Find number of contiguous significant voxels
which(diff(significant_vox)==1) 

#CONCLUSIONS:
#There is an expected amount of detected differences "by chance" (i.e. num_voxels * thr)
#However there isn't more than 2 contiguous voxels 
