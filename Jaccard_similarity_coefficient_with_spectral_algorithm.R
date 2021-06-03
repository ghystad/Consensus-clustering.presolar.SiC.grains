
#Author: Grethe Hystad
#Email: ghystad@pnw.edu or grethehystad2@gmail.com

#Jaccard Similarity Coefficient, Clustering method: Spectral Clustering
#Based on the paper by Hennig, C (2007) Cluster-wise assessment of cluster stability.
#Computational Statistics & Data Analysis 52:258-271
#https://www.sciencedirect.com/science/article/abs/pii/S0167947306004622?via%3Dihub

library(mclust)       #for clustering with mixture of Gaussians
library(kknn)         #for spectral clustering
library(dplyr)        #for data manipulation in data frames  
library(ggplot2)      #for creating graphics
library(RcppHungarian)#for the Hungarian algorithm

###################################################################

#Data preparation

#Read in the data
data=read.delim("presolar_SiC_data.txt",header=TRUE)
attach(data)

data = data%>% mutate(X29Si.28Si = ((d.29Si.28Si./1000)+1)*0.0506331 
,X30Si.28Si = ((d.30Si.28Si./1000)+1)*0.0334744 )

ID=PGD.ID
Type=PGD.Type

colnames(data) = c("ID", "Type", "Meteorite", "Technique", "X12C.13C", "X14N.15N","d.29Si.28Si","d.30Si.28Si", "X29Si.28Si","X30Si.28Si")

#Omit rows with Type "C" and Type "U"
sub=which(Type=="C" | Type=="U")
data2=data[-sub,]

#Omit rows with NA 
data2=na.omit(data2)

#Scale the data
l.data=sapply(data2[,c(5,6,9,10)], log10)
l.data.scaled=scale(l.data, center = TRUE, scale = TRUE)

#Select the original data
data3=data2%>% select(1:8)

#Combine the original data with the scaled data
data3=cbind(data3,l.data.scaled)
colnames(data3) = c("ID", "Type", "Meteorite", "Technique", "X12C.13C", "X14N.15N","d.29Si.28Si","d.30Si.28Si","scaled.X.12C.13C", "scaled.X.14N.15N", "scaled.X.29Si.28Si","scaled.X.30Si.28Si")

Type=data3$Type

#Number of observations
K=dim(data3)[1]
#Add a unique number to each observation
ID2=seq(1, K,1)
data3 = mutate(data3, ID2)
head(data3)

ncluster=7 #number of cluster

set.seed(1)

# The function, clusterCorrespondence, is used to solve the cluster label correspondence problem by 
# applying the Hungarian algorithm. The clusterCorrespondence function is based on modification of code
# given by Roberto Rösler (Matching clustering solutions using the ‘Hungarian method’, Data*Science+R, November 19, 2012)
#https://www.r-bloggers.com/2012/11/matching-clustering-solutions-using-the-hungarian-method/
clusterCorrespondence = function(clusteringA, clusteringB){    
	unA = unique(clusteringA)  # Clusters in A
    	unB = unique(clusteringB)  # Clusters in B
    	neA = length(clusteringA)  # Number of elements in A
    	neB = length(clusteringB)  # Number of elements in B
		if (length(unA) != length(unB)) {
       	 stop("number of clusters do not match")
		}
		if (neA != neB) {
        	stop("number of elements do not match")
     		}

    	neC = length(unA) #Number of clusters in A
    	index = c(1:neA)  #Indexed by the elements in A

    	# Creating the "assignment" matrix
    	assign.Matrix = matrix(numeric(neC*neC), nrow = neC)
    	for (i in 1:neC){
       	 indexA = index[clusteringA == i]
			for (j in 1: neC){
            		indexB = index[clusteringB == j]
            		indexCommon = length(intersect(indexA, indexB))
            		assign.Matrix[i,j] = (length(indexA)- indexCommon)+(length(indexB)- indexCommon)
			}
   	 }


    	#Hungarian algorithm. Assignment of rows to columns 
     	optimalAssign2 = HungarianSolver(assign.Matrix)
      optimalAssign = optimalAssign2$pairs[,2]
     	return(optimalAssign)
}

# Using the labels from the mixture of normal distributions as the reference cluster labels
overallbase=Mclust(data3[c(9:12)],model="VVV",G=ncluster) 
overallbase.clustering=overallbase$classification

#Spectral clustering
base=specClust(data3[c(9:12)], centers=ncluster, nn = 130,method = "random-walk")
clusteringB2=base$cluster


#BagClust 1

B=500 #Number of bootstrap samples

# Matrix to keep track of the observations that occur for each bootstrap sample
IN=array(NA, c(B,K)) 

# Matrix to keep track of cluster label for each observation for each bootstrap sample
clust=array(NA, c(B,K))

#Apply the clustering procedure for B bootstrap samples. 
#Keep only the unique samples
for (i in 1:B){
	newdata=data3%>% sample_n(K,replace=TRUE)
	newdata=unique(newdata)
	sc=specClust(newdata[c(9:12)],centers=ncluster, nn = 130, method = "random-walk")
	clusteringA=sc$cluster
	IN[i,newdata$ID2]= 1
	clusteringB=clusteringB2[newdata$ID2]
	
	# match the cluster labels from clusteringB, which are the original cluster labels
      Correspondence=clusterCorrespondence(clusteringA,clusteringB)
      clusterA = Correspondence[clusteringA]
      clust[i,newdata$ID2]=clusterA
}

#Count the number of times each observation was assigned to each cluster
M=array(NA,c(K,(ncluster)))
for (j in 1:K){
	for(i in 1:ncluster){
		M[j,i]=c(length(which(clust[,j]==i)))
	}
}

#Cluster assignment for each observation based on majority vote
M.vote=numeric(K)   #majority vote
for (l in 1:K){
	maj=which(M[l,]==max(M[l,]))
		if (length(maj) > 1 ){
		M.vote[l]=sample(maj,1)     #In case there is a tie, we randomly sample among the ties 
		} else{
		M.vote[l]=maj
		}
}

# Number of winning cluster assignment for each observation
Winning=numeric(K)                
for (j in 1:K){
	Winning[j]=max(M[j,])
}

# Number of times each observation occurs over the B bootstrap samples
IN2=numeric(K)
for (j in 1:K){
	IN2[j]=sum(IN[,j], na.rm=TRUE)
}

#Cluster vote for each observation
CV=Winning/IN2                          
#Put all NA equal to 0
CV=replace(CV,is.na(CV),0)
 
length(which(CV>=0.90))  #Number of observations with cluster vote greater than or equal to 0.90

# Match the labels from the reference clustering, overallbase.clustering, to have the same
#labeling across different clustering methods
Correspondence=clusterCorrespondence(M.vote, overallbase.clustering)
clusternew = Correspondence[M.vote]


#Jaccard Coefficient
clu=numeric(ncluster)
for (k in 1:ncluster){
	M=100
	gammaboot=numeric(M)
	for (j in 1:M){
		newdata=unique(data3%>% sample_n(K,replace=T))
		cl=specClust(newdata[c(9:12)], centers=ncluster, nn = 130, method = "random-walk")
		jaccard=numeric(ncluster)
		for(i in 1:ncluster){
			nk=data3[which(clusternew==k),]
			Cstar=intersect(nk,newdata)
			delta=newdata[which(cl$cluster==i),]
			d=dim(intersect(Cstar,delta))[1]
			u=dim(union(Cstar,delta))[1]
			jaccard[i] = d/u
			}
			gammaboot[j] = max(jaccard)
		}
	clu[k]=mean(gammaboot[gammaboot>0])

}




