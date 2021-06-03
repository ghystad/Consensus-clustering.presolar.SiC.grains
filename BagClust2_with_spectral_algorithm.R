#Author: Grethe Hystad
#Email: ghystad@pnw.edu or grethehystad2@gmail.com

#BagClust 2, Clustering method: Spectral Clustering
#BagClust 2 is based on the algorithm provided in the paper:
#Dudoit, S and Fridlyand, J (2003) Bagging to improve the accuracy of a clustering procedure.
#Bioinformatics 19: 1090-1099
#https://academic.oup.com/bioinformatics/article/19/9/1090/284978
################################################################################################################

library(mclust)        #for clustering with mixture of normal distributions 
library(cluster)       #for PAM clustering
library(kknn)          #for spectral clustering 
library(dplyr)         #for data manipulation in data frames  
library(ggplot2)       #for creating graphics
library(RcppHungarian) #for the Hungarian algorithm

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

L=500 # Number of iterations

#Create the connectivity matrix
M=matrix(rep(0,K*K),nrow=K)

#Create the indicator matrix
IN=matrix(rep(0,K*K),nrow=K)

for (j in 1:L){
	newdata=data3%>% sample_n(as.integer(K),replace=T)
	newdata=unique(newdata)
	mm=specClust(newdata[c(9:12)], centers=ncluster, nn = 130, method = "random-walk")
	cl=mm$cluster
	
	#Part of the code for the creation of the consensus matrix is based on modification of code by
	#Alessandra Cabassi and Paul D W Kirk, consensus-cluster.R 
	#In coca: Cluster-of-Clusters Analysis, 2020
	#https://rdrr.io/cran/coca/src/R/consensus-cluster.R 
	
	ML=matrix(rep(0,K*K),nrow=K)
	Id_matrix=matrix(rep(0,K*K),nrow=K)
      for (i in 1:ncluster){
      	ML[newdata$ID2,newdata$ID2] = ML[newdata$ID2,newdata$ID2]+ crossprod(t(as.numeric(cl == i)))
	}

      Id_matrix[newdata$ID2,newdata$ID2] = 1

	#Number of times each pair of observations occur in the same cluster
	M  = M+ML
	#Number of times each pair of observations occur in the same sample
	IN = IN+Id_matrix
}

#Consensus matrix
cons.Matrix = M/IN 
#replace NA with 0
cons.Matrix=replace(cons.Matrix,is.na(cons.Matrix),0)

#Dissimilarity matrix
D=1-cons.Matrix
#Make sure the diagonal elements are zero
diag(D)=0

# Cluster with PAM
pm=pam(D,ncluster)
table(pm$clustering,Type)

plot(silhouette(pm))
summary(pm)
si=silhouette(pm)
silhouette(pm)[,3]

stable=which(silhouette(pm)[,3]>0.5)  #Observations with silhouette width greater than 0.5
length(stable)/K  #Proportion of observations with silhouette width greater than 0.5

silhouette(pm)[stable,]
unstable=which(silhouette(pm)[,3]<=0.5)
silhouette(pm)[unstable,]

# Match the labels from the reference clustering, overallbase.clustering, to have the same
#labels across different clustering methods
Correspondence=clusterCorrespondence(pm$clustering, overallbase.clustering)
clusternew = Correspondence[pm$clustering]

table(clusternew,Type) 

