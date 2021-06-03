#Author: Grethe Hystad
#Email: ghystad@pnw.edu or grethehystad2@gmail.com

#Core Cluster Identification, clustering method: mixture of normal distributions with nine clusters
#Based on paper by Henelius, A; Puolamäki, K; Boström, H; and Papapetrou, P (2016)
#Clustering with Confidence: Finding Clusters with Statistical Guarantees.
#arXiv: https://arxiv.org/abs/1612.08714 
#################################################################################

library(mclust)        #for clustering with mixture of normal distributions
library(dplyr)         #for data manipulation in data frames  
library(ggplot2)       #for creating graphics
library(igraph)        #for network graphs
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

ncluster=9 #number of cluster

set.seed(1)

# The function, clusterCorrespondence2, is used to solve the cluster label correspondence problem by 
# applying the Hungarian algorithm. The clusterCorrespondence2 function is based on modification of code
# given by Roberto Rösler (Matching clustering solutions using the ‘Hungarian method’, Data*Science+R, November 19, 2012)
#https://www.r-bloggers.com/2012/11/matching-clustering-solutions-using-the-hungarian-method/
clusterCorrespondence2 = function(clusteringA, clusteringB){    
	unA = unique(clusteringA)  # Clusters in A
    	unB = unique(clusteringB)  # Clusters in B
    	neA = length(clusteringA)  # Number of elements in A
    	neB = length(clusteringB)  # Number of elements in B
		if (neA != neB) {
        	stop("number of elements do not match")
     		}
	
    	neC = length(unA) #Number of clusters in A
    	index = c(1:neA)  #Indexed by the elements in A
	neD = length(unB) #Number of clusters in B

    	# Creating the "assignment" matrix
    	assign.Matrix = matrix(numeric(neC*neD), nrow = neC)
    	for (i in 1:neC){
       	 indexA = index[clusteringA == i]
			for (j in 1: neD){
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
overallbase=Mclust(data3[c(9:12)],model="VVV",G=7) 
overallbase.clustering=overallbase$classification

#Clustering with a mixture of normal distributions
base=Mclust(data3[c(9:12)], model="VVV",G=ncluster)
clusternew=base$classification

# Match the cluster labels from clusternew
Correspondence2=clusterCorrespondence2(overallbase.clustering, clusternew)

# Cluster labels from the overall.base clustering matched to clusternew 
Correspondence2

#Relabel clusters in clusternew
index1=which(clusternew==1)
index2=which(clusternew==2)
index3=which(clusternew==3)
index4=which(clusternew==4)
index5=which(clusternew==5)
index6=which(clusternew==6)
index7=which(clusternew==7)
index8=which(clusternew==8)
index9=which(clusternew==9)
clusternew[index4]=3
clusternew[index7]=6
clusternew[index8]=7
clusternew[index9]=4
clusternew[index3]=8
clusternew[index6]=9


#Create the consensus matrix
L=500 # Number of iterations

#Create the connectivity matrix
M=matrix(rep(0,K*K),nrow=K)

#Create the indicator matrix
IN=matrix(rep(0,K*K),nrow=K)

for (j in 1:L){
	newdata=data3%>% sample_n(as.integer(K),replace=T)
	newdata=unique(newdata)
	mm=Mclust(newdata[c(9:12)], model="VVV",G=ncluster)
	cl=mm$classification
	
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

Core.cluster=numeric(ncluster)  #Elements in each core cluster
Core.cluster2=numeric(ncluster) #Number of elements in each core cluster
for(k in 1:ncluster){
	cl1.index=which(clusternew==k)
	cons=cons.Matrix[cl1.index,cl1.index]
      K1=length(cl1.index)
	cons.A=matrix(rep(0,K1*K1),nrow=K1)

	for(i in 1: K1){
		for (j in i:K1){
			if(cons[i,j]>=0.90) {
			cons.A[i,j]=1
			cons.A[j,i]=1
		}
			else{
				cons.A[i,j]=0
				cons.A[j,i]=0
			}
		}
	}

	g=graph_from_adjacency_matrix(cons.A, mode = c("undirected"))
	lg=largest_cliques(g)[1]          #Choose the first largest maximal clique
	Core.cluster[k]=lg                #Elements in largest maximal clique
	Core.cluster2[k]=clique.number(g) #Number of elements in the largest maximal clique

}
sum(Core.cluster2)       #Total number of observations in core clusters
sum(Core.cluster2)/K     #Proportion of observations in core clusters

cluster.index=numeric(ncluster)  #Elements in each cluster
cluster.index2=numeric(ncluster) #Number of elements in each cluster
for (k in 1:ncluster){
	cluster.index[k]=list(which(clusternew==k))
	cluster.index2[k]=length(which(clusternew==k))
	}

#Proportion of observations in core cluster for each cluster
Core.cluster2/cluster.index2


Core.cluster3=numeric(ncluster)       # Observations in core clusters
for (k in 1:ncluster){
	cl1.clust = which(clusternew==k)
	Core.cluster3[k]=list(cl1.clust[unlist(Core.cluster[k])])
	}

M.core=unlist(Core.cluster3)  

# Number of different types of grains
 lZ=length(which(Type=="Z"))
 lY=length(which(Type=="Y"))
 lX=length(which(Type=="X"))
 lN=length(which(Type=="N"))
 lM=length(which(Type=="M"))
 lAB=length(which(Type=="AB"))

 
#Proportion of types of grains in core clusters
length(which(Type[M.core]=="Z"))/lZ
length(which(Type[M.core]=="Y"))/lY
length(which(Type[M.core]=="X"))/lX
length(which(Type[M.core]=="N"))/lN
length(which(Type[M.core]=="M"))/lM
length(which(Type[M.core]=="AB"))/lAB


#Core clusters added to data frame
data.core=data3 %>% slice(M.core)
un1=unlist(Core.cluster3[1])
un2=unlist(Core.cluster3[2])
un3=unlist(Core.cluster3[3])
un4=unlist(Core.cluster3[4])
un5=unlist(Core.cluster3[5])
un6=unlist(Core.cluster3[6])
un7=unlist(Core.cluster3[7])
un8=unlist(Core.cluster3[8])
un9=unlist(Core.cluster3[9])

clustercore=c(rep(1,length(un1)),rep(2,length(un2)), rep(3,length(un3)),rep(4,length(un4)),rep(5,length(un5)),rep(6,length(un6)),rep(7,length(un7)),rep(8,length(un8)),rep(9,length(un9)))
data.core= data.core%>% mutate(clustercore=clustercore)

#Figures

#Plot core clusters in 12C/13C and 14N/15N plot
pdf("C:/Users/ghystad/Documents/Clustering/SiC_clustering/Publication_MNRAS/Figures_ready/fig8g.pdf",family="Helvetica")
coordProj(data.core[c(5:6)], classification=clustercore,col=c("#A6CEE3", "#E31A1C", "#6A3D9A","#FF7F00", "#B15928","#1F78B4","#33A02C", "#CAB2D6","#FB9A99"),symbols=c(20,17,4,9,18,0,1,3,5), log ="xy", xaxt='n',  yaxt='n', ann=F, xlim=c(1,10^4),ylim=c(1,20000))
legend("topright", c("1","2","3","4","5","6","7", "8", "9"),col=c("#A6CEE3", "#E31A1C", "#6A3D9A","#FF7F00", "#B15928","#1F78B4","#33A02C", "#CAB2D6","#FB9A99"),pch=c(20,17,4,9,18,0,1,3,5),cex=1.2, pt.cex=1)
axis(1, at = c(10^0, 10^1, 10^2, 10^3, 10^4),
     c(expression(10^0),expression(10^1),expression(10^2),
       expression(10^3), expression(10^4)), cex.lab=1.5, cex.axis=1.5)

axis(2, at = c(10^1, 10^2, 10^3, 10^4),
     c(expression(10^1),expression(10^2),
       expression(10^3), expression(10^4)), cex.lab=1.5, cex.axis=1.5)
title(xlab=bquote(""^12*"C/"^13*"C"), line=2.5, cex.lab=1.8)
title(ylab=bquote(""^14*"N/"^15*"N"), line=2.3, cex.lab=1.8)
C12_C13_S = 89 
N14_15_S = 440
N14_15_E = 272 
abline(v=C12_C13_S) 
segments(0.1, N14_15_S, x1 = 5500, y1 = N14_15_S, col = par("fg"), lty = par("lty"), lwd = par("lwd"))
segments(0.1, N14_15_E, x1 = 5500, y1 = N14_15_E, col = par("fg"), lty = 5, lwd = par("lwd"))
text(2700,550,label="Solar", col="black",cex=1.4)
text(2000,330,label="Terrestrial", col="black",cex=1.4)
text(70,1.4,label="Solar", col="black",cex=1.4, srt=90)
text(1.1,16000,label="(g)", col="black",cex=1.8)
dev.off()

#Clusters in 29Si/28Si and 30Si/28Si plot
pdf("C:/Users/ghystad/Documents/Clustering/SiC_clustering/Publication_MNRAS/Figures_ready/fig8h.pdf",family="Helvetica")
coordProj(data.core[c(8:7)], what="classification",classification=clustercore,col=c("#A6CEE3", "#E31A1C", "#6A3D9A","#FF7F00", "#B15928","#1F78B4","#33A02C", "#CAB2D6","#FB9A99"),symbols=c(20,17,4,9,18,0,1,3,5), xaxt='n',  yaxt='n', ann=F, xlim=c(-800,900),ylim=c(-700,600))
legend("topright", c("1","2","3","4","5","6","7", "8", "9"),col=c("#A6CEE3", "#E31A1C", "#6A3D9A","#FF7F00", "#B15928","#1F78B4","#33A02C", "#CAB2D6","#FB9A99"),pch=c(20,17,4,9,18,0,1,3,5),cex=1.2, pt.cex=1)
axis(1, at = c(-800,-600,-400,-200, 0, 200, 400,600,800), cex.lab=1.5, cex.axis=1.5)
axis(2, at = c(-800,-600,-400,-200,0,200,400,600,800),cex.lab=1.5, cex.axis=1.5)
title(xlab=expression(paste(delta, "("^{30}, "Si/"^{28}, "Si)(\u2030)") ),line=2.8, cex.lab=1.7)
title(ylab=expression(paste(delta, "("^{29}, "Si/"^{28}, "Si)(\u2030)") ),line=2.1, cex.lab=1.7)
Si29_28_0 = 0.0506331
Si30_28_0 = 0.0334744
abline(v=Si30_28_0) 
segments(-890, Si29_28_0, x1 = 950, y1 = Si29_28_0, col = par("fg"), lty = par("lty"), lwd = par("lwd"))
text(-43,-660,label="Solar", col="black",cex=1.4,srt=90)
text(650,33,label="Solar", col="black",cex=1.4)
text(-780,570,label="(h)", col="black",cex=1.8)
dev.off()

#Bar plot of core clusters
Type[Type=="M"]="MS"
Type2=Type[M.core]
pdf("C:/Users/ghystad/Documents/Clustering/SiC_clustering/Publication_MNRAS/Figures_ready/fig8i.pdf",family="Helvetica", width=10.5, height=10.5)
counts = table(Type2,clustercore)
data.freq=as.data.frame(counts)
ggplot(data=data.freq, aes(x=clustercore, y=Freq, fill=Type2)) +
  geom_bar(stat="identity", width=0.8) +
scale_y_continuous(breaks = seq(0,575, 100), limits = c(0, 575)) +
labs(title="", 
         x="Cluster", y = "Frequency")+
scale_fill_manual(values=c("#A6CEE3", "#CAB2D6", "#6A3D9A", "#1F78B4", "#B15928", "#33A02C"))+
coord_cartesian(ylim=c(0,575))+
geom_text(x=0.8, y=570, label="(i)",size=11)+
labs(x="Cluster",fill="Grain Type")+
theme_classic()+
theme(aspect.ratio = 1,
        axis.text = element_text(size=30),
        axis.title = element_text(size=30),
        legend.title=element_text(size=28),
        legend.text=element_text(size=28),
        panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(),
  	  panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"))
dev.off()













