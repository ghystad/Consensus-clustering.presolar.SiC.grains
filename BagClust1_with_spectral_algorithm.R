#Author: Grethe Hystad
#Email: ghystad@pnw.edu or grethehystad2@gmail.com

#BagClust 1, Clustering method: Spectral Clustering
#BagClust 1 is based on the algorithm provided in the paper:
#Dudoit, S and Fridlyand, J (2003) Bagging to improve the accuracy of a clustering procedure.
#Bioinformatics 19: 1090-1099
#https://academic.oup.com/bioinformatics/article/19/9/1090/284978


################################################################################################################

library(mclust)       #for clustering with mixture of normal distributions
library(kknn)         #for spectral clustering
library(dplyr)        #for data manipulation in data frames  
library(ggplot2)      #for creating graphics
library(RcppHungarian)#for the Hungarian algorithm


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

table(clusternew,Type) 

#Proportion of observations with CV greater than or equal to 0.90 for each cluster
Stability.index=numeric(ncluster)
for(k in 1:ncluster){
	Stability.index[k] = length((which(clusternew[which(CV>=0.90)]==k)))/length(which((clusternew==k)))
}
clusters=1:ncluster
cbind(clusters,Stability.index)

data3$majority.vote=clusternew  # added to data3
data3$cluster.vote=CV

#Observations with cluster vote greater than or equal to 0.90
d.stabil=data3%>%filter(CV>=0.90)

Unique.Type=sort(unique(Type))

#Proportion of observations with cluster vote greater than or equal to 0.90 for each grain type (6 different grain types)
Stability.Type=numeric(6)
for (i in 1:6){
	Stability.Type[i]=length(which(d.stabil[c(2)]==Unique.Type[i]))/length(which(data3[,2]==Unique.Type[i]))
}
cbind(Unique.Type,Stability.Type)

#Figures

#Clusters in 12C/13C and 14N/15N plot
pdf("C:/Users/ghystad/Documents/Clustering/SiC_clustering/Publication_MNRAS/Figures_ready/fig5a.pdf",family="Helvetica")
coordProj(data3[c(5:6)], what="classification",classification=clusternew,col=c("#A6CEE3", "#E31A1C", "#6A3D9A","#FF7F00", "#B15928","#1F78B4","#33A02C"),symbols=c(20,17,4,9,18,0,1), log ="xy", xaxt='n',  yaxt='n', ann=F)
legend("topright", c("1","2","3","4","5","6","7"),col=c("#A6CEE3", "#E31A1C", "#6A3D9A","#FF7F00", "#B15928","#1F78B4","#33A02C"),pch=c(20,17,4,9,18,0,1),cex=1.2, pt.cex=1)
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
segments(1, N14_15_S, x1 = 15000, y1 = N14_15_S, col = par("fg"), lty = par("lty"), lwd = par("lwd"))
segments(1, N14_15_E, x1 = 15000, y1 = N14_15_E, col = par("fg"), lty = 5, lwd = par("lwd"))
text(7300,550,label="Solar", col="black",cex=1.4)
text(5000,330,label="Terrestrial", col="black",cex=1.4)
text(70,5,label="Solar", col="black",cex=1.4, srt=90)
text(1.5,16000,label="(a)", col="black",cex=1.8)
dev.off()

#Clusters in 29Si/28Si and 30Si/28Si plot
pdf("C:/Users/ghystad/Documents/Clustering/SiC_clustering/Publication_MNRAS/Figures_ready/fig5b.pdf",family="Helvetica")
coordProj(data3[c(8:7)], what="classification",classification=clusternew,col=c("#A6CEE3", "#E31A1C", "#6A3D9A","#FF7F00", "#B15928","#1F78B4","#33A02C"),symbols=c(20,17,4,9,18,0,1), xaxt='n',  yaxt='n', ann=F, xlim=c(-800,900),ylim=c(-700,600))
legend("topright", c("1","2","3","4","5","6","7"),col=c("#A6CEE3", "#E31A1C", "#6A3D9A","#FF7F00", "#B15928","#1F78B4","#33A02C"),pch=c(20,17,4,9,18,0,1),cex=1.2, pt.cex=1)
axis(1, at = c(-800,-600,-400,-200, 0, 200, 400,600,800), cex.lab=1.5, cex.axis=1.5)
axis(2, at = c(-800,-600,-400,-200,0,200,400,600,800),cex.lab=1.5, cex.axis=1.5)
title(xlab=expression(paste(delta, "("^{30}, "Si/"^{28}, "Si)(\u2030)") ),line=2.8, cex.lab=1.7)
title(ylab=expression(paste(delta, "("^{29}, "Si/"^{28}, "Si)(\u2030)") ),line=2.1, cex.lab=1.7)
Si29_28_0 = 0.0506331
Si30_28_0 = 0.0334744
abline(v=Si30_28_0) 
segments(-890, Si29_28_0, x1 = 970, y1 = Si29_28_0, col = par("fg"), lty = par("lty"), lwd = par("lwd"))
text(-43,-660,label="Solar", col="black",cex=1.4,srt=90)
text(860,33,label="Solar", col="black",cex=1.4)
text(-780,570,label="(b)", col="black",cex=1.8)
dev.off()

#Bar plot
Type[Type=="M"]="MS"
pdf("C:/Users/ghystad/Documents/Clustering/SiC_clustering/Publication_MNRAS/Figures_ready/fig5c.pdf",family="Helvetica", width=10.5, height=10.5)
counts = table(Type,clusternew)
data.freq=as.data.frame(counts)
ggplot(data=data.freq, aes(x=clusternew, y=Freq, fill=Type)) +
  geom_bar(stat="identity", width=0.8) +
scale_y_continuous(breaks = seq(0,575, 100), limits = c(0, 575)) +
labs(title="", 
         x="Cluster", y = "Frequency")+
scale_fill_manual(values=c("#A6CEE3", "#CAB2D6", "#6A3D9A", "#1F78B4", "#B15928", "#33A02C"))+
coord_cartesian(ylim=c(0,575))+
geom_text(x=0.8, y=570, label="(c)",size=11)+
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
                                                              
#Stable clusters, in 12C/13C and 14N/15N plot
pdf("C:/Users/ghystad/Documents/Clustering/SiC_clustering/Publication_MNRAS/Figures_ready/fig5d.pdf",family="Helvetica")
coordProj(d.stabil[c(5:6)], what="classification",classification=clusternew[which(CV>=0.90)],col=c("#A6CEE3", "#E31A1C", "#6A3D9A","#FF7F00", "#B15928","#1F78B4","#33A02C"),symbols=c(20,17,4,9,18,0,1), log ="xy", xaxt='n',  yaxt='n', ann=F,xlim=c(1,10^4),ylim=c(1,20000))
legend("topright", c("1","2","3","4","5","6","7"),col=c("#A6CEE3", "#E31A1C", "#6A3D9A","#FF7F00", "#B15928","#1F78B4","#33A02C"),pch=c(20,17,4,9,18,0,1),cex=1.2, pt.cex=1)
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
segments(0.1, N14_15_S, x1 = 15000, y1 = N14_15_S, col = par("fg"), lty = par("lty"), lwd = par("lwd"))
segments(0.1, N14_15_E, x1 = 15000, y1 = N14_15_E, col = par("fg"), lty = 5, lwd = par("lwd"))
text(3100,550,label="Solar", col="black",cex=1.4)
text(2000,330,label="Terrestrial", col="black",cex=1.4)
text(70,1.4,label="Solar", col="black",cex=1.4, srt=90)
text(1.1,16000,label="(d)", col="black",cex=1.8)
dev.off()

                                                                
#Stable clusters in 29Si/28Si and 30Si/28Si plot
pdf("C:/Users/ghystad/Documents/Clustering/SiC_clustering/Publication_MNRAS/Figures_ready/fig5e.pdf",family="Helvetica")
coordProj(d.stabil[c(8:7)], what="classification",classification=clusternew[which(CV>=0.90)],col=c("#A6CEE3", "#E31A1C", "#6A3D9A","#FF7F00", "#B15928","#1F78B4","#33A02C"),symbols=c(20,17,4,9,18,0,1), xaxt='n',  yaxt='n', ann=F, xlim=c(-800,900),ylim=c(-700,600))
legend("topright", c("1","2","3","4","5","6","7"),col=c("#A6CEE3", "#E31A1C", "#6A3D9A","#FF7F00", "#B15928","#1F78B4","#33A02C"),pch=c(20,17,4,9,18,0,1),cex=1.2, pt.cex=1)
axis(1, at = c(-800,-600,-400,-200, 0, 200, 400,600,800), cex.lab=1.5, cex.axis=1.5)
axis(2, at = c(-800,-600,-400,-200,0,200,400,600,800),cex.lab=1.5, cex.axis=1.5)
title(xlab=expression(paste(delta, "("^{30}, "Si/"^{28}, "Si)(\u2030)") ),line=2.8, cex.lab=1.7)
title(ylab=expression(paste(delta, "("^{29}, "Si/"^{28}, "Si)(\u2030)") ),line=2.1, cex.lab=1.7)
Si29_28_0 = 0.0506331
Si30_28_0 = 0.0334744
abline(v=Si30_28_0) 
segments(-890, Si29_28_0, x1 = 970, y1 = Si29_28_0, col = par("fg"), lty = par("lty"), lwd = par("lwd"))
text(-43,-660,label="Solar", col="black",cex=1.4,srt=90)
text(860,33,label="Solar", col="black",cex=1.4)
text(-780,570,label="(e)", col="black",cex=1.8)
dev.off()



#Bar plot for stable clusters
pdf("C:/Users/ghystad/Documents/Clustering/SiC_clustering/Publication_MNRAS/Figures_ready/fig5f.pdf",family="Helvetica", width=10.5, height=10.5)
Type2=Type[which(CV>=0.90)]
clusternew2=clusternew[which(CV>=0.90)]
counts = table(Type2,clusternew2)
data.freq=as.data.frame(counts)
ggplot(data=data.freq, aes(x=clusternew2, y=Freq, fill=Type2)) +
  geom_bar(stat="identity", width=0.8) +
scale_y_continuous(breaks = seq(0,575, 100), limits = c(0, 575)) +
labs(title="", 
         x="Cluster", y = "Frequency")+
scale_fill_manual(values=c("#A6CEE3", "#CAB2D6", "#6A3D9A", "#1F78B4", "#B15928", "#33A02C"))+
coord_cartesian(ylim=c(0,575))+
geom_text(x=0.8, y=570, label="(f)",size=11)+
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


#Zoomed Clusters in 12C/13C and 14N/15N plot
pdf("C:/Users/ghystad/Documents/Clustering/SiC_clustering/Publication_MNRAS/Figures_ready/fig9a.pdf",family="Helvetica")
coordProj(data3[c(5:6)], what="classification",classification=clusternew,col=c("#A6CEE3", "#E31A1C", "#6A3D9A","#FF7F00", "#B15928","#1F78B4","#33A02C"),symbols=c(20,17,4,9,18,0,1), log ="xy", xaxt='n',  yaxt='n', ann=F,xlim=c(20, 400),ylim=c(100,20000))
legend("topright", c("1","2","3","4","5","6","7"),col=c("#A6CEE3", "#E31A1C", "#6A3D9A","#FF7F00", "#B15928","#1F78B4","#33A02C"),pch=c(20,17,4,9,18,0,1),cex=1.2, pt.cex=1)
axis(1, at = c(1, 20, 100, 400), cex.lab=1.5, cex.axis=1.5)
axis(2, at = c(10^1, 10^2, 10^3, 10^4),
     c(expression(10^1),expression(10^2),
       expression(10^3), expression(10^4)), cex.lab=1.5, cex.axis=1.5)
title(xlab=bquote(""^12*"C/"^13*"C"), line=2.5, cex.lab=1.8)
title(ylab=bquote(""^14*"N/"^15*"N"), line=2.3, cex.lab=1.8)
C12_C13_S = 89 
N14_15_S = 440
N14_15_E = 272 
abline(v=C12_C13_S) 
segments(1, N14_15_S, x1 = 4000, y1 = N14_15_S, col = par("fg"), lty = par("lty"), lwd = par("lwd"))
segments(1, N14_15_E, x1 = 4000, y1 = N14_15_E, col = par("fg"), lty = 5, lwd = par("lwd"))
text(350,500,label="Solar", col="black",cex=1.4)
text(310,320,label="Terrestrial", col="black",cex=1.4)
text(82,17000,label="Solar", col="black",cex=1.4, srt=90)
text(20,18500,label="(a)", col="black",cex=1.8)
dev.off()



#Zoomed Clusters in 29Si/28Si and 30Si/28Si plot
pdf("C:/Users/ghystad/Documents/Clustering/SiC_clustering/Publication_MNRAS/Figures_ready/fig9b.pdf",family="Helvetica")
coordProj(data3[c(8:7)], what="classification",classification=clusternew,col=c("#A6CEE3", "#E31A1C", "#6A3D9A","#FF7F00", "#B15928","#1F78B4","#33A02C"),symbols=c(20,17,4,9,18,0,1), xaxt='n',  yaxt='n', ann=F, xlim=c(-100,280),ylim=c(-100,200))
legend("topright", c("1","2","3","4","5","6","7"),col=c("#A6CEE3", "#E31A1C", "#6A3D9A","#FF7F00", "#B15928","#1F78B4","#33A02C"),pch=c(20,17,4,9,18,0,1),cex=1.2, pt.cex=1)
axis(1, at = c(-200,-100, 0, 100, 200,300), cex.lab=1.5, cex.axis=1.5)
axis(2, at = c(-300,-200,-100,0,100,200,300),cex.lab=1.5, cex.axis=1.5)
title(xlab=expression(paste(delta, "("^{30}, "Si/"^{28}, "Si)(\u2030)") ),line=2.8, cex.lab=1.7)
title(ylab=expression(paste(delta, "("^{29}, "Si/"^{28}, "Si)(\u2030)") ),line=2.1, cex.lab=1.7)
Si29_28_0 = 0.0506331
Si30_28_0 = 0.0334744
abline(v=Si30_28_0) 
segments(-860, Si29_28_0, x1 = 750, y1 = Si29_28_0, col = par("fg"), lty = par("lty"), lwd = par("lwd"))
text(7,187,label="Solar", col="black",cex=1.4,srt=90)
text(265,8,label="Solar", col="black",cex=1.4)
text(-100,195,label="(b)", col="black",cex=1.8)
dev.off()


