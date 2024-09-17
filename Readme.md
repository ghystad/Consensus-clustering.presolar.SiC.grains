# The R-code is created for the paper

## Evaluation of the classification of presolar silicon carbide grains using consensus clustering with resampling methods: an assessment of the confidence of grain assignments
Grethe Hystad^1, Asmaa Boujibar^2, Nan Liu^3, Larry R. Nittler^2, and Robert M. Hazen^2

*1.Department of Mathematics and Statistics, Purdue University Northwest, 2200 169th Street Hammond, IN 46323-2094, U.S.A.*
*2.Earth and Planets Laboratory, Carnegie Institution for Science, 5251 Broad Branch Rd, NW, Washington, DC 20015, U.S.A.*
*3.Department of Physics, Washington University in St. Louis, One Brookings Drive, St. Louis, MO 63130, U.S.A.*  

The paper is using several consensus clustering methods:

BagClust 1 and BagClust 2 are based on the algorithm provided in the paper
Dudoit, S and Fridlyand, J (2003) Bagging to improve the accuracy of a clustering procedure.
Bioinformatics 19: 1090-1099
https://academic.oup.com/bioinformatics/article/19/9/1090/284978

Core Cluster Identification method is based on the paper by
Henelius, A; Puolamäki, K; Boström, H; and Papapetrou, P (2016)
Clustering with Confidence: Finding Clusters with Statistical Guarantees.
arXiv: https://arxiv.org/abs/1612.08714

Jaccard Similarity Coefficient for testing stability of clusters is based on the paper by
Hennig, C (2007) Cluster-wise assessment of cluster stability.
Computational Statistics & Data Analysis 52:258-271
https://www.sciencedirect.com/science/article/abs/pii/S0167947306004622?via%3Dihub

The paper is using spectral clustering, clustering with a mixture of normal distributions, and clustering with a mixture of t-distributions.

I here provide the R-code for using spectral clustering with seven clusters and clustering with a mixture of normal distributions using nine clusters.

In order to obtain the results for clustering with a mixture of t-distributions using seven clusters, replace the clustering methods with 

base=teigen(data3[c(9:12)],models="UUUU",Gs=ncluster, scale=FALSE, parallel.cores=TRUE)
base$classification

The data used for the paper are from the database PGD_SiC_2020-01-30
Reference for database: 
Stephan, T.; Bose, M; and Boujibar, A et al., The presolar grain database reloaded -  Silicon Carbide, Lunar and Planetary Science Conference, 2020, Abstract #2140
https://presolar.physics.wustl.edu/presolar-grain-database/

The application is released under GNU GPL version 3 license, noting that the following R-libraries and websites are licensed as follows:

mclust: GPL-2, GPL-3 (expanded from GPL-2)
RcppHungarian: GPL-2, GPL-3 (expanded from GPL-2)
teigen:  GPL ≥ 2
kknn: GPL ≥ 2
cluster: GPL ≥ 2
dplyr: MIT
ggplot2: MIT
R-bloggers: MIT or Apache Licenses

Cite the code: [![DOI](https://zenodo.org/badge/373536377.svg)](https://zenodo.org/doi/10.5281/zenodo.13773052)
