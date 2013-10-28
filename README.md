Taxonomic, Phylogenetic and Functional Diversity in South American Hummingbirds
======

Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW

Coauthors: Boris A. Tinoco, Juan L. Parra, PhD, Leone M. Brown, PhD, Gary Stiles, PhD, Jim A. McGuire, PhD, Catherine H. Graham, PhD

The code is broken into pieces due to the massive computational size of the pairwise comparisons. 

The data is made available in the InputData folder with the exception of the env variables, which are freely avaialble from

www.worldclim.org, and exceed the file size for github

Step 1 Observed Environmental and Spatial Dissimilarity
----------------

Read in assemblage information and extract spatial and environmental information

Code Block: DimDivEnv.R

Step 2 Compute Dimensions of Betadiversity Metrics
------------------------------

While this can be done in pieces, it is greatly advised to run on a cluster, since the code will run 47 serially. On 400 cores, the code takes 7 hrs. 
Thanks to the Xsede supercomputing cluster "Stampede" for their access to allow this size of computation

Code Block: DimDivClusterAll.R

Step 3 Create Figures from the Cluster Output
--------------------------------

FinalData.csv - Contains all the environmental information, spatial information, and observed betadiversity for all 23871 unique comparisons (219*218/2)

FinalDataNull.csv - Contains the combination of betadiversity dimension deliniations (ie. High Taxonomic, High Phylogenetic, Low Trait) for all unique assemblage comparisons

and create the final figures.

Code Block: PostCluster_Figures.R

Step 4 Permutation Tests
--------------------------------------
Bootstrapping of Environmental and Spatial variables for each of the eight combinations of betadiversity dimensions and predictors of betadiversity

Credit to Tony Ives and Matt Helmus for their suggestions of how best to create this null model.

Code Block: IvesCoreTest.R


Appendix 1 of comparison of betadiversity metrics

Code Block: CompareMetrics.R

Step 5

