# SWIPE-NMF

This readme file accompanies our manuscript "Integrative construction of regulatory region networks in 127 human reference epigenomes by matrix factorization "( https://www.biorxiv.org/content/10.1101/217588v1.abstract)

Despite large experimental and computational efforts aiming to dissect the mechanisms underlying disease risk, mapping cis-regulatory elements to target genes remains a challenge. Here, we introduce a matrix factorization framework to integrate physical and functional interaction data of genomic segments. The framework was used to predict a regulatory network of chromatin interaction edges linking more than 20,000 promoters and 1.8 million enhancers across 127 human reference epigenomes, including edges that are present in any of the input datasets. Our network integrates functional evidence of correlated activity patterns from epigenomic data and physical evidence of chromatin interactions. An important contribution of this work is the representation of heterogeneous data with different qualities as networks. We show that the unbiased integration of independent data sources suggestive of regulatory interactions produces meaningful associations supported by existing functional and physical evidence, correlating with expected independent biological features.



AllFunctions11.R file contains all the functions needed to conduct SWIPE_NMF

GenerateInput file convert raw genomic or functional data into appropriate matrices that can be taken into SWIPE_NMF

NMFFunction is the matrix factorization function used to integrate different data types in a "intermediate integration manner"

RunNMFCluster is the code run the integration method, 
