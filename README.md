# RNA-seq-tutorial-for-gene-differential-expression-analysis
This tutorial is created for educational purposes. 

# About the RNA-Seq analysis
The R script performs several steps in RNAseq gene differential expression analysis, including filtering, preprocessing, visualization, clustering, and Enrichment. For the analysis, several R Bioconductor packages are required to be installed (Installation commands are provided in the script. However, users can also refer to the Bioconductor website for detailed instructions). 

# Required data files
You should have a raw count and annotation/metadata file for running this analysis. Raw count files are usually obtained from tools such as featureCount, Rsem etc.

# Bioconductor packages to be installed
 DESeq2
 
 edgeR
 
 biomaRt (Very useful for gene filtering and annotations)
 
 PCAtools (PCA detailed analysis)
 
 ReactomePA (enrichment analysis)

# [RNA-Seq-DGE.R](https://github.com/amarinderthind/RNA-seq-tutorial-for-gene-differential-expression-analysis/blob/master/RNA-Seq-DGE.R) is the R script.

# Note:
PCA and Enrichment analysis is based on  Deseq2. However, users may be interested in considering only those genes that are commonly differentially expressed between DEseq2 and EdgeR.  

# Contact
In case you have any query please feel free to contact thind.amarinder@gmail.com for any other queries.
Amarinder Singh Thind created this repository that integrates many Bioconductor packages and All rights reserved.
