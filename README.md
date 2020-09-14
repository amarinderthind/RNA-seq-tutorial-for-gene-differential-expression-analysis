# RNASeq tutorial for gene differential expression analysis
This tutorial is created for educational purposes. 

[DOI 10.13140/RG.2.2.10955.62242/1]( http://doi.org/10.13140/RG.2.2.10955.62242/1)

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
#  [RNA-Seq-DGE.rmd](https://github.com/amarinderthind/RNA-seq-tutorial-for-gene-differential-expression-analysis/blob/master/RNA-Seq-DGE.rmd) used to create output of the script shown in [the PDF file here](https://github.com/amarinderthind/RNA-seq-tutorial-for-gene-differential-expression-analysis/blob/master/RNA-Seq-DGE.pdf).

No significant enrichment found from the demo example, so enrichments plots are empty or commented. 

# Note:
Enrichment analysis is based on  Deseq2. However, users may be interested in considering only those genes that are commonly differentially expressed between DEseq2 and EdgeR.  
If data obtained by different batch processing please consider ~batch (batch effect in the design matrix). 

# Contact
In case you have any query please feel free to contact thind.amarinder@gmail.com for any other queries.
Amarinder Singh Thind created this repository that integrates many Bioconductor packages and All rights reserved.
