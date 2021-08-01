# RNASeq tutorial for gene differential expression analysis
This tutorial is created for educational purposes. 

[DOI 10.13140/RG.2.2.26508.13443]( http://doi.org/10.13140/RG.2.2.26508.13443)

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
  [RNA-Seq-DGE.rmd](https://github.com/amarinderthind/RNA-seq-tutorial-for-gene-differential-expression-analysis/blob/master/RNA-Seq-DGE.rmd) used to create output of the script shown in [the PDF file here](https://github.com/amarinderthind/RNA-seq-tutorial-for-gene-differential-expression-analysis/blob/master/RNA-Seq-DGE.pdf).

No significant enrichment found from the demo example, so enrichments plots are empty or commented. 

# Note:
If data obtained by different batch processing please consider ~batch (batch effect in the design matrix). 


![alt text](RNA-Seq-DGE.pdf)

# Reading material or relavant articles
[Explore about different normalization methods here](https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html)

[Other emerging bulk RNASeq applications](https://doi.org/10.1093/bib/bbab259)

[EdgeR](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)

[Deseq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

[PCAtools](https://bioconductor.org/packages/release/bioc/vignettes/PCAtools/inst/doc/PCAtools.html)

[Reactome Pathway Analysis](https://www.bioconductor.org/packages/release/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html)

# Contact
In case you have any query please feel free to contact thind.amarinder@gmail.com for any other queries.
Amarinder Singh Thind created this repository that integrates many Bioconductor packages and All rights reserved.
