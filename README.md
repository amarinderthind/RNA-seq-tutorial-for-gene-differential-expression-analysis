### RNASeq tutorial for gene differential expression analysis and Functional enrichment analysis

This tutorial was created for educational purposes and was presented at a Workshop organised by Dollar Education.

Interested in exploring more applications of the RNASeq, Read here https://ro.uow.edu.au/test2021/3578/ 

Want to adjust for tumour purity Check  https://www.nature.com/articles/s41587-022-01440-w

### About the RNA-Seq analysis
The R script performs several steps in RNAseq gene differential expression analysis, including filtering, preprocessing, visualization, clustering, and Enrichment. For the analysis, several R Bioconductor packages are required to be installed (Installation commands are provided in the script. However, users can also refer to the Bioconductor website for detailed instructions). 

### Main Bioconductor packages to be installed
 DESeq2, edgeR, biomaRt (Very useful for gene filtering and annotations), PCAtools (PCA detailed analysis), ReactomePA (enrichment analysis)

### [RNA-Seq-DGE.R](https://github.com/amarinderthind/RNA-seq-tutorial-for-gene-differential-expression-analysis/blob/master/RNA-Seq-DGE.R) is the R script.
  [RNA-Seq-DGE.rmd](https://github.com/amarinderthind/RNA-seq-tutorial-for-gene-differential-expression-analysis/blob/master/RNA-Seq-DGE.rmd) used to create an output of the script shown in [the PDF file here](https://github.com/amarinderthind/RNA-seq-tutorial-for-gene-differential-expression-analysis/blob/master/RNA-Seq-DGE.pdf).

No significant enrichment was found from the demo example, so enrichment plots are empty or commented. 

### Note:
If data is obtained by different batch processing please consider ~batch (batch effect in the design matrix) or use CombatSeq as defined in the new script. 


### Required data files
You should have a raw count and annotation/metadata file for running this analysis. Raw count files are usually obtained from tools such as featureCount, Rsem etc.

```
setwd("/Users/Path") #Path_to_working_directory

rawcount<-read.table ("RawCount_input.csv",header=TRUE,  sep=",",  row.names=1)

## Replace NAs by zero and change the input to the required format
rawcount <- round(rawcount) 
rawcount[is.na(rawcount)] <- 0
```

Loading and filtering Data annotation  file

```
anno <-read.table ("Annotation_of_samples_12_Samples_ALL.csv",header=TRUE,  sep=",", row.names = 1) ##, In this case, we have 3 columns (a) sample (b) Condition (c) batch
#rownames(anno) <- anno$sample  ##add row names as sample name (if not already), because PCA function check row names of anno == col of the data matrix

table(anno$Condition)

library(tidyverse)
library('dplyr') ##HAS COUNT FUNCTION

### In case you want to consider the subset of samples based on some condition (when multiple e.g. >3 )
#anno <- anno %>% 
#  as.data.frame %>%
#  filter(anno$Condition =='Condition_A' |anno$Condition =='Condition_B' | anno$Condition == 'Condition_C' )  %>%
#  arrange(Condition)  	#Arrange rows by padj values 
```
sort anno based on condition ## sorted conditional sample representation in heatmap 

```
anno <- anno %>% 
  as.data.frame %>%
  arrange(Condition) 
```

PCA plot for pre-DE investigation 

```
library(PCAtools)

anno <- anno[match(colnames(rawcount), anno$Sample),] ## reordering anno rows with colnames of rawcount
lograwcount <- as.matrix(log2(rawcount +1))  ## log transformation of rawcount for PCA plot 

 top1000.order <- head(order(matrixStats::rowVars(lograwcount), decreasing = TRUE), 1000) ## taking top 1000 genes having the highest variance selected from all the genes in the input
  p <- PCAtools::pca(mat = lograwcount[top1000.order,], metadata = anno, removeVar = 0.01) ## performing PCA

  biplot(p,                                                       #visualization of PCA plot
       lab = paste0(p$metadata$Sample),
        colby = 'Batch',  #Sample #Batch #Condition #sex
        hline = 0, vline = 0,
        legendPosition = 'right',
         encircle = T )
  
  screeplot(p, axisLabSize = 18, titleLabSize = 22) ## This plot shows how much variation in the data is explained by which PC component.
  pairsplot(p) ## draw various combinations of the PCA plot
```




### Plots 

#### PCA plots for initial QC

 <p align="center">
<img src="https://user-images.githubusercontent.com/45668229/166874766-39a3a488-f97e-44b9-b704-659415aba683.png" width=48% height="400">&nbsp; &nbsp; &nbsp; &nbsp;
<img src="https://user-images.githubusercontent.com/45668229/166874917-03255c28-b586-4a26-9b24-20e5dbcc2299.png" width=48% height="400">
 
</p>

#### Dot plot and HeatMap


<p align="center">
<img src="https://user-images.githubusercontent.com/45668229/151488929-7f5c2517-935d-472c-96fc-c91e0afe2642.png" width=45% height="500">&nbsp; &nbsp; &nbsp; &nbsp;
<img src="https://user-images.githubusercontent.com/45668229/151488668-0722347f-6768-47db-8ea1-7fb6d42b2e8c.png" width=45% height="500">
 
<img src="https://user-images.githubusercontent.com/45668229/151488762-172ce41c-d5d5-46d2-b1b2-977f91db9365.png" width=85% height="700">&nbsp; &nbsp; &nbsp; &nbsp;
</p>
### Over enrichment analysis (ORA) and Gene set enrichment analysis (GSEA)

<p align="center">
<img src="https://user-images.githubusercontent.com/45668229/166875037-6c2256c6-86a6-4d20-bbcb-7b22cb47ab09.png" width=58% height="400">
<img src="https://user-images.githubusercontent.com/45668229/166875114-89edac5a-b946-43eb-b18e-f02d361b7220.png" width=40% height="400">
</p>

### Reading material or relevant articles
[Explore about different normalization methods here](https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html)

[Other emerging bulk RNASeq applications](https://doi.org/10.1093/bib/bbab259)

[EdgeR](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)

[Deseq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

[PCAtools](https://bioconductor.org/packages/release/bioc/vignettes/PCAtools/inst/doc/PCAtools.html)

[Reactome Pathway Analysis](https://www.bioconductor.org/packages/release/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html)

#### Contact
In case you have any queries please feel free to contact thind.amarinder@gmail.com for any other queries.
