### RNASeq tutorial for gene differential expression analysis and Funcrional enrichment analysis
#### (Updated on 15 Oct 2022)
This tutorial is created for educational purposes and was presentated on Workshop organised by Dollar education.

Interested in exploring more applications of the RNASeq, read here more https://ro.uow.edu.au/test2021/3578/ 
### About the RNA-Seq analysis
The R script performs several steps in RNAseq gene differential expression analysis, including filtering, preprocessing, visualization, clustering, and Enrichment. For the analysis, several R Bioconductor packages are required to be installed (Installation commands are provided in the script. However, users can also refer to the Bioconductor website for detailed instructions). 

### Required data files
You should have a raw count and annotation/metadata file for running this analysis. Raw count files are usually obtained from tools such as featureCount, Rsem etc.

### Main Bioconductor packages to be installed
 DESeq2, edgeR, biomaRt (Very useful for gene filtering and annotations), PCAtools (PCA detailed analysis), ReactomePA (enrichment analysis)

### [RNA-Seq-DGE.R](https://github.com/amarinderthind/RNA-seq-tutorial-for-gene-differential-expression-analysis/blob/master/RNA-Seq-DGE.R) is the R script.
  [RNA-Seq-DGE.rmd](https://github.com/amarinderthind/RNA-seq-tutorial-for-gene-differential-expression-analysis/blob/master/RNA-Seq-DGE.rmd) used to create output of the script shown in [the PDF file here](https://github.com/amarinderthind/RNA-seq-tutorial-for-gene-differential-expression-analysis/blob/master/RNA-Seq-DGE.pdf).

No significant enrichment found from the demo example, so enrichments plots are empty or commented. 

### Note:
If data obtained by different batch processing please consider ~batch (batch effect in the design matrix) or use CombatSeq as defined in the new script. 

Plots 
![PCA plot showing Batch effect before normalization](https://user-images.githubusercontent.com/45668229/166874766-39a3a488-f97e-44b9-b704-659415aba683.png)
![PCA plot after batch effect correction using CombatSeq](https://user-images.githubusercontent.com/45668229/166874917-03255c28-b586-4a26-9b24-20e5dbcc2299.png)

![top 20 DEG with norm count](https://user-images.githubusercontent.com/45668229/151488929-7f5c2517-935d-472c-96fc-c91e0afe2642.png)

![HeatMap](https://user-images.githubusercontent.com/45668229/151488668-0722347f-6768-47db-8ea1-7fb6d42b2e8c.png)
![Volcano_plot](https://user-images.githubusercontent.com/45668229/151488762-172ce41c-d5d5-46d2-b1b2-977f91db9365.png)

![Enrichment ORA test](https://user-images.githubusercontent.com/45668229/166875037-6c2256c6-86a6-4d20-bbcb-7b22cb47ab09.png)
![Enrichment GSEA](https://user-images.githubusercontent.com/45668229/166875114-89edac5a-b946-43eb-b18e-f02d361b7220.png)



### Reading material or relavant articles
[Explore about different normalization methods here](https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html)

[Other emerging bulk RNASeq applications](https://doi.org/10.1093/bib/bbab259)

[EdgeR](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)

[Deseq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

[PCAtools](https://bioconductor.org/packages/release/bioc/vignettes/PCAtools/inst/doc/PCAtools.html)

[Reactome Pathway Analysis](https://www.bioconductor.org/packages/release/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html)

#### Contact
In case you have any query please feel free to contact thind.amarinder@gmail.com for any other queries.
