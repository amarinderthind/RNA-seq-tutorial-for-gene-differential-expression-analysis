### RNAseq using Deseq2 and Functional enrichment Analysis ####
### Date : 18-19 April, 2022

##### Install packages, if not done before 

# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("biomaRt")
# BiocManager::install('PCAtools')
# BiocManager::install('EnhancedVolcano')

###################### load the raw count matrix #######################

setwd("/Users/athind/Dropbox/RNAseq_using_DEseq2-april16/") #Path_to_working_directory

rawcount<-read.table ("RawCount_input.csv",header=TRUE,  sep=",",  row.names=1)

## Replace NAs by zero and changing the input to required format
rawcount <- round(rawcount) 
rawcount[is.na(rawcount)] <- 0


## Discard genes, which are expressed in less than 20% of all samples (considering that we have 2 conditions in total )
##  % selection is based on no. of conditions  ## more condition means less %
keep <- rowSums(rawcount > 0) >= round(ncol(rawcount)*.20)  
rawcount <- rawcount[keep,]


###################### Data annotation  #################################

anno <-read.table ("Annotation_of_samples_12_Samples_ALL.csv",header=TRUE,  sep=",", row.names = 1) ##In this case we have 3 coulmns (a) sample (b) Condition (c) batch
#rownames(anno) <- anno$sample  ##add rownames as sample name (if not already), because pca function check rownames of anno == col of data matrix

table(anno$Condition)

library(tidyverse)
library('dplyr') ##HAS COUNT FUNCTION

### incase want to consider subset of samples based on some condition (when multiple e.g. >3 )
#anno <- anno %>% 
#  as.data.frame %>%
#  filter(anno$Condition =='Condition_A' |anno$Condition =='Condition_B' | anno$Condition == 'Condition_C' )  %>%
#  arrange(Condition)  	#Arrange rows by padj values 


## sort anno based on condition ## good representation in heatmap 
anno <- anno %>% 
  as.data.frame %>%
  arrange(Condition) 



##############################################################
############ PCA plot for pre DE investigation ##############


library(PCAtools)

anno <- anno[match(colnames(rawcount), anno$Sample),] ## reordering anno rows with colnames of rawcount
lograwcount <- as.matrix(log2(rawcount +1))  ## log transformation of rawcount for PCA plot 

 top1000.order <- head(order(matrixStats::rowVars(lograwcount), decreasing = TRUE), 1000) ## taking top 1000 genes having highest variance selected from all the genes in the input
  p <- PCAtools::pca(mat = lograwcount[top1000.order,], metadata = anno, removeVar = 0.01) ## performing PCA

  biplot(p,                                                       #visualization of PCA plot
       lab = paste0(p$metadata$Sample),
        colby = 'Batch',  #Sample #Batch #Condition #sex
        hline = 0, vline = 0,
        legendPosition = 'right',
         encircle = T )
  
  screeplot(p, axisLabSize = 18, titleLabSize = 22) ## this plot shows how much variation in the data is explained by which PC component.
  pairsplot(p) ## draw various combinations of the PCA plot

##############################################################
################# Lets check combat normalization ############
############## SVA #####################
  
  #BiocManager::install("sva")
  
  library('sva')
  rawcount <- as.matrix(rawcount)
  adjusted_counts <- ComBat_seq(rawcount, batch=anno$Batch, group=anno$Condition) ##In ComBat-seq, user may specify biological covariates, whose signals will be preserved in the adjusted data. I
  
  nor_set <- as.matrix(log2(adjusted_counts+1)) ## log transformation of adjusted count
  top1000.order <- head(order(matrixStats::rowVars(nor_set), decreasing = TRUE), 1000)
  pp <- PCAtools::pca(mat =nor_set[top1000.order,] , metadata = anno, removeVar = 0.01)
  biplot(pp,
          lab = paste0(p$metadata$Sample),
          #colby = 'Batch',   #Batch_log', #Condition
          colby = 'Condition',
          hline = 0, vline = 0,
          legendPosition = 'right',encircle = T)

   

  ##### Do we suppose to remove any default sample/s  #########
   
   ### subset raw and conditional data for defined pairs
   ##### Removing sample number 7 ########## 
   
    
   anno <- anno[!(anno$Sample == 'sample_7' | anno$Sample == 'sample_8'),]
   rawcount <- as.data.frame(rawcount)
   rawcount <- rawcount[,names(rawcount) %in% anno$Sample]
    
   ### Go back to PCA plot and check what happned 
   ### perform combat normalization again after removal of sample
   rawcount <- as.matrix(rawcount)
   adjusted_counts <- ComBat_seq(rawcount, batch=anno$Batch, group=anno$Condition) ##In ComBat-seq, user may specify biological covariates, whose signals will be preserved in the adjusted data. I
   
  

############################### Create DESeq2 datasets #############################
library(DESeq2)


##dds <- DESeqDataSetFromMatrix(countData = rawcount, colData = anno, design = ~Condition )   ##rawcount ## simpledesign
## dds <- DESeqDataSetFromMatrix(countData = rawcount, colData = anno, design =  ~Batch+Condition )  ###USE this one if you have extra col in anno data with Batch info
dds = DESeq2::DESeqDataSetFromMatrix(countData = adjusted_counts, colData = anno, design = ~ Condition)  ##https://github.com/zhangyuqing/ComBat-seq/issues/7

##When considering batch effects in group design of Deseq2, it takes into account the mean differences across batch, 
##not necessarily the variance differences. ComBat-Seq is designed to address both mean and variance batch effects.
###In theory, no, you do not need to include batch as a covariate any more. However, you can always try both and evaluate the results. 
## https://github.com/zhangyuqing/ComBat-seq/issues/7



#View(counts(dds))

dds <- estimateSizeFactors(dds)

vst <- vst(dds, blind=TRUE)  ### Transform counts for data visualization #options (1) vst (2) rld

##normalized_counts <- counts(dds, normalized=TRUE)  ## extract normalization count after executing Deseq2 for visualization purpose
normalized_counts <- as.data.frame(assay(vst))

plotPCA(vst, intgroup="Condition")  ### Plot PCA 

## Run DESEQ2
dds <- DESeq(dds)

##ensure your data is a good fit for the DESeq2 model
plotDispEsts(dds)

################# contrast based  comparison ##########################

# Define conditions (for contrast) that you want to compare if you have more than one #control #case

firstC<-"Condition_A"       #case1 #case2 #case3 etc          
SecondC <-"Condition_B"     
p.threshold <- 0.05   ##define threshold for filtering

#In case of multiple comparisons ## we need to change the contrast for every comparision
contrast<- c("Condition",firstC,SecondC)

res <- results(dds, contrast=contrast)  ## extract result dataframe 
View(as.data.frame(res))

### Valcono plot
library(EnhancedVolcano)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                #pCutoff = 1e-05,
                #FCcutoff = 1,
                y = 'pvalue')   ## Default cut-off for log2FC is >|2| and for P value is 10e-6. USE  pCutoff = 10e-6, FCcutoff = 2.0 

#?EnhancedVolcano


res$threshold <- as.logical(res$padj < p.threshold)  #Threshold defined earlier#creating col with logic 

nam <- paste('down_in',firstC, sep = '_')   
#res$nam <- as.logical(res$log2FoldChange < 0)
res[, nam] <- as.logical(res$log2FoldChange < 0)  #adding extra information, in which condition it's down or up 

genes.deseq <- row.names(res)[which(res$threshold)]   ### list of gene with Padjust < defined threshold
genes_deseq2_sig <- as.data.frame(res[which(res$threshold),])


########### Plots normalized count of top 20 genes ## sorted based on padjust and filter by |logFC| >=1

res$gene <- row.names(res)
View(as.data.frame(res))

# Order results by padj values

#library(dplyr)
library(tidyverse)

top20 <- res %>% 
  as.data.frame %>%
  arrange(padj) %>% 	#Arrange rows by padj values
  filter(abs(log2FoldChange) >=1.5) %>%   #filter based on logFC
  pull(gene) %>% 		#Extract character vector of ordered genes
  head(n=20) 		#Extract the first 20 genes

top20_norm <- as.data.frame(normalized_counts[rownames(normalized_counts) %in% top20,])

top20_norm_v2 <- top20_norm ## will use later for heatmap

top20_norm <- (top20_norm+1) ## in later step to remove infinity bias due to log

top20_norm$gene <-  row.names(top20_norm)  
top20_norm <- top20_norm %>% 
  pivot_longer(!gene, names_to = "samplename", values_to = "normalized_counts") # Gathering the columns to have normalized counts to a single column

# Create tibbles including row names #for anno
mov10_meta <- anno %>% 
  rownames_to_column(var="samplename") %>% 
  as_tibble()

top20_norm <- inner_join(mov10_meta, top20_norm)  ## we are merging anno of 20 with previously extract format to draw dot plot

################
## plot using ggplot2

ggplot(top20_norm) +
  geom_point(aes(x = gene, y = normalized_counts, color = Condition)) +
  ## scale_y_log10() +  ##want to scale it or not??
  xlab("Genes") +
  ylab("log 10 CPM Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes with abs(logFC) =>1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))


##################
heat_colors <- colorRampPalette(c("blue",'white','red'))(n=40)

### Run pheatmap
library(pheatmap) 
pheatmap(top20_norm_v2 , 
         color = heat_colors, 
         cluster_rows = T,
         cluster_cols = F,
         show_rownames = T,
         annotation_col = anno[,c(1,3)], 
         border_color = NA, 
         fontsize = 10, 
         scale = "row",   ## for VST_count you may not need Scaling
         fontsize_row = 10, 
         height = 20)



file <- paste('Deseq2_',firstC,'_v_',SecondC,'_results_significant_padj',p.threshold,'.csv',sep = '') 
all_results <- paste('Deseq2_',firstC,'_v_',SecondC,'_all_results.csv',sep = '')

res <- as.data.frame(res)
View(res)

write.table(res,all_results,sep = ",")  ## no LogFC threshold




#####################################################################
############## functional Enrichment analysis #######################
####################################################################


### Option (1) continue from DGE analysis or (2) upload data from saved file of DEG 
#genes_deseq2_sig <- read.csv("Deseq2_case1_v_Control_results_significant_padj0.05.csv") 

###################### Extract various types of gene ids from Biomart ########################
library("biomaRt")

#new_config <- httr::config(ssl_verifypeer = FALSE) ############For certificate error
#httr::set_config(new_config, override = FALSE)     ############For certificate error

### define the mart for h_sapiens

#ensembl_mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")  ## either this or following line
ensembl_mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "asia") ## takes little bit time

###### following lines is extarcting other alternate names of the hugo gene symobols 
###### in this case entrez gene ids
View(as.data.frame(res))
genes.entrezid <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"), filters = "hgnc_symbol", values = res$gene, mart = ensembl_mart)

#genes.entrezid = as.data.frame(genes.entrezid) ## if not defined as df

res <- as.data.frame(res)
merged <- merge(res, genes.entrezid, by.x= "gene", by.y="hgnc_symbol")


######### Rank all genes based on their fold change #########

#BiocManager::install("clusterProfiler", force = TRUE)
#BiocManager::install("pathview", force = TRUE)
#BiocManager::install("enrichplot", force = TRUE)

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# SET THE DESIRED ORGANISM HERE ###https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
organism = "org.Hs.eg.db"  ## search other organism annotations here http://bioconductor.org/packages/release/BiocViews.html#___OrgDb

#BiocManager::install(organism, character.only = TRUE, force = TRUE)

library(organism, character.only = TRUE)


#We will take the log2FoldChange value from previously saved significant results file
#Deseq2_case1_v_Control_results_significant.csv


# we want the log2 fold change 
original_gene_list <- merged$log2FoldChange
print(original_gene_list)

# name the vector
names(original_gene_list) <- merged$entrezgene_id
print(original_gene_list)

# omit any NA values  ## excluding genes where no entrez ids available ## you may loose some information 
gene_list1<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list1 = sort(gene_list1, decreasing = TRUE)

gene_list1
#### Gene Set Enrichment  of Gene Ontology #####

library(stats)

keytypes(org.Hs.eg.db)

gse <- gseGO(geneList=gene_list1, 
             ont ='ALL', #"ALL", #### Try GO with all different ont methods parameter  ## BP = Biological Processes, CC= Cellular component, MF = Molecular functions
             keyType = "ENTREZID", 
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")
#?gseGO

# require(DOSE)
 view(as.data.frame(gse))
 #dotplot(gse, showCategory=10, split=".sign", orderBy = "X")
 gseaplot(gse, geneSetID="GO:0030198")
 
### Lets exlore other functions with a sample dataset and see what analysis we can do with
## a list of differentially expressed genes
###### geneList dataset of DOSE package #######
 data(geneList)
 print(geneList)
 
 
 gsecc <- gseGO(geneList=geneList, ont="CC", OrgDb=org.Hs.eg.db, verbose=F)
 view(as.data.frame(gsecc)) ## or use ## head(summary(gsecc))
 gseaplot(gsecc, geneSetID="GO:0000775")

##GO Enrichment Analysis of a gene set. 
##Given a vector of genes, enrichGO function will return the 
##enrichment GO categories after FDR control.



#library(clusterProfiler)
#library(org.Hs.eg.db)
library(enrichplot)
library(GOSemSim)
library(ggnewscale)
library(DOSE)

View(as.data.frame(geneList))
gene <- names(geneList)[abs(geneList) > 2] 
gene
## you may to consider additional filter based on P-adjusted Values
## Last time we created df of significant DE gene i.e. 'genes_deseq2_sig' 
## you can apply the above filter of logFC 
##something like '
#   subdf <- genes_deseq2_sig[abs(genes_deseq2_sig$log2FoldChange) > 1.5,] 
#  gene_list <- subdf$log2FoldChange
#   names(gene_list) <- subdf$hgnc_symbol    ### you may want to consider entrez gene conversion


#?enrichGO
ego <- enrichGO(gene  = gene,  
                universe      = names(geneList), ##background genes
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH", # "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

View(as.data.frame(ego))
##### Visualization of enrichGO ######

d <- godata('org.Hs.eg.db', ont="BP") #prepare GO DATA for measuring semantic similarity ## required for next step
ego2 <- pairwise_termsim(ego, method="Wang", semData = d) #enrichment result #method of calculating the similarity between nodes #GOSemSimDATA object
emapplot(ego2)
emapplot_cluster(ego2)
view(as.data.frame(ego2))
#cnetplot(ego2, categorySize="pvalue", foldChange=gene_list)


###In the following example, we selected fold change above 1 as the differential genes 
##and analyzing their disease association.

#### enrich DO (Disease Ontology) #####
##http://bioconductor.org/packages/release/bioc/html/DOSE.html

library(ggupset)

data(geneList) ## loaded the same data 
gene = names(geneList)[abs(geneList) > 1.5] ## applied theshold 
head(gene)
X = enrichDO(gene,ont = "DO", 
              pvalueCutoff=0.05, 
              pAdjustMethod = "BH", 
              universe = names(geneList),
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)

#The readable is a logical parameter, 
#indicates whether the entrezgene IDs will mapping to gene symbols or not
head(X)
View(as.data.frame(X))

#setReadable function helps to convert entrezgene IDs to gene symbols
X <- setReadable(X, 'org.Hs.eg.db')
View(as.data.frame(X))

## Visualization of enrichDO results ##
barplot(X, showCategory=15)
dotplot(X)

#gene may belong to multiple annotation categories, 
#cnetplot function to extract the complex association between genes and diseases
cnetplot(X, categorySize="pvalue", foldChange=geneList)

#upsetplot is an alternative to 
#cnetplot for visualizing the complex association between genes and diseases.

upsetplot(X)


###### KEGG Enrichment Analysis #######

library(clusterProfiler)

## KEGG pathway over-representation analysis (ORA)
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

view(as.data.frame(kk))
#?enrichKEGG ## you can defined background using 'universe'  ## universe = names(geneList),



## KEGG module ORA (over-representation analysis)
#KEGG Module is a collection of manually defined function units. In some situation,
#KEGG Modules have a more straightforward interpretation

mkk <- enrichMKEGG(gene = gene,
                   organism = 'hsa',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
View(as.data.frame(mkk))                   

## KEGG module GSEA (gene set enrichment analysis) ##

mkk2 <- gseMKEGG(geneList = geneList,
                 organism = 'hsa',
                 pvalueCutoff = 1)
View(as.data.frame(mkk2)) 


####### Visualize enriched KEGG pathways #########

## To view the KEGG pathway,use the browseKEGG function,
#which will open a web browser and highlight enriched genes.

browseKEGG(kk, 'hsa04110')

###use the pathview() function from the pathview to visualize enriched KEGG 
##pathways identified by the clusterProfiler package

library("pathview")
hsa04110 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa04114",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))



 


















