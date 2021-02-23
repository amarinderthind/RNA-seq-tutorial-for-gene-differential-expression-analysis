### Author##############
## Amarinder Singh Thind


# Install and load packages

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("DESeq2")
#BiocManager::install("edgeR")
#BiocManager::install("biomaRt")
#BiocManager::install('PCAtools')
#BiocManager::install('EnhancedVolcano')

library(edgeR)
library(DESeq2)
library("biomaRt")

###################### load the raw count matrix #######################

setwd("./") #Path_to_working_directory

rawcount<-read.table ("RawGeneCounts.tsv",header=TRUE,  sep="\t",  row.names=1)

######################  Filter for coding genes (In case want to filter non-coding Genes) ########################

mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
all_coding_genes <- getBM(attributes = c( "hgnc_symbol"), filters = c("biotype"), values = list(biotype="protein_coding"), mart = mart)
rawcount <- rawcount[row.names(rawcount) %in%  all_coding_genes$hgnc_symbol,]

###################### Data annotation  #################################

anno <-read.table ("Annotation_of_samples.csv",header=TRUE,  sep=",") ##In this case Two coulmns (a) sample (b) Condition
rownames(anno) <- anno$sample

# Define conditions (for contrast) that you want to compare if you have more than one #control #case
# This is pair-wise comparison, so only consider one pair at one time

firstC<-"case1"       #case1 #case2 #case3 etc          
SecondC <-"Control"     
p.threshold <- 0.05   ##define threshold for filtering


### subset raw and conditional data for defined pairs

anno <- anno[(anno$Condition ==firstC |anno$Condition ==SecondC),]

anno <- anno[anno$sample %in% names(rawcount),]
rawcount <- rawcount[,names(rawcount) %in% anno$sample]

############################### Create DESeq2 datasets #############################

dds <- DESeqDataSetFromMatrix(countData = rawcount, colData = anno, design = ~Condition )

## Run DESEQ2
dds <- DESeq(dds)

################# contrast based  comparison ##########################

#In case of multiple comparisons ## we need to change the contrast for every comparision
contrast<- c("Condition",firstC,SecondC)

res <- results(dds, contrast=contrast)

### Valcono plot
library(EnhancedVolcano)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')   ## Default cut-off for log2FC is >|2| and for P value is 10e-6. USE  pCutoff = 10e-6, FCcutoff = 2.0 


res$threshold <- as.logical(res$padj < p.threshold)  #Threshold defined earlier

nam <- paste('down_in',firstC, sep = '_')
#res$nam <- as.logical(res$log2FoldChange < 0)
res[, nam] <- as.logical(res$log2FoldChange < 0)

genes.deseq <- row.names(res)[which(res$threshold)]
genes_deseq2_sig <- res[which(res$threshold),]



file <- paste('Deseq2_',firstC,'_v_',SecondC,'_results_significant_padj',p.threshold,'.csv',sep = '')
all_results <- paste('Deseq2_',firstC,'_v_',SecondC,'_all_results.csv',sep = '')

write.table(genes_deseq2_sig,file,sep = ",")
write.table(res,all_results,sep = ",")

################### PCA and Heat-MAp Plots ############################

## Varinace transformation vst or rlog
vsd <- vst(dds, blind=FALSE)   #Variance type (a) Vst or (b) rlog
#rld <- rlog(dds, blind=FALSE) 

library(ggplot2)
###### PCA with design consideration ###
pcaData <- plotPCA(vsd, intgroup=c("Condition", "sample"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=sample, shape=Condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

##heatmap
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
library('pheatmap')
sampleDistMatrix <- as.matrix(sampleDists)

rownames(sampleDistMatrix) <- paste(vsd$Condition, vsd$sample, sep="-")

colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)



############## edgeR  ##########################


dge <- DGEList(counts=rawcount, group=anno$Condition)

# Normalize by total count
dge <- calcNormFactors(dge, method = "TMM")

# filter out lowly expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep,keep.lib.sizes=FALSE]  # It is recommended to recalculate the library sizes of the DGEList object after the filtering,
                                        #although the downstream analysis is robust to whether this is done or not.

# You can also filter the expression matrix based on the treatment factors of scientific interest 
#keep <- filterByExpr(y, group=Condition)


## PCA ## for more details, please visit following link
##https://bioconductor.org/packages/release/bioc/vignettes/PCAtools/inst/doc/PCAtools.html
library(PCAtools)

cpmlog <- cpm(dge, log = TRUE, prior.count = 1)

p <-pca(cpmlog, metadata = anno, removeVar = 0.1) ## -- removing the lower 10% of variables based on variance
#biplot(p)
plotloadings(p)
 
biplot(p,
       lab = paste0(p$metadata$sample),
       colby = 'Condition',
       hline = 0, vline = 0,
       legendPosition = 'right')


# Create the contrast matrix
design.mat <- model.matrix(~ 0 + dge$samples$group)
colnames(design.mat) <- levels(dge$samples$group)
design.mat

# Estimate dispersion parameter for GLM
dge <- estimateGLMCommonDisp(dge, design.mat)
dge <- estimateGLMTrendedDisp(dge, design.mat) 
dge<- estimateGLMTagwiseDisp(dge,design.mat)
# Plot mean-variance
#plotBCV(dge)


# Model fitting 
##  EdgeR glmLRT vs glmQLFTest ## https://support.bioconductor.org/p/84291/

##  both of the methods will work for your data set, the QL F-test is probably the better choice. 
##There are some situations where the QL F-test doesn't work well - for example, if you don't have replicates,
##you'd have to supply a fixed dispersion, which defeats the whole point of modelling estimation uncertainty.
##Another situation is where the dispersions are very large and the counts are very small, whereby some of the approximations in the QL framework seem to fail.

fit.edgeR <- glmQLFit(dge, design.mat)  #glmFit

# Differential expression

contrasts.edgeR <- makeContrasts(case1 - Control, levels=design.mat)    ##FirstC-SecondC ##Define 

qlf.edgeR <-glmQLFTest(fit.edgeR, contrast=contrasts.edgeR)  # glmLRT

##### DGE at padjust 0.05

# Access results tables
edgeR_results <- qlf.edgeR$table
sig.edgeR <- decideTestsDGE(qlf.edgeR, adjust.method="BH", p.value = p.threshold)
#View(sig.edgeR) 
significant_table <- edgeR_results[which(sig.edgeR != 0),]
significant_table$gene <- row.names(significant_table)
genes.edgeR <- row.names(edgeR_results)[which(sig.edgeR != 0)]

edgeR_results$genes <- row.names(edgeR_results)


file_sigTab <- paste('edgeR_',firstC,'_v_',SecondC,'_results_significant_padj',p.threshold,'.csv',sep = '')
file_allRes <- paste('edgeR_',firstC,'_v_',SecondC,'_all_results.csv',sep = '')

write.table(significant_table,file_sigTab,sep = ",")
write.table(edgeR_results,file_allRes,sep = ",")


################# Overlapped genes between deseq2 and edgeR  ##########

library(gplots)

venn(list(edgeR = genes.edgeR, DESeq2 = genes.deseq))
overlapped_genes <- intersect(genes.deseq,genes.edgeR)


file_common <- paste('Common_DEG_deseq2_edgeR_',firstC,'_v_',SecondC,'.csv',sep = '')
write.table(overlapped_genes,file_common,sep = ",", row.names = F)

############ Quick enrichment analysis ##################
#BiocManager::install("ReactomePA")

library(ReactomePA)

all <- overlapped_genes   ## retreive EntrezGene id's

genes=getBM(attributes = c("hgnc_symbol", "entrezgene_id"), filters = "hgnc_symbol", values = all, bmHeader = T, mart = mart)

genes1 <- genes$`NCBI gene (formerly Entrezgene) ID` 

#?enrichPathway #pvalueCutoff=0.02, #pAdjustMethod = "BH", qvalueCutoff = 0.01,
x <- enrichPathway(gene=genes1,  pvalueCutoff=0.05,readable=T)

#head(as.data.frame(x))
barplot(x, showCategory=10)
dotplot(x, showCategory=10)
emapplot(x)
cnetplot(x, categorySize="pvalue", foldChange=genes1)
emapplot(x, color="pvalue")
viewPathway("Extracellular matrix organization", readable=TRUE, foldChange=genes1)   ## it's an example


