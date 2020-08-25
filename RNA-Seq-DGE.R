### Author##############
## Amarinder Singh Thind


# Install and load packages

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("edgeR")
BiocManager::install("biomaRt")
BiocManager::install('PCAtools')

library(edgeR)
library(DESeq2)
library("biomaRt")

###################### load the raw count matrix #######################

setwd("Path_to_working_directory")

rawcount<-read.table ("RawGeneCounts.tsv",header=TRUE,  sep="\t",  row.names=1)

######################  Filter for coding genes ########################

mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
all_coding_genes <- getBM(attributes = c( "hgnc_symbol"), filters = c("biotype"), values = list(biotype="protein_coding"), mart = mart)
rawcount <- rawcount[row.names(rawcount) %in%  all_coding_genes$hgnc_symbol,]

######################  Filter low count gene  #########################

keep <- rowSums(cpm(rawcount)>1) >= 5   ## depends case to case and on the number of samples
rawcount<- rawcount[keep,]

###################### Data annotation  #################################

anno <-read.table ("Annotation_of_samples.csv",header=TRUE,  sep=",") ##In this case Two coulmns (a) sample (b) Condition
rownames(anno) <- anno$sample

# Define conditions that you want to compare if you have more than one #control #case
# This is pair-wise comparison, so only consider one pair at one time

firstC<-"Case1"       #case1 #case2 #case3 etc          
SecondC <-"Control"     
p.threshold <- 0.05   ##define threshold for filtering


### subset raw and conditional data for defined pairs

anno <- anno[(anno$Condition ==firstC |anno$Condition ==SecondC),]
rawcount <- rawcount[,names(rawcount) %in% anno$sample]

############################### Create DESeq2 datasets #############################

dds <- DESeqDataSetFromMatrix(countData = rawcount, colData = anno, design = ~Condition )

## Run DESEQ2
dds <- DESeq(dds)

################# contrast based  comparison ##########################

#In case of multiple comparisons ## we need to change the contrast for every comparision
contrast<- c("Condition",firstC,SecondC)

res <- results(dds, contrast=contrast)
res$threshold <- as.logical(res$padj < p.threshold)  #Threshold defined earlier

nam <- paste('down_in',firstC, sep = '_')
#res$nam <- as.logical(res$log2FoldChange < 0)
res[, nam] <- as.logical(res$log2FoldChange < 0)

genes.deseq <- row.names(res)[which(res$threshold)]

genes_deseq2_sig <- res[which(res$threshold),]

file <- paste('Deseq2_',firstC,'_v_',SecondC,'_results_significant_padj0.05.csv',sep = '')
all_results <- paste('Deseq2_',firstC,'_v_',SecondC,'_all_results.csv',sep = '')

write.table(genes_deseq2_sig,file,sep = ",")
write.table(res,all_results,sep = ",")

################### PCA and Heat-MAp Plots ############################

library(PCAtools)

cpmcount <- cpm(rawcount)
p <-pca(cpmcount, metadata = anno, removeVar = 0.1)
#biplot(p)
plotloadings(p)
 
biplot(p,
       lab = paste0(p$metadata$sample),
       colby = 'Condition',
       hline = 0, vline = 0,
       legendPosition = 'right')


## Varinace transformation vst or rlog
vsd <- vst(dds, blind=FALSE)   #Variance type (a) Vst or (b) rlog
#rld <- rlog(dds, blind=FALSE) 

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


############ Quick enrichment analysis ##################
BiocManager::install("ReactomePA")
library(ReactomePA)

all <- genes.deseq   ## retreive EntrezGene id's

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


############## edgeR  ##########################


dge <- DGEList(counts=rawcount, group=anno$Condition)

# Normalize by total count
dge <- calcNormFactors(dge, method = "TMM")

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

fit.edgeR <- glmFit(dge, design.mat)

# Differential expression

contrasts.edgeR <- makeContrasts(Case1 - Control, levels=design.mat)    ##FirstC-SecondC ##Define 


lrt.edgeR <- glmLRT(fit.edgeR, contrast=contrasts.edgeR)

##### DGE at padjust 0.05

# Access results tables
edgeR_results <- lrt.edgeR$table
sig.edgeR <- decideTestsDGE(lrt.edgeR, adjust.method="BH", p.value = p.threshold)
#View(sig.edgeR) 
significant_table <- edgeR_results[which(sig.edgeR != 0),]
significant_table$gene <- row.names(significant_table)
genes.edgeR <- row.names(edgeR_results)[which(sig.edgeR != 0)]

edgeR_results$genes <- row.names(edgeR_results)


file_sigTab <- paste('edgeR_',firstC,'_v_',SecondC,'_results_significant_padj0.05.csv',sep = '')
file_allRes <- paste('edgeR_',firstC,'_v_',SecondC,'_all_results.csv',sep = '')

write.table(significant_table,file_sigTab,sep = ",")
write.table(edgeR_results,file_allRes,sep = ",")


################# Overlapped genes between deseq2 and edgeR  ##########

library(gplots)

venn(list(edgeR = genes.edgeR, DESeq2 = genes.deseq))
overlapped_genes <- intersect(genes.deseq,genes.edgeR)


file_common <- paste('Common_DEG_deseq2_edgeR_',firstC,'_v_',SecondC,'.csv',sep = '')
write.table(overlapped_genes,file_common,sep = ",", row.names = F)

## Save session info
sessionInfo()
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
