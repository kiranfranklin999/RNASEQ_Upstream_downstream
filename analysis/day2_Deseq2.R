# script to perform differential gene expression analysis using DESeq2 package

##set working directory
setwd("D:/Bulk-RNA_training")

## install required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("biomaRt")
BiocManager::install("EnhancedVolcano")
BiocManager::install("airway")

install.packages("pheatmap")
install.packages("tidyverse")


# load libraries
library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
library(biomaRt)
library(pheatmap)

###Step 1: Import Data
# read in counts data
counts_data <- read.csv('counts_data.csv')
head(counts_data)

# read in sample info
colData <- read.csv('sample_info.csv')

# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))

# are they in the same order?
all(colnames(counts_data) == rownames(colData))

# Step 2: construct a DESeqDataSet object ----------

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                       colData = colData,
                       design = ~ cellLine + dexamethasone)

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# set the factor level
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

###Step 3: Quality control
##vst normalization for PCA
vsd = vst(dds) 
pca_plot=DESeq2::plotPCA(vsd, intgroup='dexamethasone')
pca_plot+ geom_text(size=2.4,aes(label=colData$cellLine),vjust=2)

## for labels as same names
#pca_plot+ geom_text(size=2.4,aes(label=name),vjust=2)
#pca_plot

###Step 4: Run DESeq 
dds <- DESeq(dds)
res <- results(dds, contrast = c('dexamethasone', 'untreated', 'treated'))

# Explore Results ----------------

summary(res)


###Step 5: Filtering DEGs
sig_genes = as.data.frame(res[which(res$padj < 0.05), ])

upregulated = sig_genes[which(sig_genes$log2FoldChange > 0.6), ]
downregulated = sig_genes[which(sig_genes$log2FoldChange < -0.6), ]



#####Step 6: Gene ID conversion
sig_genes$EnsemblID<- rownames(sig_genes)
gene_list= sig_genes$EnsemblID

ensembl_human=useDataset("hsapiens_gene_ensembl",mart = useMart("ensembl"))

# ensembl_human <- useMart(biomart = "ensembl", 
#                          dataset = "hsapiens_gene_ensembl",
#                          host = "useast.ensembl.org")

map<- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                    filters = "ensembl_gene_id",
                    values=gene_list,
                    mart = ensembl_human)
colnames(map)[1]= "EnsemblID"
deg_1= as_tibble(sig_genes) %>% left_join(map, by=("EnsemblID"), multiple="all")
head(deg_1)


###Step 7 : Visualizations
##Volcano plot
vplot= EnhancedVolcano(deg_1,
                       lab = deg_1$hgnc_symbol,
                       x = 'log2FoldChange',
                       y = 'padj',
                       pCutoff = 0.05,
                       FCcutoff = 1,                
                       title = '',
                       subtitle = "Treated Vs Untreated",                     ##
                       #caption = "|Log2FC| > 0.58; padj < 0.05",
                       pointSize = 3.0,
                       colAlpha = 0.7,
                       legendLabels=c('Not sig.','Log2FC','padj',
                                      'padj & Log2FC'),
                       #drawConnectors = TRUE,
                       labSize = 4.0)
vplot


##heatmap of top 50 genes

##selecting top 50 degs
deg_sort = deg_1[order(abs(deg_1$log2FoldChange), decreasing = T),]
top_50= deg_sort[1:50,]

####converting counts to normalized counts and applying scaling
mat= counts(dds, normalized=T)[top_50$EnsemblID,]
mat.z=t(apply(mat,1,scale))
colnames(mat.z)= colnames(mat)
#rownames(mat.z)= top_50$hgnc_symbol


annot_col= as.data.frame(colData[,2])
rownames(annot_col)= rownames(colData)
colnames(annot_col)= "dexamethasone"
pheatmap(mat.z, cluster_cols = T, fontsize_row = 5, annotation_col = annot_col) -> hplot

###Step 8: Exporting results

write.csv(deg_1, file="DEG_Result.csv", row.names = F)
png("Heatmap.png", width = 600, height = 800)
hplot
dev.off()
png("Volcano_plot.png", width = 800, height = 800)
vplot
dev.off()