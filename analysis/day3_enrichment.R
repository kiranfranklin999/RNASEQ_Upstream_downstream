###set working directory
setwd("D:/Bulk-RNA_training")

###installation of required packages
BiocManager::install("reactome.db")
BiocManager::install("ReactomePA")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DOSE")
BiocManager::install("enrichplot")

### loading libraries
library(clusterProfiler)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(tidyverse)
library(DOSE)
library(enrichplot)

####input result object (res) of DESeq2 or use existing result object from your environment
#saveRDS(res, file="res.RDS")
res= readRDS("res.RDS")
head(res)

###converting ensemblID to entrez_ID/ncbi IDs
clusterProfiler::bitr(rownames(res), 'ENSEMBL', 'ENTREZID', OrgDb = org.Hs.eg.db) -> ncbi_list
#res$EnsemblID[-which(ncbi_list$ENSEMBL %in% res$EnsemblID)] -> unmapped

colnames(ncbi_list)[1]<- "EnsemblID"
res$EnsemblID<- rownames(res)
as_tibble(res) %>% left_join(ncbi_list, by="EnsemblID") -> res_mapped
head(res_mapped)

###Gene ontology 
##extracting significant genes
sig_genes = as.data.frame(res_mapped[which(res_mapped$padj < 0.05),])
res2<- sig_genes[which(abs(sig_genes$log2FoldChange)>=1),]
list2<- res2$EnsemblID

enGO=enrichGO(list2, keyType = "ENSEMBL", OrgDb = org.Hs.eg.db, ont = "BP", readable = T)  ##can change parameter "BP" to "ALL", "MF" or "CC"

enGO2<- enGO[which(enGO@result$p.adjust<= 0.05),]


###Gene ontology (Gene Set Enrichment analysis)

##creating named gene list for all genes
ngenes<- res_mapped$stat
names(ngenes)<- res_mapped$ENTREZID
ngenes<- sort(ngenes, decreasing = T)
ngenes[-which(is.na(names(ngenes))==T)] -> ngenes
ngenes[unique(names(ngenes))]-> ngenes

###Running GSEA for GO
gse_res= gseGO(geneList = ngenes, OrgDb = org.Hs.eg.db,keyType = "ENTREZID", ont = "BP")

gsea_mapped <- setReadable(gse_res, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

######################################################################################################################################


###Pathway enrichment 

##pathway enrichment using Reactome

###overrepresentation
enp<- enrichPathway(gene=res2$ENTREZID, pvalueCutoff = 0.05, readable=TRUE, organism = "human")
View(enp@result)

#enp2<- enp[which(enp@result$p.adjust<= 0.05),]

##GSEA
enp_gsea<- gsePathway(ngenes,organism = "human") #uses Reactome

enp_gsea_kegg<- gseKEGG(ngenes, organism = "hsa") ##using KEGG
enp_gsea <- setReadable(enp_gsea, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

######################################################################################################################################


###vizualisation

###network plot of genes as associated GO terms

##relevant terms from ORA biological process
as_tibble(enGO@result) %>% filter(str_detect(Description,"cortico")) -> relevant_go_ora

as_tibble(gse_res@result) %>% filter(str_detect(Description,"cortico")) -> relevant_go_gsea

x <- enGO[enGO$Description %in% relevant_go_ora$Description, asis = TRUE]

netplot=cnetplot(x, circular = TRUE, colorEdge = TRUE, cex_label_gene=0.5)
netplot



##dotplot of enriched GO terms
dotplot(enGO, showCategory= 10)


y <- as.data.frame(enGO)
head(y)
ggplot(y[c(250,360,449,565, 590, 600),], # you can replace the numbers to the row number of pathway of your interest
       aes(x = Count, y = Description)) + 
  geom_point(aes(size = Count, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("GO biological process enrichment")



######################################################################################################################################

####Optional

###Disease enrichment
##over representation
enDO<- enrichDO(res2$ENTREZID,pvalueCutoff = 0.05, readable = T, ont = "DO")

###gsea
gsea_DO<- gseDO(ngenes)

View(enDO@result)