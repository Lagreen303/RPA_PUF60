# For data management
install.packages('tidyverse')
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
# For visualisation
install.packages('pheatmap')
install.packages("DOSE")
install.packages("enrichplot")
install.packages("ggupset")

library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
  library(enrichplot)
library(edgeR)
#Load results object from DGE script!

# Load your normalized count matrix from a text file
puf60count <- read.table("count_files/genename_counts.txt", header = TRUE, sep = "\t", row.names = 1)
puf60_group <- ifelse(colnames(puf60count) %in% c("P002", "P001", "SOT103", "SOT149", "SOT211", "SOT325"), 
                      "PUF60 Variant", "No Variant")
puf60.metaData <- data.frame(Sample = colnames(puf60count), Group = puf60_group)
group <- factor(puf60.metaData$Group)
dge_list <- DGEList(counts = puf60count, group = group) 

# Normalize and filter (check you have not already done this)
#dge_list <- calcNormFactors(dge_list)
#keep <- filterByExpr(dge_list)
#dge_list <- dge_list[keep, ]

dge_list <- estimateDisp(dge_list)
fit <- glmFit(dge_list)
lrt <- glmLRT(fit)
results <- topTags(lrt, n = Inf)

write.csv(results, "DGE_results.csv")

#################

df=read.csv("DGE_results.csv", header=T)
gene_list<-df$logFC
names(gene_list)<-df$X
gene_list<-sort(gene_list, decreasing = T)

gse <- gseGO(geneList = gene_list,
             ont = "ALL",
             keyType = "SYMBOL",
             minGSSize = 10,
             maxGSSize = 500,
             pvalueCutoff = 0.05,
             verbose=TRUE,
             OrgDb = org.Hs.eg.db)
Alldot<-dotplot(gse, showCategory=10, split=".sign") + ggtitle("All PUF60 Variants")

cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)

pairwisegse<-pairwise_termsim(gse)
emapplot(pairwisegse, showCategory = 100, repel=T)
ridgeplot(gse, label_format = 70) + labs(x = "Enrichment Distribution")

gseaplot(gse, by = "all", title = gse$Description[2], geneSetID = 2)

##KEGG GSEA

# Load your normalized count matrix from a text file
keggcount <- read.table("count_files/genename_counts.txt", header = TRUE, sep = "\t", row.names = 1)
kegg_group <- ifelse(colnames(keggcount) %in% c("P002", "P001", "SOT103", "SOT149", "SOT211", "SOT325"), 
                      "PUF60 Variant", "No Variant")
kegg.metaData <- data.frame(Sample = colnames(keggcount), Group = kegg_group)
group <- factor(kegg.metaData$Group)
dge_list <- DGEList(counts = keggcount, group = group) 

# Normalize and filter (check you have not already done this)
#dge_list <- calcNormFactors(dge_list)
#keep <- filterByExpr(dge_list)
#dge_list <- dge_list[keep, ]

dge_list <- estimateDisp(dge_list)
fit <- glmFit(dge_list)
lrt <- glmLRT(fit)
results <- topTags(lrt, n = Inf)

write.csv(results, "DGE_ENSEMBL_results.csv")

df=read.csv("DGE_ENSEMBL_results.csv", header=T)
df$X <- sub("\\..*", "", df$X)

gene_list<-df$logFC
names(gene_list)<-df$X
gene_list<-sort(gene_list, decreasing = T)

ids<-bitr(names(gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
df2 = df[df$X %in% dedup_ids$ENSEMBL,]
df2$Y = dedup_ids$ENTREZID

kegg_gene_list<- df2$logFC
names(kegg_gene_list)<-df2$Y
kegg_gene_list<-na.omit(kegg_gene_list)
kegg_gene_list=sort(kegg_gene_list, decreasing = T)

kegg_organism<- "hsa"
kk2<- gseKEGG(geneList = kegg_gene_list,
              organism = kegg_organism,
              nPerm = 10000,
              minGSSize = 3,
              maxGSSize = 800,
              pvalueCutoff = 0.05,
              pAdjustMethod = "none",
              keyType = "ncbi-geneid")

#KEGG plots
allkegg<-dotplot(kk2, showCategory = 10, title = "All PUF60 Variants" , split=".sign")

pairwisekk2<-pairwise_termsim(kk2)
emapplot(pairwisekk2,showCategory = 11, repel=T)

cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)

ridgeplot(kk2) + labs(x = "enrichment distribution")

gseaplot(kk2, by = "all", title = kk2$Description[2], geneSetID = 2)

#pathwaay pic
library(pathview)

dme2<-pathview(gene.data = kegg_gene_list, pathway.id = "hsa05340", species = kegg_organism,)

######################
# Load necessary libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(enrichplot)

# Read the differential gene expression results
df <- read.csv("DGE_results.csv", header = TRUE)
gene_list <- df$logFC
names(gene_list) <- df$X
gene_list <- sort(gene_list, decreasing = TRUE)

# Convert gene symbols to Entrez IDs
gene_list_entrez <- bitr(names(gene_list), fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)
gene_list <- gene_list[gene_list_entrez$SYMBOL]
names(gene_list) <- gene_list_entrez$ENTREZID

# Perform GSEA for KEGG pathways
gse_kegg <- gseKEGG(geneList = gene_list,
                    organism = "hsa",
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 0.05,
                    verbose = TRUE)

# Visualize the results with a dot plot
dotplot(gse_kegg, showCategory = 10, split = ".sign")
