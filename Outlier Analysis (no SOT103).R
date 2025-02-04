library(tidyverse)
library(edgeR)
library(factoextra)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(reshape2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
#Analysis WITHOUT P001 AND SOT103
#file names
file_names <- c("counts/P001.out.tab", "counts/P002.out.tab", 
                "counts/SOT111.out.tab", "counts/SOT112.out.tab", "counts/SOT149.out.tab", 
                "counts/SOT211.out.tab", "counts/SOT325.out.tab", "counts/SOT344.out.tab")
sot103count <- NULL  

# Loop through each file and extract count
for (file in file_names) {
  # Read the count file
  count_data <- read.table(file, header = FALSE, sep = "\t", row.names = NULL)
  
  # Extract count
  fourth_column <- count_data[, 4]  # Get the fourth column
  
  # Combine counts into the main matrix
  if (is.null(sot103count)) {
    sot103count <- fourth_column 
  } else {
    sot103count <- cbind(sot103count, fourth_column) 
  }
}

# Set row names to gene names
rownames(sot103count) <- count_data[, 1]

# Set column names 
colnames(sot103count) <- gsub(".out.tab", "", basename(file_names)) 

# Remove N_multimapping etc
unwanted_rows <- c("N_ambiguous", "N_multipmapping", "N_noFeature", "N_unmapped")
sot103count <- sot103count[!rownames(sot103count) %in% unwanted_rows, ]

sot103count <- sot103count[-1, ]

#save
write.table(sot103count, "Outliers/no_sot103count_Counts.txt", sep = "\t", quote = FALSE, col.names = NA)

##NORMALISATION AND FILTERING##
raw_data <- read.table("Outliers/no_sot103count_Counts.txt", header = T)
colSums(raw_data)

raw_filter <- rowSums(raw_data > 60) >= 3
filtered_data <- raw_data[raw_filter, ]
dge1 <- DGEList(counts = filtered_data, remove.zeros = TRUE)

dge1 <-calcNormFactors(dge1)
normalised_counts <- cpm(dge1)

write.table(normalised_counts, "Outliers/filtered_sot103count_Counts.txt", sep = "\t", quote = FALSE, col.names = NA)

##Gene name Swap##
# Load your count matrix from a text file
count_matrix <- read.table("Outliers/filtered_sot103count_Counts.txt", header = TRUE, sep = "\t", row.names = 1)

# Load the mapping file with Ensembl IDs and Gene Names
mapping_file <- read.table("gene_name_reference", header = FALSE, sep = "\t") 
colnames(mapping_file) <- c("EnsemblID", "GeneName")

count_df <- as.data.frame(count_matrix)

# Add a column for Ensembl IDs
count_df$EnsemblID <- rownames(count_df)

# Merge the count data with the mapping file
merged_df <- merge(count_df, mapping_file, by = "EnsemblID", all.x = TRUE)

# Create a new column for final gene names
merged_df <- merged_df %>%
  mutate(FinalGeneName = ifelse(duplicated(GeneName) | duplicated(GeneName, fromLast = TRUE), EnsemblID, GeneName))

# Remove the GeneName column and keep only FinalGeneName and counts
final_counts <- merged_df[, c("FinalGeneName", colnames(count_df)[-ncol(count_df)])]

# Set row names to FinalGeneName
rownames(final_counts) <- final_counts$FinalGeneName
final_counts$FinalGeneName <- NULL

# Write the final counts to a new file
write.table(final_counts, "Outliers/labelled_sot103count_counts.txt", sep = "\t", row.names = TRUE)

##########################################################################################################
##PCA##
sot103count <- read.table("Outliers/labelled_sot103count_counts.txt", header = TRUE, row.names=1)
sot103count.svd <- 
  sot103count %>% 
  t() %>% 
  prcomp(scale = T)

fviz_eig(sot103count.svd, addlabels = T)
fviz_pca_ind(sot103count.svd, repel = T)

fviz_pca_var(sot103count.svd, select.var = list(cos2 = 10), repel = T)


# Define the groups based on your sample names
sot103count_group <- ifelse(colnames(sot103count) %in% c("P002", "P001", "SOT103", "SOT149", "SOT211", "SOT325"), 
                         "PUF60 Variant", "No Variant")

# Create a new metadata data frame if needed
sot103count.metaData <- data.frame(Sample = colnames(sot103count), Group = sot103count_group)

# Generate the PCA biplot
fviz_pca_biplot(sot103count.svd, 
                select.var = list(name = c("TMEM1", "NOP1")),
                habillage = sot103count.metaData$Group, 
                label = "all", 
                repel = T,
                invisible = "quali",
                addEllipses = F) + 
  ggtitle("PUF60 Group Separation")

sot103PCA<-fviz_pca_biplot(sot103count.svd, 
                select.var = list(name = c("TMEM1", "NOP1")),
                habillage = sot103count.metaData$Group, 
                label = "all", 
                repel = TRUE,
                invisible = "quali",
                addEllipses = FALSE) + 
  ggtitle("A (SOT103 Excluded)") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1.5, "lines")
  ) +
  scale_color_manual(values = c("red", "blue"))

##DGE## 

# Create a DGEList object with a different name
dge_list <- DGEList(counts = sot103count, group = sot103count_group) 

# Estimate dispersion
dge_list <- estimateDisp(dge_list)

# Fit the model
fit <- glmFit(dge_list)

# Perform likelihood ratio test
lrt <- glmLRT(fit)

# Extract results
results <- topTags(lrt, n = Inf)
significant_genes <- results$table[abs(results$table$logFC) >2& results$table$PValue < 0.05 & results$table$FDR < 0.05, ]  
write.table(significant_genes, "Outliers/sot103count_significantgenes.txt", sep = "\t", row.names = TRUE)

#Sample Correlation heatmap (Change correlation power (44) as needed for visualisation)
melt_corr <-cor(sot103count) |>
  melt()
melt_corr$raised_cor <- (melt_corr$value)^6 

ggplot(data = melt_corr, aes(x=Var1, y=Var2, fill = raised_cor)) +
  geom_tile(color = "white") + scale_fill_gradient2(low = "blue",high="red",mid="white", midpoint = 0.5, limit =c(0,1), space ="Lab", name="Correllation")+
  theme_minimal()

###Volcano Plot###

#label whether genes are up or downregulated
results$table$diffexpressed<- "NO"
results$table$diffexpressed[results$table$logFC > 2 & results$table$FDR <0.05] <-"UP"
results$table$diffexpressed[results$table$logFC < -2 & results$table$FDR <0.05] <- "DOWN"
results$table$diffexpressed[results$table$logFC > -1 & results$table$logFC < 1 & results$table$FDR <0.05] <- "NO"

#label significant genes
results$table$genelabel <- ""
results$table$genelabel <- ifelse(rownames(results$table) == "VWCE"
                                  |rownames(results$table) == "ENSG00000276241"
                                  |rownames(results$table) == "OVCH1-AS1"
                                  |rownames(results$table) == "VSIG10"
                                  |rownames(results$table) == "RPL18AP3"
                                  |rownames(results$table) == "LINC02531"
                                  |rownames(results$table) == "TPRG1"
                                  |rownames(results$table) == "XK"
                                  |rownames(results$table) == "NEAT1"
                                  |rownames(results$table) == "LINC02506", rownames(results$table), "")


ggplot(as.data.frame(results$table)) +
  geom_point(aes(x = logFC, y = -log10(FDR), color = diffexpressed), size=2) +
  geom_text_repel(aes(label = results$table$genelabel, max.overlaps = 1, x = logFC, y = -log10(FDR))) +
  scale_colour_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-2, 2), linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed')+
  theme_classic() +
  theme(legend.position = "none")+
  labs(title = "A (SOT103 Excluded)", x = "Log Fold Change", y = "-Log10 FDR") +
  theme(plot.title = element_text(face = "bold"))

##GSEA##
gene_list<-results$table$logFC
names(gene_list)<-rownames(results$table)
gene_list<-sort(gene_list, decreasing = T)

gse <- gseGO(geneList = gene_list,
             ont = "ALL",
             keyType = "SYMBOL",
             minGSSize = 10,
             maxGSSize = 500,
             pvalueCutoff = 0.05,
             verbose=TRUE,
             OrgDb = org.Hs.eg.db )
s103dot<-dotplot(gse, showCategory=10, split=".sign") + ggtitle("SOT103 Excluded")

cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)

pairwisegse<-pairwise_termsim(gse)
emapplot(pairwisegse, showCategory = 30, repel=T)
ridgeplot(gse, label_format = 70) + labs(x = "Enrichment Distribution")

gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)


