# Load  libraries
library(edgeR)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(reshape2)

# Load your normalized count matrix from a text file
puf60count <- read.table("count_files/genename_counts.txt", header = TRUE, sep = "\t")

# Load your metadata
puf60_group <- ifelse(colnames(puf60count) %in% c("P002", "P001", "SOT103", "SOT149", "SOT211", "SOT325"), 
                      "PUF60 Variant", "No Variant")

# Create a new metadata data frame if needed
puf60.metaData <- data.frame(Sample = colnames(puf60count), Group = puf60_group)

# Set up grouping information from the metadata
group <- factor(puf60.metaData$Group)

# Create a DGEList object with a different name
dge_list <- DGEList(counts = puf60count, group = group) 

# Normalize and filter (check you have not already done this)
dge_list <- calcNormFactors(dge_list)
keep <- filterByExpr(dge_list)
dge_list <- dge_list[keep, ]

# Estimate dispersion
dge_list <- estimateDisp(dge_list)

# Fit the model
fit <- glmFit(dge_list)

# Perform likelihood ratio test
lrt <- glmLRT(fit)

# Extract results
results <- topTags(lrt, n = Inf)
significant_genes <- results$table[abs(results$table$logFC) >2& results$table$PValue < 0.05 & results$table$FDR < 0.05, ]  

#Sample Correlation heatmap (Change correlation power (44) as needed for visualisation)
melt_corr <-cor(puf60count) |>
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
results$table$diffexpressed[results$table$logFC > -2 & results$table$FDR <0.05] <- "NO"

#label significant genes
results$table$genelabel <- ""
results$table$genelabel <- ifelse(rownames(results$table) == "ENSG00000276241"
                                  |rownames(results$table) == "RPL18AP3"
                                  |rownames(results$table) == "OVCH1-AS1"
                                  |rownames(results$table) == "VWCE"
                                  |rownames(results$table) == "TPRG1"
                                  |rownames(results$table) == "TMSB4XP4", rownames(results$table), "")

#Plot volcano plot
ggplot(as.data.frame(results$table)) +
  geom_point(aes(x = logFC, y = -log10(FDR), color = diffexpressed), size=2) +
  geom_text_repel(aes(label = results$table$genelabel, max.overlaps = 1, x = logFC, y = -log10(FDR))) +
  scale_colour_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-2, 2), linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed')+
  theme_classic() +
  theme(legend.position = "none")+
  labs(title = "PUF60 Differential Gene Expression Volcano Plot", x = "Log Fold Change", y = "-Log10 FDR")

# Heatmap of Significant Genes
heatmap_data <- dge_list$counts[rownames(significant_genes), ]
heatmap_data <- cpm(dge_list, normalized.lib.sizes = TRUE)  
pheatmap(heatmap_data, scale = "row", clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
         clustering_method = "complete", show_rownames = F, show_colnames = T, main = "Heatmap of Gene Expression Across Samples", treeheight_row = 0, treeheight_col = 50)

