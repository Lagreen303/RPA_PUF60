library(OUTRIDER)

ods<-readRDS("carolina/outmodel_DRAGEN")
highlight_samples <- c("P001", "P002", "SOT149", "SOT103", "SOT211", "SOT325")
control_samples <- c("SOT111", "SOT112", "SOT344")

plotExpressionRank(ods, "PUF60", basePlot = , groups = highlight_samples, log = TRUE)

ods <- counts(ods, normalized = TRUE, minE = 0.5)
plot_gene_counts(ods, gene_of_interest = "PUF60", samples_to_highlight = highlight_samples )

counts_matrix <- ods@assays@data$counts


#Extract PUF60 results (change path to outrider dataset)
library(dplyr)
ods<-readRDS("carolina/outmodel_DRAGEN")
res <- results(ods, all = TRUE)
puf60res <- res %>% 
  filter((geneID == "PUF60"))

samplepuf60res <- puf60res %>% 
  filter(sampleID %in% c(highlight_samples, control_samples)) %>% 
  select(sampleID, zScore, normcounts, meanCorrected)

write.csv(samplepuf60res, "PUF60_Zscores.csv")
