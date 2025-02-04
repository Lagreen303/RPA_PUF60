#PCA
library(tidyverse)
library(factoextra)

puf60count <- read.table("count_files/genename_counts.txt", header = TRUE, row.names=1)
puf60.svd <- 
  puf60count %>% 
  t() %>% 
  prcomp(scale = T)

summary(puf60.svd)



fviz_eig(puf60.svd, addlabels = T)
fviz_pca_ind(puf60.svd, repel = T)

fviz_pca_var(puf60.svd, select.var = list(cos2 = 10), repel = T)


# Define the groups based on your sample names
puf60_group <- ifelse(colnames(puf60count) %in% c("P002", "P001", "SOT103", "SOT149", "SOT211", "SOT325"), 
                      "PUF60 Variant", "No Variant")

# Create a new metadata data frame if needed
puf60.metaData <- data.frame(Sample = colnames(puf60count), Group = puf60_group)

# Generate the PCA biplot
fviz_pca_biplot(puf60.svd, 
                select.var = list(name = c("TMEM1", "NOP1")),
                habillage = puf60.metaData$Group, 
                label = "all", 
                repel = T,
                invisible = "quali",
                addEllipses = F) + 
  ggtitle("PUF60 Group Separation")

# Generate the PCA biplot with bigger font and red and blue colors
fviz_pca_biplot(puf60.svd, 
                select.var = list(name = c("TMEM1", "NOP1")),
                habillage = puf60.metaData$Group, 
                label = "all", 
                repel = TRUE,
                invisible = "quali",
                addEllipses = FALSE) + 
  ggtitle("PUF60 Group Separation") +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1.5, "lines")
  ) +
  scale_color_manual(values = c("red", "blue"))
