raw_data <- read.table("combined_PUF60_Counts.txt", header = T)
colSums(raw_data)

raw_filter <- rowSums(raw_data > 60) >= 3
filtered_data <- raw_data[raw_filter, ]

plotMDS(filtered_data)
boxplot(filtered_data)

dge1 <- DGEList(counts = filtered_data, remove.zeros = TRUE)

dge1 <-calcNormFactors(dge1)
normalised_counts <- cpm(dge1)

plotMDS(dge1)

filtered_normalised_counts <- write.table(normalised_counts, "filtered_PUF60_counts.txt", sep = "\t", quote = FALSE, col.names = NA)
