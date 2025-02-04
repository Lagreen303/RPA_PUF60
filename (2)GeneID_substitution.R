#Script to swap Gene IDs for Gene names
# Load necessary libraries
library(dplyr)

# Load your count matrix from a text file
count_matrix <- read.table("filtered_PUF60_counts.txt", header = TRUE, sep = "\t", row.names = 1)

# Load the mapping file with Ensembl IDs and Gene Names
mapping_file <- read.table("gene_names", header = FALSE, sep = "\t") 
colnames(mapping_file) <- c("EnsemblID", "GeneName")

count_df <- as.data.frame(count_matrix)

# Add a column for Ensembl IDs
count_df$EnsemblID <- rownames(count_df)

# Merge the count data with the mapping file
merged_df <- merge(count_df, mapping_file, by = "EnsemblID", all.x = TRUE)

# Remove the EnsemblID column and keep only GeneNames and counts
final_counts <- merged_df[, c("GeneName", colnames(count_df)[-ncol(count_df)])]
write.table(final_counts, "genename_counts.txt", sep = "\t", row.names = T)

rownames(final_counts) <- final_counts$GeneName

cc<- read.table("genename_counts.txt")

######################################
# Load necessary libraries
library(dplyr)

# Load your count matrix from a text file
count_matrix <- read.table("filtered_PUF60_counts.txt", header = TRUE, sep = "\t", row.names = 1)

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
write.table(final_counts, "genename_counts.txt", sep = "\t", row.names = TRUE)
