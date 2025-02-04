#file names
file_names <- c("counts/P001.out.tab", "counts/P002.out.tab", "counts/SOT103.out.tab", 
                "counts/SOT111.out.tab", "counts/SOT112.out.tab", "counts/SOT149.out.tab", 
                "counts/SOT211.out.tab", "counts/SOT325.out.tab", "counts/SOT344.out.tab")

#empty matrix
combined_Counts <- NULL  

# Loop through each file and extract count
for (file in file_names) {
  # Read the count file
  count_data <- read.table(file, header = FALSE, sep = "\t", row.names = NULL)
  
  # Extract count
  fourth_column <- count_data[, 4]  # Get the fourth column
  
  # Combine counts into the main matrix
  if (is.null(combined_Counts)) {
    combined_Counts <- fourth_column 
  } else {
    combined_Counts <- cbind(combined_Counts, fourth_column) 
  }
}

# Set row names to gene names
rownames(combined_Counts) <- count_data[, 1]

# Set column names 
colnames(combined_Counts) <- gsub(".out.tab", "", basename(file_names)) 

# Remove N_multimapping etc
unwanted_rows <- c("N_ambiguous", "N_multipmapping", "N_noFeature", "N_unmapped")
combined_Counts <- combined_Counts[!rownames(combined_Counts) %in% unwanted_rows, ]

combined_Counts <- combined_Counts[-1, ]

# Save 
write.table(combined_Counts, "combined_PUF60_Counts.txt", sep = "\t", quote = FALSE, col.names = NA)

 