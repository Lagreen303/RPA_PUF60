library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(pheatmap)
library(ggpubr)

#SKIPPING EVENTS
SE_rmats <- read.table('rmats/SE.MATS.JCEC.txt', header = T) %>%
  dplyr::select(., c('ID', 'geneSymbol', 'IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_2', 'FDR', 'IncLevel1', 'IncLevel2', 'IncLevelDifference'))

SE_rmats <- SE_rmats %>%
  #Split the replicate read counts that are separated by commas into different columns
  separate(., col = IJC_SAMPLE_1, into = c('IJC_P1', 'IJC_P2', 'IJC_S103', 'IJC_S149', 'IJC_S211', 'IJC_S325'), sep = ',', remove = T, convert = T) %>%
  separate(., col = SJC_SAMPLE_1, into = c('SJC_P1', 'SJC_P2', 'SJC_S103', 'SJC_S149', 'SJC_S211', 'SJC_S325'), sep = ',', remove = T, convert = T) %>%
  separate(., col = IJC_SAMPLE_2, into = c('IJC_S111', 'IJC_S112', 'IJC_S344'), sep = ',', remove = T, convert = T) %>%
  separate(., col = SJC_SAMPLE_2, into = c('SJC_S111', 'SJC_S112', 'SJC_S344'), sep = ',', remove = T, convert = T)

filtered_SE_rmats <- SE_rmats[abs(SE_rmats$IncLevelDifference) >= 0.2 & SE_rmats$FDR <= 0.05, ]

##ALTERNATIVE 3' 
A3_rmats <- read.table('rmats/A3SS.MATS.JCEC.txt', header = T) %>%
  dplyr::select(., c('ID', 'geneSymbol', 'IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_2', 'FDR', 'IncLevel1', 'IncLevel2', 'IncLevelDifference'))

A3_rmats <- A3_rmats %>%
  separate(., col = IJC_SAMPLE_1, into = c('IJC_P1', 'IJC_P2', 'IJC_S103', 'IJC_S149', 'IJC_S211', 'IJC_S325'), sep = ',', remove = T, convert = T) %>%
  separate(., col = SJC_SAMPLE_1, into = c('SJC_P1', 'SJC_P2', 'SJC_S103', 'SJC_S149', 'SJC_S211', 'SJC_S325'), sep = ',', remove = T, convert = T) %>%
  separate(., col = IJC_SAMPLE_2, into = c('IJC_S111', 'IJC_S112', 'IJC_S344'), sep = ',', remove = T, convert = T) %>%
  separate(., col = SJC_SAMPLE_2, into = c('SJC_S111', 'SJC_S112', 'SJC_S344'), sep = ',', remove = T, convert = T)

filtered_A3_rmats <- A3_rmats[abs(A3_rmats$IncLevelDifference) >= 0.2 & A3_rmats$FDR <= 0.05, ]

#ALTERNATIVE 5'
A5_rmats <- read.table('rmats/A5SS.MATS.JCEC.txt', header = T) %>%
  dplyr::select(., c('ID', 'geneSymbol', 'IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_2', 'FDR', 'IncLevel1', 'IncLevel2', 'IncLevelDifference'))

A5_rmats <- A5_rmats %>%
  separate(., col = IJC_SAMPLE_1, into = c('IJC_P1', 'IJC_P2', 'IJC_S103', 'IJC_S149', 'IJC_S211', 'IJC_S325'), sep = ',', remove = T, convert = T) %>%
  separate(., col = SJC_SAMPLE_1, into = c('SJC_P1', 'SJC_P2', 'SJC_S103', 'SJC_S149', 'SJC_S211', 'SJC_S325'), sep = ',', remove = T, convert = T) %>%
  separate(., col = IJC_SAMPLE_2, into = c('IJC_S111', 'IJC_S112', 'IJC_S344'), sep = ',', remove = T, convert = T) %>%
  separate(., col = SJC_SAMPLE_2, into = c('SJC_S111', 'SJC_S112', 'SJC_S344'), sep = ',', remove = T, convert = T)

filtered_A5_rmats <- A5_rmats[abs(A5_rmats$IncLevelDifference) >= 0.2 & A5_rmats$FDR <= 0.05, ]

#INTRON RETENTION
RI_rmats <- read.table('rmats/RI.MATS.JCEC.txt', header = T) %>%
  dplyr::select(., c('ID', 'geneSymbol', 'IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_2', 'FDR', 'IncLevel1', 'IncLevel2', 'IncLevelDifference'))

RI_rmats <- RI_rmats %>%
  separate(., col = IJC_SAMPLE_1, into = c('IJC_P1', 'IJC_P2', 'IJC_S103', 'IJC_S149', 'IJC_S211', 'IJC_S325'), sep = ',', remove = T, convert = T) %>%
  separate(., col = SJC_SAMPLE_1, into = c('SJC_P1', 'SJC_P2', 'SJC_S103', 'SJC_S149', 'SJC_S211', 'SJC_S325'), sep = ',', remove = T, convert = T) %>%
  separate(., col = IJC_SAMPLE_2, into = c('IJC_S111', 'IJC_S112', 'IJC_S344'), sep = ',', remove = T, convert = T) %>%
  separate(., col = SJC_SAMPLE_2, into = c('SJC_S111', 'SJC_S112', 'SJC_S344'), sep = ',', remove = T, convert = T)

filtered_RI_rmats <- RI_rmats[abs(RI_rmats$IncLevelDifference) >= 0.2 & RI_rmats$FDR <= 0.05, ]

#MUTUALLY EXCLUSIVE EXON
MXE_rmats <- read.table('rmats/MXE.MATS.JCEC.txt', header = T) %>%
  dplyr::select(., c('ID', 'geneSymbol', 'IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_2', 'FDR', 'IncLevel1', 'IncLevel2', 'IncLevelDifference'))

MXE_rmats <- MXE_rmats %>%
  separate(., col = IJC_SAMPLE_1, into = c('IJC_P1', 'IJC_P2', 'IJC_S103', 'IJC_S149', 'IJC_S211', 'IJC_S325'), sep = ',', remove = T, convert = T) %>%
  separate(., col = SJC_SAMPLE_1, into = c('SJC_P1', 'SJC_P2', 'SJC_S103', 'SJC_S149', 'SJC_S211', 'SJC_S325'), sep = ',', remove = T, convert = T) %>%
  separate(., col = IJC_SAMPLE_2, into = c('IJC_S111', 'IJC_S112', 'IJC_S344'), sep = ',', remove = T, convert = T) %>%
  separate(., col = SJC_SAMPLE_2, into = c('SJC_S111', 'SJC_S112', 'SJC_S344'), sep = ',', remove = T, convert = T)

filtered_MXE_rmats <- MXE_rmats[abs(MXE_rmats$IncLevelDifference) >= 0.2 & MXE_rmats$FDR <= 0.05, ]

#summary table
# Calculate the mean IncLevelDifference for each event
mean_A3 <- mean(filtered_A3_rmats$absIncLevelDifference, na.rm = TRUE)
mean_A5 <- mean(filtered_A5_rmats$IncLevelDifference, na.rm = TRUE)
mean_MXE <- mean(filtered_MXE_rmats$IncLevelDifference, na.rm = TRUE)
mean_RI <- mean(filtered_RI_rmats$IncLevelDifference, na.rm = TRUE)
mean_SE <- mean(filtered_SE_rmats$IncLevelDifference, na.rm = TRUE)

#Number of events
num_A3 <- nrow(filtered_A3_rmats)
num_A5 <- nrow(filtered_A5_rmats)
num_MXE <- nrow(filtered_MXE_rmats)
num_RI <- nrow(filtered_RI_rmats)
num_SE <- nrow(filtered_SE_rmats)

#proportion of increased events
# Calculate the proportion of positive to negative events
perc_pos_A3 <- (sum(filtered_A3_rmats$IncLevelDifference > 0) / num_A3) * 100
perc_pos_A5 <- (sum(filtered_A5_rmats$IncLevelDifference > 0) / num_A5) * 100
perc_pos_MXE <- (sum(filtered_MXE_rmats$IncLevelDifference > 0) / num_MXE) * 100
perc_pos_RI <- (sum(filtered_RI_rmats$IncLevelDifference > 0) / num_RI) * 100
perc_pos_SE <- (sum(filtered_SE_rmats$IncLevelDifference > 0) / num_SE) * 100

# Create a summary table
summary_table <- data.frame(
  EventType = c("A3", "A5", "MXE", "RI", "SE"),
  EventNum = c(num_A3, num_A5, num_MXE, num_RI, num_SE),
  MeanIncLevelDifference = c(mean_A3, mean_A5, mean_MXE, mean_RI, mean_SE),
  PercPosEvents = c(perc_pos_A3, perc_pos_A5, perc_pos_MXE, perc_pos_RI, perc_pos_SE)
)

#Barcharts
filtered_A3_rmats$EventID <- paste(filtered_A3_rmats$ID, filtered_A3_rmats$geneSymbol, sep="_")
A3bar <- ggplot(filtered_A3_rmats, 
       aes(x=EventID, y=IncLevelDifference, fill=EventID)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "A3 Events",
       x = "Splicing Event",
       y = "Inclusion Level Difference") +
  theme_bw() + #for blank use theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) | For IDs on charts use theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  guides(fill = FALSE)

filtered_A5_rmats$EventID <- paste(filtered_A5_rmats$ID, filtered_A5_rmats$geneSymbol, sep="_")
A5bar<-ggplot(filtered_A5_rmats, 
       aes(x=EventID, y=IncLevelDifference, fill=EventID)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "A5 Events",
       x = "Splicing Event",
       y = "Inclusion Level Difference") +
  theme_bw() + #for blank use theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) | For IDs on charts use theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  guides(fill = FALSE)

filtered_MXE_rmats$EventID <- paste(filtered_MXE_rmats$ID, filtered_MXE_rmats$geneSymbol, sep="_")
MXEbar<-ggplot(filtered_MXE_rmats, 
       aes(x=EventID, y=IncLevelDifference, fill=EventID)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "MXE Events",
       x = "Splicing Event",
       y = "Inclusion Level Difference") +
  theme_bw() + #for blank use theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) | For IDs on charts use theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  guides(fill = FALSE)

filtered_RI_rmats$EventID <- paste(filtered_RI_rmats$ID, filtered_RI_rmats$geneSymbol, sep="_")
RIbar<-ggplot(filtered_RI_rmats, 
       aes(x=EventID, y=IncLevelDifference, fill=EventID)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "RI Events",
       x = "Splicing Event",
       y = "Inclusion Level Difference") +
  theme_bw() + #for blank use theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) | For IDs on charts use theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  guides(fill = FALSE)

filtered_SE_rmats$EventID <- paste(filtered_SE_rmats$ID, filtered_SE_rmats$geneSymbol, sep="_")
SEbar<-ggplot(filtered_SE_rmats, 
       aes(x=EventID, y=IncLevelDifference, fill=EventID)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "SE Events",
       x = "Splicing Event",
       y = "Inclusion Level Difference") +
  theme_bw() + #for blank use theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) | For IDs on charts use theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  guides(fill = FALSE)

ggarrange(A5bar, A3bar, MXEbar, RIbar, SEbar + rremove("x.text"),
          labels = c("A", "B", "C", "D", "E"),
          ncol = 3, nrow = 2)

ggarrange(Alldot, s103dot, p1s103dot,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)

