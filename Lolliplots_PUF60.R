# Loliplots

# load libraries
library(trackViewer)
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)

# import gtf file
gtf_data <- import.gff("gencode.v47.basic.annotation.gtf.gz")

# extract exon information from gtf file
p <- gtf_data$gene_name=="PUF60" & gtf_data$type =="exon" & gtf_data$transcript_id=="ENST00000526683.6"
puf60<- gtf_data[p,]
puf60

# create GRanges object with the features we need (exons and protein domains)
features <- GRanges("chr8", IRanges(c(143818058,143817688,143816553), width = c(440,315,263), names=c("RRM1", "RMM2", "RMM3")))

puf60_dom <- append(puf60, features)
puf60_dom$featureLayerID <- rep(c("Exons","domains"), times=c(12,3))
puf60_dom$fill <- rep(c("#FF8833", "#51C6E6"), times=c(12,3))
n <- rep("Exons", times=12)
names <- c(n, "RRM1", "RRM2", "RRM3")
names(puf60_dom) <- names

puf60_dom$height[13:15] <- list(unit(0.5, "lines"), unit(0.5, "lines"), unit(0.5, "lines"))

# input variant data
SNP <- c(143816605,143818077,143818236,143821817,143817874, 143817750)
# add variant annotation
sample.gr <- GRanges("chr8", IRanges(SNP, width=1, names=c("c.1594del (P001)", "c.604-2A (SOT149)", "c.560T>A (SOT211)", "c.207+1G>T (SOT103)", "c.805_806del (SOT325)", "c.850dup (P002)")))
sample_col<- c("blue","red","red", "red", "blue","blue" )
sample.gr$color <- sample_col
sample.gr$border <- sample(c("black", "black"), length(SNP), replace=TRUE)

legend <- list(
  labels = c("In-frame", "Frameshift"),
  fill = c("red", "blue"),  
  border = "black",
  cex = 1.2,  # Increase the text size
  pch = 21
)
# make the lolliplot
lolliplot(sample.gr, puf60_dom, legend=legend, legend.pos="bottom")
lolliplot(sample.gr, puf60_dom, legend=legend)

