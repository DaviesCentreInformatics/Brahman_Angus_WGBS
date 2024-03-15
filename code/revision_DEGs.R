# Load libraries
library(DESeq2)
library(Glimma)
library(edgeR)
library(pheatmap)
library(RColorBrewer)
library(vsn)
library(ggplot2)

rm(list=ls())
gc()
# Set working directory
getwd()
wd_path = "/Users/callummacphillamy/PhD/gigascience_revisions/WGBS_RNA"
setwd(wd_path)
assertthat::are_equal(getwd(), wd_path)

# Load count data
counts <- read.csv('./all.aligned2Angus.RNA.txt',
                   sep = '\t', row.names = 1)

# Rename the columns so they look nicer. Replace X with F
# Remove .sorted.bam from the end of the column names
colnames(counts) <- gsub('.sorted.bam', '', colnames(counts))
colnames(counts) <- gsub('X', 'F', colnames(counts))

head(counts)


# Load the samples
samples <- read.table("./sampleinfo.csv", sep=",",
                      header=TRUE,
                      stringsAsFactors = FALSE)
head(samples)

# Convert columns in the samples dataframe to factors
samples$Genetics <- as.factor(samples$Genetics)
samples$Sex <- as.factor(samples$Sex)
samples$Batch <- as.factor(samples$Batch)
samples$Sire <- as.factor(samples$Sire)
samples$Dam <- as.factor(samples$Dam)

# Ensure sample names in counts matrix match the order in the samples dataframe
table(colnames(counts) == samples$Sample)
colorder <- c(as.character(samples$Sample))
countsdata <- counts[, colorder]
table(colnames(countsdata) == samples$Sample)

# Create a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = countsdata,
                              colData = samples,
                              design = ~0 + Genetics + Sex + Batch)
dds_mat <- DESeqDataSetFromMatrix(countData = countsdata,
                                  colData = samples,
                                  design = ~0 + Dam + Sex + Batch)

dds_pat <- DESeqDataSetFromMatrix(countData = countsdata,
                                  colData = samples,
                                  design = ~0 + Sire + Sex + Batch)

# Transform the data for plotting
vsd <- vst(dds, blind=FALSE)
vsd_mat <- vst(dds_mat, blind=FALSE)
vsd_pat <- vst(dds_pat, blind=FALSE)

# Plot sample-to-sample distances
# sampleDists <- dist(t(assay(vsd)))
# sampleDistMatrix <- as.matrix(sampleDists)
# sampleDistMatrix
# rownames(sampleDistMatrix) <- paste(vsd$Genetics, vsd$Sex, vsd$Batch, sep="-")
# colnames(sampleDistMatrix) <- NULL
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# pheatmap(sampleDistMatrix,
#          clustering_distance_rows=sampleDists,
#          clustering_distance_cols=sampleDists,
#          col=colors)
# par(mfrow = c(1,2))
# pcaData <- plotPCA(vsd, intgroup=c("Genetics","Sex", "Batch"), returnData=TRUE)
# percentVar <- round(100 * attr(pcaData, "percentVar"))

# ggplot(pcaData, aes(PC1, PC2, color=Genetics, shape=Batch)) +
#     geom_point(size=3) +
#     xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#     ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#     coord_fixed() 

# glimmaMDS(dds, top=500, groups=dds$Genetics, 
#           main="MDS plot of top 500 genes")


# Pre-filtering
keep <- rowSums(counts(dds) >= 10) >= 3
dim(dds)
dds <- dds[keep, ]
dim(dds)

# Perform Differential Expression Analysis
dds <- DESeq(dds)

dim(subset(results(dds, contrast=c("Genetics", "BTBT", "BIBI" ), alpha=0.05), padj < 0.05))

res <- results(dds, contrast=c("Genetics", "BTBT", "BIBI"), alpha=0.05)
summary(res)

# Save only the significant genes
write.csv(subset(results(dds, contrast=c("Genetics", "BTBT", "BIBI"), alpha=0.05), padj < 0.05), file="BTvBI_DESeq2_aligned2Angus.csv")

# Pre-filtering
keep <- rowSums(counts(dds_mat) >= 10) >= 6
dim(dds_mat)
dds_mat <- dds_mat[keep, ]
dim(dds_mat)

# Perform Differential Expression Analysis
dds_mat <- DESeq(dds_mat)

dim(subset(results(dds_mat, contrast=c("Dam", "BT", "BI"), alpha=0.05), padj < 0.05))

res <- results(dds_mat, contrast=c("Dam", "BT", "BI"), alpha=0.05)
summary(res)

# Save only the significant genes
write.csv(subset(results(dds_mat, contrast=c("Dam", "BT", "BI"), alpha=0.05), padj < 0.05), file="maternal_BTvBI_DESeq2_aligned2Angus.csv")

# Pre-filtering
keep <- rowSums(counts(dds_pat) >= 10) >= 6
dim(dds_pat)
dds_mat <- dds_pat[keep, ]
dim(dds_pat)

# Perform Differential Expression Analysis
dds_pat <- DESeq(dds_pat)


dim(subset(results(dds_pat, contrast=c("Sire", "BT", "BI"), alpha=0.05), padj < 0.05))

res <- results(dds_pat, contrast=c("Sire", "BT", "BI"), alpha=0.05)
summary(res)

# Save only the significant genes
write.csv(subset(results(dds_pat, contrast=c("Sire", "BT", "BI"), alpha=0.05), padj < 0.05), file="paternal_BTvBI_DESeq2_aligned2Angus.csv")
