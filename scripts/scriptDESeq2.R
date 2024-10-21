# Load DESeq2 library
library(DESeq2)

# Read in the sample information
sampleInfo <- read.csv("sample_info.txt", header = TRUE)

# Replace '/' in sample names to match count data
sampleInfo$Sample <- gsub("/", ".", sampleInfo$Sample)

# Read in the count data
countData <- read.table("FeatureCounts/all_featureCounts_default.txt", header = TRUE, row.names = 1, sep = "\t")

# Keep gene id and counts columns
countData <- as.matrix(countData[, !colnames(countData) %in% c("Chr", "Start", "End", "Strand", "Length")])

# Reorder sample names to match the sample data
sampleNames <- sampleInfo[,'Sample']
countData <- countData[,sampleInfo[,'Sample']]

# Convert variables to factors
sampleInfo$Zfp335 <- as.factor(sampleInfo$Zfp335)
sampleInfo$DietArm <- as.factor(sampleInfo$DietArm)
sampleInfo$Sex <- as.factor(sampleInfo$Sex)

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleInfo, design = ~ Zfp335 + DietArm + Sex)
dds$group <- factor(paste0(dds$Zfp335, dds$DietArm)) # add mutant-statin interaction
design(dds) <- ~ group + Sex # add sex covariate

# Perform DESeq2 normalization and differential expression analysis
dds <- DESeq(dds)

# View available contrasts
resultsNames(dds)

# Get differential expression results
res <- results(dds, contrast=c("group", "wildstatinArm", "mutantstatinArm"))
res

# Summarize the results
summary(res)

# Filter results with padj threshold
filtered_res <- res[!is.na(res$padj) & res$padj <= 0.1, ]
filtered_res

gene_ids <- rownames(filtered_res)
ensembl_ids <- gsub("\\..*", "", gene_ids)
write.table(ensembl_ids, file = "gene_list.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)


# VISUALIZATIONS
# Plot counts
plotCounts(dds, gene = which.min(res$padj), intgroup = "Zfp335")

# Visualize results
plotMA(res)

# Load package to normalize expression data
library(vsn)

# Perform transformation
vsd <- vst(dds, blind=FALSE)
meanSdPlot(assay(vsd))

filtered_vsd <- assay(vsd)[rownames(filtered_res),]

# Load package to map ensembl id's
library(biomaRt)

# Connect to the Ensembl database using the biomaRt package
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Query the conversion table for your gene list
conversion_result <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                           filters = "ensembl_gene_id", 
                           values = ensembl_ids, 
                           mart = ensembl)

# Use ensembl id if gene name not present
library(dplyr)

# Replace empty strings in "external_gene_name" with corresponding values from "ensembl_gene_id"
conversion_result <- conversion_result %>%
  mutate(external_gene_name = ifelse(external_gene_name == "", ensembl_gene_id, external_gene_name))

rownames(filtered_vsd) <- conversion_result$external_gene_name

# Load pheatmap library
library(pheatmap)

# Construct annotation dataframe
df <- as.data.frame(colData(dds)[,c("Zfp335","DietArm","Sex")])

# Plot heatmap
heatmap_plot <- pheatmap(filtered_vsd, cluster_rows = TRUE,
         show_rownames = TRUE, show_colnames = FALSE,
         cluster_cols = TRUE, annotation_col = df, )

# Save heatmap
library(ggplot2)
ggsave("heatmap_0.1.png", plot = heatmap_plot)

# Plot PCA
pcaData <- plotPCA(vsd, intgroup = c("Zfp335","Sex"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=pcaData$Zfp335, shape=pcaData$Sex)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

# Calculate Pearson correlation coefficients between significantly DEGs
cor_matrix <- cor(countData[rownames(filtered_res), ])

# Threshold the correlation matrix to identify significant interactions
cor_threshold <- 0.8  # Adjust as needed

# Convert correlation matrix to adjacency matrix
adj_matrix <- ifelse(abs(cor_matrix > cor_threshold), 1, 0)

library(igraph)
#Construct gene regulatory network
GRN <- graph_from_adjacency_matrix(adj_matrix, mode = "directed", weighted = TRUE)

# Visualize the network using igraph
plot(GRN, layout = layout_with_fr)
