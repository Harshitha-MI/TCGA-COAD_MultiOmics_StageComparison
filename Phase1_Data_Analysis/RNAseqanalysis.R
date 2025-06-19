
# Load libraries
library(SummarizedExperiment)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(tibble)
library(EnhancedVolcano)
library(dplyr)


# Load the data
rna_data <- readRDS("RNAseq_STARCounts.rds")

# Extract unstranded counts and gene names
counts <- assay(rna_data, "unstranded")
gene_symbols <- rowData(rna_data)$gene_name
# Changing ENsemble id to gene names
rownames(counts) <- gene_symbols
#filter out NA gene names and duplicate gene names
counts <- counts[!is.na(rownames(counts)) & !duplicated(rownames(counts)), ]

# Extract patient IDs (TCGA-XX-YYYY)
patient_ids <- substr(colnames(counts), 1, 12)
table(patient_ids)
length(patient_ids)
colData_rna <- colData(rna_data)
table(colData_rna$sample_type)

metadata <- colData_rna[, c("barcode", "patient", "sample_type", "ajcc_pathologic_stage")]

# Standardize staging info
metadata$ajcc_pathologic_stage <- as.character(metadata$ajcc_pathologic_stage)

# Convert subtypes like "Stage IA", "Stage IB" → "StageI", same for Stage IV
metadata$ajcc_pathologic_stage_grouped <- ifelse(
  grepl("^Stage IV", metadata$ajcc_pathologic_stage), "StageIV",
  ifelse(grepl("^Stage I", metadata$ajcc_pathologic_stage), "StageI", NA)
)
# Filter metadata to only tumor samples in Stage I and IV
metadata_filtered <- metadata[metadata$sample_type != "Solid Tissue Normal" &
                              metadata$ajcc_pathologic_stage_grouped %in% c("StageI","StageIV"),]

colnames(metadata_filtered)[colnames(metadata_filtered) == "barcode"] <- "SampleID"
colnames(metadata_filtered)[colnames(metadata_filtered) == "sample_type"] <- "SampleType"
counts_filtered <- counts[, metadata_filtered$SampleID]

# Create DESeq2 object using the grouped stage info
dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = metadata_filtered,
                              design = ~ ajcc_pathologic_stage_grouped)

# performing filteration: removing low expression genes
keep <- rowSums(counts(dds) >= 10) >= (0.2 * ncol(dds))  # at least 10 counts in at least 20% of samples
dds <- dds[keep, ]
dds$ajcc_pathologic_stage_grouped <- relevel(dds$ajcc_pathologic_stage_grouped, ref = "StageI")
dds <- DESeq(dds)
dds


res <- results(dds, contrast = c("ajcc_pathologic_stage_grouped", "StageIV", "StageI"))
colSums(is.na(res))
head(res)
# Remove rows with NA in log2FoldChange or padj
res_cleaned <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]

# Filter for significant DE genes
res_filtered <- res_cleaned[abs(res_cleaned$log2FoldChange) >= 1 & res_cleaned$padj < 0.05, ]

# Separate upregulated and downregulated
upregulated_genes <- res_filtered[res_filtered$log2FoldChange > 0, ]
downregulated_genes <- res_filtered[res_filtered$log2FoldChange < 0, ]

# Print counts
cat("Upregulated genes: ", nrow(upregulated_genes), "\n")
cat("Downregulated genes: ", nrow(downregulated_genes), "\n")
write.csv(as.data.frame(res_filtered), "DEG_StageIV_vs_StageI.csv")

EnhancedVolcano(res,
                lab = rownames(res), # Gene labels
                x = 'log2FoldChange', # X-axis: log2 fold change
                y = 'padj', # Y-axis: adjusted p-value (FDR)
                pCutoff = 0.05, # FDR cutoff of 0.05
                FCcutoff = 1, # Fold-change cutoff of 4 (log2FC = 2)
                pointSize = 2.0, # Size of points
                labSize = 3.0, # Size of labels
                col = c('lightskyblue4', 'steelblue1', 'rosybrown1', 'hotpink1'), # Custom colors
                title = "Volcano Plot", # Plot title
                subtitle = "FDR < 0.05, |log2FC| > 2" # Plot subtitle
)


