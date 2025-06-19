# installing required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("vroom")
BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment")
BiocManager::install(c("sesame", "sesameData"))

# Load TCGAbiolinks
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(sesameData)
library(sesame)
library(vroom)

# Define project
project_id <- "TCGA-COAD"

# ===============================
# FUNCTIONS
# ===============================

# Fetch clinical data
clinical <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical")

# Create a named vector: names = barcodes, values = stage labels
barcode_labels <- clinical$submitter_id
names(barcode_labels) <- clinical$ajcc_pathologic_stage
head(barcode_labels)

clinical$stage_group <- case_when(
  grepl("^Stage IV", clinical$ajcc_pathologic_stage, ignore.case = TRUE) ~ "StageIV",
  grepl("^Stage III", clinical$ajcc_pathologic_stage, ignore.case = TRUE) ~ "StageIII",
  grepl("^Stage II", clinical$ajcc_pathologic_stage, ignore.case = TRUE) ~ "StageII",
  grepl("^Stage I", clinical$ajcc_pathologic_stage, ignore.case = TRUE) ~ "StageI",
  TRUE ~ clinical$ajcc_pathologic_stage
)
unique(clinical$stage_group)
table(clinical$stage_group, useNA="ifany")

# Function to get data.types
get_data_types <- function(project_id, category) {
  query <- GDCquery(project = project_id, data.category = category)
  unique(query$results[[1]]$data_type)
}

get_data_types(project_id, "Proteome Profiling")

# Function to get barcodes
get_barcodes <- function(project_id, category, type, workflow = NULL, platform = NULL, clinical_data = NULL, stage_filter = NULL) {
  params <- list(
    project = project_id,
    data.category = category,
    data.type = type
  )
  
  # Add workflow.type only if not NULL
  if (!is.null(workflow)) {
    params$workflow.type <- workflow
  }
  
  # Add platform only if not NULL
  if (!is.null(platform)) {
    params$platform <- platform
  }
  
  # Call GDCquery with parameters dynamically
  query <- do.call(GDCquery, params)
  
  # Extract all barcodes (patient IDs)
  all_barcodes <- unique(substr(query$results[[1]]$cases, 1, 12))
  
  # If no clinical data or no filter specified, return all barcodes
  if (is.null(clinical_data) || is.null(stage_filter)) {
    return(all_barcodes)
  }
  
  # Filter clinical data for desired stage(s)
  filtered_clinical <- clinical_data[clinical_data$stage_group %in% stage_filter, ]
  
  # Get barcodes from filtered clinical data
  filtered_barcodes <- filtered_clinical$submitter_id
  
  # Return barcodes present in both query and filtered clinical
  intersect(all_barcodes, filtered_barcodes)
}

getProjectSummary(project_id)


# Get barcodes
rna_satgeI_cases <- get_barcodes(project_id,"Transcriptome Profiling", "Gene Expression Quantification", "STAR - Counts",
                          clinical_data = clinical, stage_filter = "StageI")
rna_satgeIV_cases <- get_barcodes(project_id,"Transcriptome Profiling", "Gene Expression Quantification", "STAR - Counts",
                                  clinical_data = clinical, stage_filter = "StageIV")
maf_stageI_cases <- get_barcodes(project_id, "Simple Nucleotide Variation", "Masked Somatic Mutation",
                          clinical_data = clinical, stage_filter = "StageI")
maf_stageIV_cases <- get_barcodes(project_id, "Simple Nucleotide Variation", "Masked Somatic Mutation",
                                  clinical_data = clinical, stage_filter = "StageIV")
meth_stageI_cases <- get_barcodes(project_id, "DNA Methylation", "Methylation Beta Value", 
                                   platform = "Illumina Human Methylation 450",
                           clinical_data = clinical, stage_filter = "StageI")
meth_stageIV_cases <- get_barcodes(project_id, "DNA Methylation", "Methylation Beta Value",
                           platform = "Illumina Human Methylation 450",
                           clinical_data = clinical, stage_filter = "StageIV")
cnv_satgeI_cases <- get_barcodes(project_id,"Copy Number Variation", "Copy Number Segment", 
                                   clinical_data = clinical, stage_filter = "StageI")
cnv_satgeIV_cases <- get_barcodes(project_id,"Copy Number Variation", "Copy Number Segment", 
                          clinical_data = clinical, stage_filter = "StageIV")
proteome_stageI_cases <- get_barcodes(project_id, "Proteome Profiling", "Protein Expression Quantification", 
                                       clinical_data = clinical, stage_filter = "StageI")
proteome_stageIV_cases <- get_barcodes(project_id, "Proteome Profiling", "Protein Expression Quantification", 
                                        clinical_data = clinical, stage_filter = "StageIV")


# Intersect common cases
common_stageI_cases <- Reduce(intersect, list(
  rna_satgeI_cases, 
  maf_stageI_cases, 
  meth_stageI_cases, 
  cnv_satgeI_cases,
  proteome_stageI_cases
))

length(common_stageI_cases)  # Number of Stage I patients with all omics
head(common_stageI_cases)

common_stageIV_cases <- Reduce(intersect, list(
  rna_satgeIV_cases, 
  maf_stageIV_cases, 
  meth_stageIV_cases, 
  cnv_satgeIV_cases,
  proteome_stageIV_cases
))

length(common_stageIV_cases)  # Number of Stage III patients with all omics
head(common_stageIV_cases)



#set.seed(123)  # for reproducibility

#subset_stageII <- sample(common_stageII_cases, 20)
#subset_stageIII <- sample(common_stageIII_cases, 20)

# Check the subsets
#length(subset_stageII)  # should be 20
#length(subset_stageIII) # should be 20

#head(subset_stageII)
#head(subset_stageIII)

# Combine all barcodes you want data for
all_barcodes <- c(common_stageI_cases, common_stageIV_cases)

# ===============================
# 1. GENE EXPRESSION (STAR - Counts)
# ===============================
query_rna <- GDCquery(project = project_id,
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts",
                      barcode = all_barcodes)

GDCdownload(query_rna)
rna_data <- GDCprepare(query_rna)

# Save for reuse
saveRDS(rna_data, file = "RNAseq_STARCounts.rds")

# ===============================
# 2. MUTATION DATA (Masked Somatic MAF)
# ===============================
query_maf <- GDCquery(project = project_id,
                      data.category = "Simple Nucleotide Variation",
                      data.type = "Masked Somatic Mutation",
                      barcode = all_barcodes)

GDCdownload(query_maf)
maf_data <- GDCprepare(query_maf)

# Save MAF
saveRDS(maf_data, file = "Mutation_MAF.rds")

# ===============================
# 3. DNA METHYLATION (450K)
# ===============================
query_meth <- GDCquery(project = project_id,
                       data.category = "DNA Methylation",
                       data.type = "Methylation Beta Value",
                       platform = "Illumina Human Methylation 450",
                       barcode = all_barcodes)

GDCdownload(query_meth, files.per.chunk = 42)
#GDCdownload(query_meth)
meth_data <- GDCprepare(query_meth)

# Save methylation
saveRDS(meth_data, file = "Methylation_450K.rds")

# ===============================
# 4. CLINICAL DATA
# ===============================
clinical_data_subset <- clinical[clinical$submitter_id %in% all_barcodes, ]
saveRDS(clinical_data_subset, file = "Clinical_Data_Subset.rds")

# ===============================
# 5. COPY NUMBER VARIATION (Optional)
# ===============================
query_cnv <- GDCquery(project = project_id,
                      data.category = "Copy Number Variation",
                      data.type = "Copy Number Segment",
                      barcode = all_barcodes)

GDCdownload(query_cnv)
cnv_data <- GDCprepare(query_cnv)

saveRDS(cnv_data, file = "CopyNumber_Segments.rds")

# ===============================
# 6. PROTEOMICS 
# ===============================
query_proteome <- GDCquery(project = project_id,
                           data.category = "Proteome Profiling", 
                           data.type = "Protein Expression Quantification", 
                           barcode = all_barcodes)

GDCdownload(query_proteome)
proteome_data <- GDCprepare(query_proteome)

saveRDS(proteome_data, file = "Proteome_exp.rds")

# ===============================
# Done!
# ===============================
message("✅ All key datasets downloaded. Phase 0 complete!")


# ===============================
# 1. ERROR FOR BARCODES FUNCTION
# ===============================
# -- Helper function
#get_barcodes <- function(project_id, category, type, workflow = NULL, platform = NULL) {
#query <- GDCquery(
#project = project_id,
#data.category = category,
#data.type = type,
#workflow.type = workflow,
#platform = platform
# )
#unique(substr(query$results[[1]]$cases, 1, 12))  # First 12 characters = patient ID
#}
#rna <- GDCquery(project = project_id,
#data.category = "Transcriptome Profiling", 
#data.type = "Gene Expression Quantification",
#workflow.type = "STAR - Counts" )

#rna_barcodes <- unique(substr(rna$results[[1]]$cases, 1, 12))
#length(rna_barcodes)

#clinical$stage_group <- case_when(
#grepl("^Stage II", clinical$ajcc_pathologic_stage, ignore.case = TRUE) ~ "Stage II",
# grepl("^Stage III", clinical$ajcc_pathologic_stage, ignore.case = TRUE) ~ "Stage III",
#TRUE ~ clinical$ajcc_pathologic_stage
#)


#clinical$stage_group <- ifelse(
#grepl("^Stage II", clinical$ajcc_pathologic_stage, ignore.case = TRUE), "Stage II",                      
#ifelse(
#grepl("^Stage III", clinical$ajcc_pathologic_stage, ignore.case = TRUE), "Stage III",                 
#clinical$ajcc_pathologic_stage                                                
#  )
#)

#unique(clinical$stage_group)

#class(clinical$ajcc_pathologic_stage)
#class(clinical$stage_group)

#clinical[grepl("^Stage II", clinical$ajcc_pathologic_stage, ignore.case = TRUE), c("ajcc_pathologic_stage", "stage_group")]

#unique_clinical <- unique(clinical$ajcc_pathologic_stage)
#print(unique_clinical)

#table(grepl("^Stage III", clinical$ajcc_pathologic_stage, ignore.case = TRUE))

# ===============================
# END OF DEBUG STEPS
# ===============================
