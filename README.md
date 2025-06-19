# 🧬 TCGA-COAD Multi-Omics Stage Comparison

This project explores multi-omics integration in colorectal cancer using TCGA-COAD data. It focuses on comparing **Stage I** (early-stage) and **Stage IV** (late-stage) tumor samples to identify molecular changes across multiple omics layers.

---

## 🔍 Objective

To develop a reproducible pipeline that integrates RNA-seq, somatic mutation, methylation, CNV, and proteomics data — enabling exploration of tumor progression in colorectal cancer.

---

## 📊 Data Modalities Used

- RNA-seq (STAR Counts)
- Somatic Mutations (Masked MAF)
- DNA Methylation (450K)
- Copy Number Variation (CNV Segments)
- Proteomics (Protein Expression)
- Clinical metadata

---

## ⚙️ Phases of Analysis

### ✅ **Phase 0: Data Acquisition**
- Retrieved clinical metadata from TCGA-COAD and grouped samples into Stage I–IV.
- Defined a barcode retrieval function per omics type and filtered for patients with complete data.
- Downloaded and saved:
  - RNA-seq
  - MAF (somatic mutation)
  - Methylation (450K)
  - CNV
  - Proteomics
  - Clinical data

### ✅ **Phase 1: RNA-seq Analysis**
- Processed STAR Count data.
- Filtered and normalized using `DESeq2`.
- Performed differential expression between Stage I and Stage IV tumor samples.
- Identified significantly up- and down-regulated genes and exported results.

---

More phases to come as additional omics layers are explored and integrated.


