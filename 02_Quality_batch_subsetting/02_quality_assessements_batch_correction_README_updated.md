#  **02_Quality Assessment and Batch Correction Module**
### *Genome Quality, Coverage Evaluation, and Data Harmonization for Biosurveillance Analysis*

This module performs **quality evaluation and harmonization** of genome-resolved metagenomic data within the GRUMB framework.  
It ensures that all abundance matrices, assemblies, and MAG-level outputs are accurate, comparable, and ready for downstream ecological and machine learning analyses.

---

## **Overview**

This module automates:

- Construction of species-level TPM matrices for batch correction  
- CLR transformation and batch correction using **limma**  
- Genome and assembly quality evaluation   
- Host contamination screening with 
- Runtime benchmarking and reproducibility profiling  
- Multivariate diagnostics (PCA, UMAP, t-SNE, PERMANOVA) pre- and post-correction  

---

##  **Directory Structure**
```

02_quality_assessements_batch_correction/
├── 01_Species_TPM matrix_for_limma.py # Builds species-level TPM matrix
├── 02_CLRtransformation_batch_correction.R # CLR transform + batch correction (limma)
├── 03_final_collect_files_downstream.sh # Merges QC outputs and batch-corrected matrices
├── 04_fastqscreen.py # Summarizes host contamination reports
├── 05_mag_quality.py # Summarizes MAG completeness & contamination
├── 06_quast.py # Aggregates QUAST assembly metrics
├── 07_runtime_analytics.py # Profiles runtime and compute performance
└── README.md
```

---

## **Workflow Summary**

### **1 Build Species-Level TPM Matrix for Batch Correction**
**Script:** `01_Species_TPM matrix_for_limma.py`

- Compiles **TPM abundance data** across all samples and environments into a single unified matrix.  
- Ensures correct sample–metadata alignment and consistent header formatting.  
- Exports both `.csv` and `.tsv` versions for R-based processing.  

*Output:**  
- `TPM_matrix.csv`  
- `TPM_matrix.tsv`  

---

### **2️ CLR Transformation and Batch Correction**
**Script:** `02_CLRtransformation_batch_correction.R`  

- Performs **Centered Log-Ratio (CLR)** transformation to normalize compositional TPM data.  
- Applies **limma::removeBatchEffect()** to remove technical variation while preserving biological structure.  
- Two correction strategies are tested:
  - **Version 1 (v1):** Corrects for *Instrument* and *Sequencing Center*  
  - **Version 2 (v2):** Adds *Project_ID* as a batch term  
- Diagnostic plots before and after correction include:
  - **PCA**, **UMAP**, and **t-SNE** embeddings  
  - **PERMANOVA (R²)** comparisons for Group and technical factors  
- Automatically exports high-resolution publication-ready figures.  

*Output Files:**
- `TPM_clr.csv`  
- `TPM_clr_batch_corrected_v1.csv`  
- `TPM_clr_batch_corrected_v2.csv`  
- `PCA_Before_Correction.png`, `UMAP_Before_Correction.png`, `tSNE_Before_Correction.png`  
- `PERMANOVA_before.csv`, `PERMANOVA_after_v1.csv`, `PERMANOVA_after_v2.csv`  
- `PERMANOVA_R2_comparison_fixed.png`, `Barplot_variance_explained.tiff`  
- `tSNE_After_Batch_Correction_v1.tiff`  

**Purpose:**  
Ensures batch harmonization across diverse metagenomic datasets while preserving true environmental differences.

---

### **3️ Final Collection of Corrected and QC Files**
**Script:** `03_final_collect_files_downstream.sh`  ### 

- These data are part of the bioinformatic outputs.  
- They are use visualization of data quality
- These outputs are used for the subsequent scripts for quality metrics evaluation


---

### **4 Host Contamination Evaluation**
**Script:** `04_fastqscreen.py`

- Summarizes **FastQ Screen** output across samples.  
- Quantifies read proportions mapping to *human*, *mouse*, *Plasmodium*, and *vector* genomes.  
- Produces contamination profiles per environment for QC reporting.  

**Output:**  
- `fastqscreen_summary.csv`  
- `fastqscreen_barplots.png`

---

### **5️ MAG Quality Assessment**
**Script:** `05_mag_quality.py`

- Parses **CheckM2** outputs for completeness and contamination metrics.  
- Categorizes MAGs into **high-quality (≥90%)** and **medium-quality (≥50%)** bins.  
- Generates per-environment summaries of genome quality.  

**Output:**  
- `mag_quality_summary.csv`  
- `MAG_quality_boxplots.png`

---

### **6️ Assembly Quality Evaluation**
**Script:** `06_quast.py`

- Aggregates and visualizes **QUAST** reports for assemblies.  
- Computes N50, L50, total contig length, and GC% per environment.  
- Performs correlation analysis (Pearson/Spearman) between metrics.  
- Generates comparative boxplots and scatterplots.  

**Output:**  
- `quast_merged_filtered.tsv`  
- `Quast_boxplots_Filtered.png`, `n50_vs_contigs_Filtered.png`, `gc_vs_n50_scatter_Filtered.png`

---

### **7️ Runtime Analytics**
**Script:** `07_runtime_analytics.py`

- Tracks average and variance in runtime across major pipeline modules (QC → Annotation).  
- Supports reproducibility and scalability documentation.  

**Output:**  
- `stepwise_runtime_per_environment.csv`  
- `stepwise_runtime_global.csv`

 **Example Runtime Summary (from Supplementary Table 5)**

| **Workflow Stage** | **Avg Runtime (min)** | **Std Dev (±)** | **Notes** |
|--------------------|----------------------:|----------------:|-----------|
| Preprocessing & QC | 8.6 | 22.3 | Trimming, filtering, host removal |
| Assembly & Binning | 22.4 | 45.5 | Scales with input complexity |
| Dereplication & QC | 18.7 | 38.9 | dRep + CheckM2 |
| Annotation (CARD/VFDB) | 9.8 | 17.3 | Functional annotation |
| Output Aggregation | 5.2 | 10.8 | Merging and cleanup |

---

## **Dependencies**

| Tool / Library | Purpose |
|----------------|----------|
| **Python ≥ 3.10** | Data wrangling and QC summaries |
| **R ≥ 4.2** | CLR transformation and limma batch correction |
| **limma**, **compositions**, **vegan** | Batch correction and ecological analysis |
| **pandas**, **numpy**, **matplotlib**, **seaborn** | Data processing & plotting |
| **umap**, **Rtsne** | Dimensionality reduction |
| **gridExtra**, **ggplot2** | Multi-panel figure composition |

---

##  **Example Usage**

```bash
# Step 1: Build TPM matrix
python3 01_Species_TPM matrix_for_limma.py

# Step 2: CLR transform and batch correction
Rscript 02_CLRtransformation_batch_correction.R

# Step 3: Collect outputs for downstream modules
bash 03_final_collect_files_downstream.sh

# Step 4–7: Run additional QC and performance diagnostics
python3 04_fastqscreen.py
python3 05_mag_quality.py
python3 06_quast.py
python3 07_runtime_analytics.py
```


## **Integration with GRUMB Framework**
This module represents Stage 2 of the GRUMB workflow, bridging raw metagenomic outputs and analytical modeling.
It ensures that only high-quality, batch-harmonized, and biologically coherent data progress to downstream ecological and machine-learning analyses.