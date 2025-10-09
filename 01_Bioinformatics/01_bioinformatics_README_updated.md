# **01_Bioinformatics Module**
### *Genome-Resolved Urban Microbiome Biosurveillance: Assembly, Binning, and Annotation Pipeline*

This module implements the **bioinformatics backbone** of the GRUMB framework, providing a reproducible, SLURM-compatible pipeline for genome-resolved metagenomic analysis.  
It supports large-scale biosurveillance studies spanning  **ambulances**, **hospital environments**, **hospital sewage**, and **public transport** (key interfaces in the One Health continuum).

---

## **Overview**

The Bioinformatics Module performs:

- **Read-level preprocessing, decontamination, and QC** using BBTools and FastQ Screen  
- **Metagenome assembly and binning** for MAG recovery  
- **Genome quality control and dereplication** (CheckM2 + dRep)  
- **Taxonomic and functional annotation** of AMR and virulence genes  
- **Automated output aggregation** for ecological and machine learning analysis  

The output is a curated collection of high-quality metagenome-assembled genomes (MAGs) and species-level abundance tables for downstream ecological modeling.

---

##  **Directory Structure**
```
01_Bioinformatics/
├── 01_actual_New_run.sh # Preprocessing, decontamination, and quality control
├── 02_Assembly_Binning.sh # Assembly → binning and MAG recovery
├── 03_global_drep_taxonomy_abundance.sh# Dereplication, taxonomy, and abundance quantification
├── 04_global_ARG.sh # AMR gene profiling using RGI (CARD)
├── 05_global_VFDB.sh # Virulence factor profiling using DIAMOND (VFDB)
├── delete.sh # Removes intermediate and large temporary files
├── Nonpareil_analysis.R # Sequencing coverage estimation and visualization
├── simplifyFastaHeaders.pl # Header normalization for downstream compatibility
├── /collect_scripts/ # Aggregation and downstream prep scripts
│ ├── 01_Global_file_abundance_processsing.sh
│ ├── 02_Merge_abundances.sh
│ ├── 03_Mapping_to_species.sh
│ ├── 04_Selecting_ARG_VF_genomes.sh
│ ├── 05_labelling_species.sh
│ └── README.md
└── README.md
```
---

##  **Workflow Summary**

### **1️ Read Preprocessing, Quality Control, and Contaminant Removal**  
**Script:** `01_actual_New_run.sh`

- Downloads reads via `pysradb` + `fasterq-dump`
- Performs trimming, adapter removal, and Phix filtering using `bbduk.sh`  
  (`qtrim=rl`, `trimq=20`, `minlen=25`)
- Repairs paired-end reads with `repair.sh`  
- Removes duplicates with `clumpify.sh`
- Merges overlapping reads using `bbmerge.sh`
- Assesses read quality using `FastQC`
- Removes host/contaminant sequences using **FastQ Screen v0.15.3**, against:
  - *Homo sapiens* (GRCh38)
  - *Mus musculus* (GRCm39)
  - *Plasmodium spp.*
  - Vector database (UniVec)
- Evaluates sequencing redundancy and coverage using **Nonpareil v3.5.5**

*Output:**
- Cleaned paired FASTQ files (`*_R1/2.decontam.paired.fq.gz`)
- Singleton files (`*_singletons.decontam.fq.gz`)
- Host-screening reports and Nonpareil coverage plots

---

### **2️ Metagenome Assembly and Genome Binning**  
**Script:** `02_Assembly_Binning.sh`

- **Assembly:** MEGAHIT (`--presets meta-large`)
- Filters contigs < 1.5 kbp
- **Binning:** MetaBAT2 (coverage-based binning with BBMap)
- **Genome quality:** CheckM2 (≥ 80 % completeness, ≤ 10 % contamination)
- **Dereplication:** dRep (95 % ANI threshold per environment and globally)

*Output:**
- Medium- and high-quality MAGs  
- CheckM2 quality summary tables  
- Environment-specific bin directories

---

### **3 Dereplication, Taxonomy, and Abundance Quantification**  
**Script:** `03_global_drep_taxonomy_abundance.sh`

- Performs global dereplication of MAGs across all environments
- Taxonomic classification with **GTDB-Tk v2.4.1 (GTDB R226)**
- Species-level abundance quantification using **CoverM (TPM, RPKM, covered_fraction)**
- Generates global abundance matrices for downstream ecological modeling

*Output:**
- Dereplicated MAG set (95 % ANI clusters)
- GTDB-based taxonomy tables
- Species-by-sample TPM matrix

---

### **4️ AMR Gene Profiling (CARD)**  
**Script:** `04_global_ARG.sh`

- Uses **RGI v6.0.4 (CARD 2025)** for AMR annotation
- Input type: `--input_type contig`
- Retains hits classified as *strict* or *perfect* with ≥ 80 % identity and coverage
- Extracts matched contig sequences and builds Bowtie2 indices for abundance quantification

*Output:**
- Annotated ARG tables per MAG  
- Mapped read coverage for ARG-bearing contigs

---

### **5 Virulence Factor Profiling (VFDB)**  
**Script:** `05_global_VFDB.sh`

- Executes **DIAMOND v2.1.10** `blastx` search against VFDB core dataset (May 2025)
- Filters hits with ≥ 80 % identity, ≥ 70 % coverage, e-value ≤ 1e−5
- Extracts corresponding nucleotide contigs for abundance quantification

**Output:**
- Annotated VF tables per MAG  
- VFDB-mapped abundance matrices  

---

### **6 Output Cleanup**  
**Script:** `delete.sh`

- Removes large intermediate files to conserve storage
- Keeps essential outputs for reproducibility 

---

### **7 Sequencing Coverage Analysis**  
**Script:** `Nonpareil_analysis.R`

- Estimates redundancy and coverage across all environments  
- Generates coverage curves for ambulances, hospital environments, sewage, and public transport  
- Confirms saturation (>95% coverage) at ~1 Gbp sequencing depth

*Output:**  
- Nonpareil curves (`nonpareil_summary_all_projects.csv`)  
- Publication-ready coverage figures  

---

## **Dependencies**

| Tool | Version | Function |
|------|----------|-----------|
| **BBTools suite** | ≥ 37.0 | Read QC and trimming |
| **FastQC** | ≥ 0.11 | Quality assessment |
| **FastQ Screen** | ≥ 0.15.3 | Host decontamination |
| **MEGAHIT** | ≥ 1.2.9 | Assembly |
| **MetaBAT2** | ≥ 2.18 | Genome binning |
| **CheckM2** | ≥ 1.0.1 | Genome quality control |
| **dRep** | ≥ 3.4.2 | Dereplication |
| **GTDB-Tk** | ≥ 2.4.1 | Taxonomic classification |
| **RGI** | ≥ 6.0.4 | ARG annotation |
| **DIAMOND** | ≥ 2.1.10 | VF annotation |
| **CoverM** | ≥ 0.6.1 | Abundance quantification |
| **Nonpareil** | ≥ 3.0 | Coverage estimation |
| **Perl 5.x** | — | FASTA header utilities |

---

##  **Example Usage**

```bash
bash 01_actual_New_run.sh PRJNA123456
bash 02_Assembly_Binning.sh
bash 03_global_drep_taxonomy_abundance.sh
bash 04_global_ARG.sh
bash 05_global_VFDB.sh
Rscript Nonpareil_analysis.R
bash collect_scripts/06_final_collect_files_downstream.sh
```



## **Runtime Benchmark**

| **Workflow Stage**         | **Avg Runtime (min)**| **Std Dev (±)** | **Notes**                              |
|----------------------------|----------------------|----------------:|----------------------------------------|
| Preprocessing & QC         | 8.6                  | 22.3            | Trimming, filtering, host removal      |
| Assembly & Binning         | 22.4                 | 45.5            | Scales with complexity                 |
| Dereplication & QC         | 18.7                 | 38.9            | dRep + CheckM2                         |
| Annotation (CARD/VFDB)     | 9.8                  | 17.3            | Functional annotation                  |
| Output Aggregation         | 5.2                  | 10.8            | Merging + cleanup                      |

- **Average total runtime:** ~23 minutes per sample (parallelized across 56 threads, 96 GB RAM)  
- **Demonstrated reproducibility and scalability across 769 metagenomes**


##  **Notes**
Fully SLURM-compatible (parallel job scheduling)
Intermediate files automatically cleaned
Environment metadata maintained for all samples
Outputs structured for ecological, ML, and simulation analysis






