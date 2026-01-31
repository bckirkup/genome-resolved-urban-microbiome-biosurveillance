#  **GRUMB : Genome-Resolved Urban Microbiome Biosurveillance**

---

 
### *A fully integrated framework for genome-resolved metagenomics, ecological modeling, and predictive risk assessment*  

**GRUMB (Genome-Resolved Urban Microbiome Biosurveillance)** is an open-source, SLURM-compatible pipeline that unites **assembly-based metagenomics, ecological modeling, and simulation-driven machine learning** for infrastructure-scale biosurveillance.  
It reconstructs **metagenome-assembled genomes (MAGs)** from raw shotgun reads, performs **taxonomic and functional annotation**, applies **batch-aware normalization**, and models **environmental stability, cross-contamination, and risk** through interpretable, uncertainty-aware ML diagnostics.  

Built for **urban and human-built environments** such as ambulances, hospital environments, hospital sewage, and  public transport systems, GRUMB transforms descriptive microbiome profiling into quantitative, predictive biosurveillance.  

---

##  What This Framework Offers  

| Capability | Description |
|-------------|-------------|
| **Genome-resolved profiling** | From raw FASTQ or SRA accessions → high-quality MAGs via MEGAHIT assembly, MetaBAT2 binning, CheckM2 QC, and GTDB-Tk taxonomy. |
| **Functional annotation** | ARG  and virulence gene screening via DIAMOND and RGI integration. |
| **Ecological modeling** | Alpha/beta diversity metrics (Shannon, Simpson, PERMANOVA, ANOSIM), indicator species identification, and WHO-pathogen network mapping. |
| **Machine learning diagnostics** | Environment classification (Random Forest / Decision Tree) with nested CV, feature importance ranking, and entropy-based uncertainty. |
| **Simulation-based inference** | Synthetic donor–recipient mixing to quantify Minimal Detectable Contamination (MDC), entropy tipping points, and source–sink dominance. |
| **Risk quantification** | Composite contamination risk index combining override frequency, entropy, richness, and variance to rank environmental vulnerability. |

Designed for **reproducibility, explainability, and cross-environment comparability**, GRUMB  operationalizes biosurveillance as a quantitative, diagnostic process.  

---

##  Pipeline Overview 

```text
             ┌──────────────────────────────────────────────────────────────┐
             │                 GRUMB : Analytical Workflow                  │
             └──────────────────────────────────────────────────────────────┘

         ┌──────────────┐
         │  Raw Reads   │     Shotgun metagenomes (FASTQ/SRA)
         └──────┬───────┘
                │
                ▼
┌────────────────────────────────────────────────────────────────────────────┐
│   MODULE 1: BIOINFORMATICS                                                 │
│  - Quality control & host filtering (BBTools, FastQ Screen)                │
│  - Assembly (MEGAHIT) + Binning (MetaBAT2)                                 │
│  - MAG quality (CheckM2) + Dereplication (dRep, 95% ANI)                   │
│  - Taxonomy (GTDB-Tk) + Functional annotation (CARD, VFDB)                 │
│  - Output: High-quality dereplicated MAGs + species-level TPM matrices     │ 
│  - Runtime: ~22 min average per sample                                     │                                                            
└────────────────────────────────────────────────────────────────────────────┘
                │
                ▼
┌────────────────────────────────────────────────────────────────────────────┐
│   MODULE 2: QUALITY ASSESSMENT & BATCH CORRECTION                          │
│  - Sequencing coverage (Nonpareil) + MAG/assembly metrics (QUAST)          │
│  - CLR transformation + limma batch correction                             │
│  - PCA, UMAP, t-SNE, PERMANOVA diagnostics                                 │
│  - Output: Batch-corrected TPM matrices, harmonized data                   │
└────────────────────────────────────────────────────────────────────────────┘
                │
                ▼
┌────────────────────────────────────────────────────────────────────────────┐
│   MODULE 3: ECOLOGY & DIVERSITY                                            │
│  - Alpha & Beta diversity (Shannon, Simpson, PERMANOVA, ANOSIM)            │
│  - Indicator species & NMDS ordination (indicspecies)                      │
│  - WHO-priority pathogen network & co-occurrence modeling                  │
│  - Output: Species prevalence, network graphs, ecological fingerprints     │
└────────────────────────────────────────────────────────────────────────────┘
                │
                ▼
┌────────────────────────────────────────────────────────────────────────────┐
│  MODULE 4: MACHINE LEARNING & ECOLOGICAL SIMULATIONS                       │
│  - Environment classification (Random Forest, DT, SVM, GB, LR)             │
│  - Feature importance (Gini, permutation, consensus ranking)               │
│  - Synthetic donor–recipient mixing & entropy-based tipping points         │
│  - Risk scoring: override freq + entropy + richness + variance             │
│  - Output: Contamination thresholds, entropy trajectories, risk indices    │
└────────────────────────────────────────────────────────────────────────────┘
                │
                ▼
┌────────────────────────────────────────────────────────────────────────────┐
│ FINAL OUTPUTS                                                              │
│  - MAGs + Annotations (CARD/VFDB)                                          │
│  - TPM matrices (raw, CLR, corrected)                                      │
│  - Ecological diversity & network plots                                    │
│  - ML benchmarks + simulation results                                      │
│  - Composite ecological risk profiles per environment                      │
└────────────────────────────────────────────────────────────────────────────┘

```
#### **Each module is fully independent and can be deployed separately or executed end-to-end on HPC clusters. For more detailed descriptions of the modules, check the individual README in the sub-folders**

## Dependencies  

| Category | Software / Library |
|-----------|--------------------|
| **Core** | Python ≥ 3.11, R ≥ 4.2, SLURM scheduler |
| **Assembly & Binning** | MEGAHIT, MetaBAT2, CheckM2, dRep |
| **Taxonomy & Annotation** | GTDB-Tk, DIAMOND, RGI |
| **QC & Metrics** | FastQC, BBTools, FastQ Screen, QUAST, Nonpareil |
| **Normalization & Ecology** | limma, vegan, indicspecies, compositions, umap, Rtsne |
| **Machine Learning & Simulation** | scikit-learn, pandas, numpy, seaborn, matplotlib, scipy, scikit-posthocs |

---

##  Example Execution  

```bash
# Bioinformatics Module
bash 01_Bioinformatics/01_actual_New_run.sh
bash 01_Bioinformatics/02_Assembly_Binning.sh
bash 01_Bioinformatics/03_global_drep_taxonomy_abundance.sh

# Quality Assessment + Batch Correction
python3 02_quality_assessements_batch_correction/01_Species_TPM\ matrix_for_limma.py
Rscript 02_quality_assessements_batch_correction/02_CLRtransformation_batch_correction.R
bash 02_quality_assessements_batch_correction/03_final_collect_files_downstream.sh

# Ecology and Machine Learning
Rscript 03_Ecology/01_Ecology_Analysis.R
python3 04_Machine_Learning_and_Ecological_Blending_Simulation/03_Contaminations_Simulations_Ecological_Scoring.py
```

## System Requirements
```
Linux (Ubuntu 20.04 +)
Python 3.11 +, R 4.2 +
≥ 96 GB RAM (typical HPC node)
SLURM or equivalent job scheduler
```

---

## License

MIT License — free to use, adapt, and cite with attribution.

---

## Version and Citation Notice  

This repository hosts **GRUMB v2.0**, the **updated and extended version** of the *Genome-Resolved Urban Microbiome Biosurveillance (GRUMB)* framework originally published in **Bioinformatics**:

> **Aminu S., Ascandari A., Benhida R., Daoud R. (2025)**  
> *GRUMB: A Genome-Resolved Metagenomic Framework for Monitoring Urban Microbiomes and Diagnosing Pathogen Risk.*  
> *Bioinformatics, Oxford University Press.*  
> [https://doi.org/10.1093/bioinformatics/btaf548](https://doi.org/10.1093/bioinformatics/btaf548)  

---

### Framework Update

**GRUMB v2.0** builds upon the original GRUMB architecture with substantial extensions:  

- Integration of **simulation-driven ecological modeling**, including synthetic microbial blending and entropy-based stability analysis.  
- Expanded **machine learning diagnostics**, featuring nested cross-validation, feature stability analysis, and predictive risk scoring.  
- Addition of **batch correction and cross-study harmonization** modules (CLR + limma).  
- Enhanced **computational scalability** and runtime tracking across  metagenomes from multiple projects.  
 
These updates were published at **Microbiome**:  


---

### Citation Guidelines  

If you use this framework, please **cite both works** as follows:  

> **Primary Framework (Bioinformatics):**  
> Aminu S., Ascandari A., Benhida R., Daoud R. (2025).  
> *GRUMB: A Genome-Resolved Metagenomic Framework for Monitoring Urban Microbiomes and Diagnosing Pathogen Risk.*  
> *Bioinformatics.* DOI: [10.1093/bioinformatics/btaf548](https://doi.org/10.1093/bioinformatics/btaf548)

> **Companion Study (Microbiome):**  
> Aminu, S., Ascandari, A., Mokhtar, M.M. et al.  
> *Genome-resolved surveillance and predictive ecological risk modeling of urban microbiomes.*  
> *Microbiome* 14, 45 (2026). https://doi.org/10.1186/s40168-025-02315-3
>  


A **Zenodo DOI** is also available for version tracking and reproducibility:  
[https://doi.org/10.5281/zenodo.15505402](https://doi.org/10.5281/zenodo.15505402)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15505402.svg)](https://doi.org/10.5281/zenodo.15505402)


###  Summary  

> GRUMB v2.0 bridges genome-resolved metagenomics with ecological modeling and simulation-based machine learning to advance predictive biosurveillance under the **One Health** framework.  
> It represents a unified analytical ecosystem, extending from genome recovery to ecological resilience modeling designed for reproducibility, transparency, and cross-environment comparability.

---

### Contact
For questions, feedback, or collaboration regarding this framework, please reach out:

Suleiman Aminu (Lead Developer),
PhD Researcher,
Department of Chemical and Biochemical Sciences,
University Mohammed VI Polytechnic (UM6P), Morocco,
suleiman.aminu@um6p.ma; saminu83@gmail.com

Abdulaziz Ascandari,
PhD Researcher,
Department of Chemical and Biochemical Sciences,
University Mohammed VI Polytechnic (UM6P), Morocco,
abdulaziz.ascandari@um6p.ma; ryandari87@gmail.com

Prof. Rachid Daoud,
Group Leader & Supervisor,
Department of Chemical and Biochemical Sciences, 
University Mohammed VI Polytechnic (UM6P), Morocco,
rachid.daoud@um6p.ma

### Acknowledgements
This work was supported by internal funding from Mohammed VI Polytechnic University (UM6P) and computational resources from the African Supercomputing Center (https://toubkal.um6p.ma).
We thank all members of the Chemical and Biochemical Sciences Department for their contributions to the GRUMB project.
