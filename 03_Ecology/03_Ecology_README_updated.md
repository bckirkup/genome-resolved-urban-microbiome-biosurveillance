---
output:
  html_document: default
  pdf_document: default
---
# *03_Ecology Module*

This module conducts **comprehensive microbial ecology analyses** integrating diversity, indicator species, and WHO-priority pathogen ecology across urban infrastructure environments.  
It includes **three major R scripts** that together profile diversity patterns, ecological markers, and priority pathogen dynamics using genome-resolved metagenomic data.

---

## **Scripts:**
```
03_Ecology/
├── 01_Diversity.R
├── 02_Indicator_Specie_analysis.R
├── 03_WHO_Priority_Prevalence.R

```

##  Scripts and Functions

## `01_Diversity.R` — Diversity, Prevalence & Environmental Ecology

## **Overview**
Performs quantitative ecological profiling of microbial communities using **genome-resolved TPM abundance matrices**.  
It computes and visualizes alpha/beta diversity, prevalence structures, and networked species sharing between environments.

## **Key Analyses**
- **Alpha Diversity (TPM-based):**
  - Shannon, Simpson, and species richness indices
  - Kruskal–Wallis and Dunn’s post-hoc comparisons across environments
  - Boxplots for Shannon, Simpson, and richness thresholds (`TPM ≥ 0.1, 1, 10`)
- **Beta Diversity (CLR-transformed):**
  - NMDS ordination (Euclidean)
  - ANOSIM, PERMANOVA, and multivariate dispersion analyses
  - Pairwise PERMANOVA with significance codes
- **Species Prevalence Ecology:**
  - Defines **Core (≥30%)**, **Secondary (≥10%)**, and **Peripheral (<10%)** prevalence classes
  - Environment-specific and global prevalence heatmaps
  - UpSet plots and overlap networks for shared taxa among environments

## **Outputs**
- `AlphaDiversity_TPM_allmetrics.csv`, `AlphaDiversity_statistics.csv`
- `PERMANOVA_CLR_Euclidean.csv`, `PERMANOVA_Pairwise_CLR_Euclidean.csv`
- `Species_Prevalence_ByEnv_all_species.csv`
- High-resolution figures:
  - `Alpha_*.tiff`, `NMDS_CLR_Euclidean_AllGroups.tiff`
  - `Prevalence_Global_relaxed.tiff`, `Heatmap_EnvSpecific_Core_Species.tiff`
  - `Network_Env_Species_Sharing.tiff`

## **Scientific Context**
Captures ecological diversity gradients, compositional heterogeneity, and shared species networks between **ambulance**, **hospital**, and **public transport** microbiomes.

---

## `02_Indicator_Specie_analysis.R` — Indicator Species Discovery

## **Overview**
Implements **Indicator Value (IndVal.g)** analysis to identify statistically significant microbial taxa associated with each environment.  
Ensures robust **Sample ID normalization** and diagnostic reporting to prevent metadata mismatches.

## **Key Analyses**
- Automated **sample name harmonization** with Unicode normalization and edit-distance diagnostics.
- **Presence/absence filtering** (`TPM > 0`) with configurable minimum presence (`min_present`, `min_per_grp`).
- **Indicator species inference** using `indicspecies::multipatt`:
  - Environment-specific and cross-environment indicators.
  - Benjamini–Hochberg correction for multiple testing.
- Visualizes top indicator species (up to 60) colored by environment.

## **Outputs**
- `indicator_DIAGNOSTICS.csv` — filtering and group summary  
- `indicator_IDs_META_normalization.csv`, `indicator_IDs_TPM_normalization.csv`  
- `indicator_SIGNIFICANT_BH.csv`, `indicator_ALL_raw.csv`  
- `Top_Indicators.tiff` — top taxa by indicator strength  

## **Scientific Context**
Defines **sentinel microbial taxa** that act as ecological markers for each environment, bridging diversity metrics with habitat specificity.

---

## `03_WHO_Priority_Prevalence.R` — WHO Priority Pathogen Ecology & Risk Networks

## **Overview**
Characterizes **WHO-priority pathogens** (Critical, High, and Other) across environments, visualizing their abundance, prevalence, and co-occurrence networks.  
Integrates WHO lists and genus-level matches for extended pathogen coverage.

## **Key Analyses**
- **WHO Pathogen Detection:**
  - Identifies pathogens from WHO priority lists and their genera.
  - Computes log-transformed mean TPM per environment.
  - Generates categorized heatmaps by WHO priority level.
- **Prevalence Classification:**
  - Core (≥80%), Secondary (50–80%), and Peripheral (<50%) species.
  - Environment-wise pathogen category counts.
- **Network Analyses:**
  - **Co-occurrence network** (Spearman |ρ| > 0.5) for all and WHO-only species.
  - WHO-pathogen **subnetworks** highlighting co-dominant taxa.
  - Centrality metrics (Degree, Betweenness, Eigenvector, Clustering Coefficient).
  - Radar plots for comparative risk profiles.
  - Community detection using Louvain modularity.

## **Outputs**
- `WHO_mean_logTPM_by_environment.csv`
- `WHO_species_prevalence_classification.csv`
- `WHO_Priority_Network_Statistics.csv`
- Publication-quality figures:
  - `WHO_heatmap_meanTPM_clean.tiff`
  - `WHO_Priority_Cooccurrence_Network_upgraded.tiff`
  - `WHO_Risk_Radar.tiff`, `WHO_Hubs_Barplot.tiff`

## **Scientific Context**
Quantifies **pathogen resilience, network centrality, and community risk**, supporting biosurveillance of high-risk microbial ecosystems in healthcare and transit environments.

---

##  **Output Summary**

| Script | Output Type | Key Files |
|--------|--------------|-----------|
| `01_Diversity.R` | Diversity, prevalence & network ecology | Alpha diversity tables, PERMANOVA, species heatmaps |
| `02_Indicator_Specie_analysis.R` | Indicator taxa per environment | Indicator value statistics, diagnostic reports |
| `03_WHO_Priority_Prevalence.R` | WHO-pathogen ecology & networks | Heatmaps, prevalence tables, co-occurrence & radar plots |

---

##  **Requirements**

- **R (≥ 4.3)**  
  Packages:  
  `vegan`, `ggplot2`, `dplyr`, `tidyr`, `tibble`, `ComplexHeatmap`, `circlize`, `indicspecies`,  
  `igraph`, `ggraph`, `FSA`, `UpSetR`, `patchwork`, `RColorBrewer`

---

## **Example Usage**

```bash
# Run ecological diversity and prevalence analyses
Rscript 01_Diversity.R

# Identify environment-specific indicator species
Rscript 02_Indicator_Specie_analysis.R

# Analyze WHO pathogen ecology and co-occurrence networks
Rscript 03_WHO_Priority_Prevalence.R

```
