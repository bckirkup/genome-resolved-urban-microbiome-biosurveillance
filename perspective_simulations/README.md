# Perspective Simulations

## Overview
This folder contains the scripts and input files used for a **Perspective paper submitted to _The ISME Journal_** entitled:  
**“Synthetic Ecological Blending and Predictive Foresight: Integrating the R²DI and GRUMB Frameworks for Global Microbial Biosurveillance.”**

The analyses integrate genome-resolved metagenomics, ecological diversity metrics, and synthetic blending simulations to evaluate microbial compositional resilience and contamination risk across built and urban environments.

---

## Contents
- **Ecology_Country.R** — Computes alpha and beta diversity metrics (Shannon, Simpson, Richness) and generates diversity plots.  
- **simulation_blending_isme_perspectives.py** — Runs synthetic donor–recipient community blending simulations, calculates entropy trajectories, and generates Minimal Detectable Contamination (MDC) and dominance matrices.  
- **Input_Files/** — TPM and metadata files required to reproduce analyses (available in the **Data** folder of the GRUMB repository).  
- **Output_Files/** — Derived CSV tables and publication-quality figures.

---

## Notes
All genome-resolved metagenomic datasets analyzed are publicly available through the **NCBI Sequence Read Archive (SRA)**, with project identifiers listed in the **GRUMB repository**.  
This code and processed outputs correspond to **GRUMB v2.0** and are permanently archived on Zenodo:  
➡️ [10.5281/zenodo.15505402](https://doi.org/10.5281/zenodo.15505402)