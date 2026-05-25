---
name: run-bioinformatics-pipeline
description: Run the GRUMB bioinformatics pipeline (Module 1) for genome-resolved metagenomics. Covers quality control, assembly, binning, taxonomy, and functional annotation. Requires HPC/SLURM.
---

# Run Bioinformatics Pipeline

## Prerequisites

- Linux (Ubuntu 20.04+)
- Python >= 3.11, R >= 4.2
- SLURM job scheduler
- >= 96 GB RAM (typical HPC node)
- Tools: MEGAHIT, MetaBAT2, CheckM2, dRep, GTDB-Tk, DIAMOND, RGI, BBTools, FastQC, FastQ Screen

## Devin Secrets Needed

None — pipeline processes local FASTQ files. HPC credentials may be needed for cluster submission.

## Pipeline Steps (Module 1: Bioinformatics)

### Step 1: Quality control and host filtering
```bash
bash 01_Bioinformatics/01_actual_New_run.sh
```
Performs QC with FastQC, host read filtering with BBTools, and contamination screening with FastQ Screen.

### Step 2: Assembly and binning
```bash
bash 01_Bioinformatics/02_Assembly_Binning.sh
```
Assembles reads with MEGAHIT and bins contigs with MetaBAT2.

### Step 3: Dereplication, taxonomy, and abundance
```bash
bash 01_Bioinformatics/03_global_drep_taxonomy_abundance.sh
```
Dereplicates MAGs at 95% ANI with dRep, assigns taxonomy with GTDB-Tk, and computes species-level TPM abundance matrices.

### Step 4: ARG screening
```bash
bash 01_Bioinformatics/04_global_ARG.sh
```
Screens for antimicrobial resistance genes using DIAMOND against the CARD database.

### Step 5: Virulence factor screening
```bash
bash 01_Bioinformatics/05_global_VFDB.sh
```
Screens for virulence factors using DIAMOND against the VFDB database.

### Step 6: Collect and merge files
```bash
bash 01_Bioinformatics/collect_scripts/01_Global_file_abundance_processsing.sh
bash 01_Bioinformatics/collect_scripts/02_Merge_abundances.sh
bash 01_Bioinformatics/collect_scripts/03_Mapping_to_species.sh
bash 01_Bioinformatics/collect_scripts/04_Selecting_ARG_VF_genomes.sh
bash 01_Bioinformatics/collect_scripts/05_labelling_species.sh
```

## Pipeline Steps (Module 2: Quality Assessment & Batch Correction)

```bash
python3 02_Quality_batch_subsetting/01_Species_TPM\ matrix_for_limma.py
Rscript 02_Quality_batch_subsetting/02_CLRtransformation_batch_correction.R
bash 02_Quality_batch_subsetting/03_final_collect_files_downstream.sh
```

Additional QC scripts:
```bash
python3 02_Quality_batch_subsetting/04_fastqscreen.py
python3 02_Quality_batch_subsetting/04_mag_quality_metrics_analysis.py
python3 02_Quality_batch_subsetting/05_mag_quality.py
python3 02_Quality_batch_subsetting/06_quast.py
Rscript 02_Quality_batch_subsetting/06_Batch_Correction_PCA_UMAP_tsne.R.R
python3 02_Quality_batch_subsetting/07_runtime_analytics.py
```

## Output

| Output | Description |
|--------|-------------|
| High-quality MAGs | Dereplicated, quality-filtered genome bins |
| Species-level TPM matrices | Raw, CLR-transformed, and batch-corrected |
| ARG/VFDB annotations | Resistance gene and virulence factor profiles |
| QC reports | FastQC, CheckM2, QUAST metrics |

## Data Directory Structure

| Path | Contents |
|------|----------|
| `Data/CARD_VFDB_DATA/` | Reference databases for functional annotation |
| `Data/Corrected_Data_withAligned_Metadata/` | Batch-corrected abundance data |
| `Data/DATA_from_HPC/` | Raw HPC output files |
| `Data/Subset_Data/` | Filtered/subset data for downstream analysis |

## Notes

- Each module is fully independent and can be deployed separately
- Scripts are designed for SLURM — adapt `sbatch` headers for your cluster
- Average runtime: ~22 min per sample for Module 1
