---
name: run-ml-analysis
description: Run the GRUMB ecology analysis and machine learning pipeline (Modules 3-4). Covers diversity analysis, indicator species, environment classification, contamination simulations, and risk scoring.
---

# Run ML & Ecology Analysis

## Prerequisites

- Python >= 3.11 with scikit-learn, pandas, numpy, scipy, seaborn, matplotlib, scikit-posthocs
- R >= 4.2 with vegan, indicspecies, compositions, umap, Rtsne
- Batch-corrected TPM matrices from Module 2 (run bioinformatics pipeline first)

## Devin Secrets Needed

None — analysis runs on local data files.

## Install Python Dependencies
```bash
pip install scikit-learn pandas numpy scipy seaborn matplotlib scikit-posthocs
```

## Install R Dependencies
```r
install.packages(c("vegan", "indicspecies", "compositions", "umap", "Rtsne", "ggplot2", "tidyverse"))
```

## Module 3: Ecology & Diversity

### Alpha/beta diversity analysis
```bash
Rscript 03_Ecology/01_R_pipeline_Ecology_Including.R
```
Computes Shannon, Simpson diversity indices, PERMANOVA, ANOSIM, and NMDS ordination.

### Indicator species analysis
```bash
Rscript 03_Ecology/02_Indicator_Specie_Analysis.R
```
Identifies environment-specific indicator taxa using the indicspecies package.

### WHO priority pathogen prevalence
```bash
Rscript 03_Ecology/03_WHO_Priority_Prevalence.R
```
Maps WHO-priority pathogens across environments with network visualization.

### Species community analysis
```bash
Rscript 03_Ecology/08_species_community_analysis.R
```

### Pathogen prevalence and diversity (Python)
```bash
python3 03_Ecology/09_pathogen_prevalence_diversity.py
```

## Module 4: Machine Learning & Ecological Simulations

### Model selection with nested cross-validation
```bash
python3 04_Machine_Learning/01_Model_selection_Holdout_Nested.py
```
Trains Random Forest, Decision Tree, SVM, Gradient Boosting, and Logistic Regression classifiers with nested CV for environment classification.

### Feature importance and stability
```bash
python3 04_Machine_Learning/02_Feature_importance_stability.py
```
Ranks features by Gini importance, permutation importance, and consensus ranking.

### Contamination simulations and ecological scoring
```bash
python3 04_Machine_Learning/03_Contaminations_Simulations_Ecological_Scoring.py
```
Runs synthetic donor-recipient mixing simulations to quantify:
- Minimal Detectable Contamination (MDC)
- Entropy-based tipping points
- Source-sink dominance thresholds
- Composite contamination risk index

### Additional ML scripts
```bash
python3 04_Machine_Learning/10_rf_model_training.py
python3 04_Machine_Learning/11_model_comparison.py
Rscript 04_Machine_Learning/12_ML_Comparison_Script.R
python3 04_Machine_Learning/13_ML_Env_Biosurveillance.py
```

## Perspective Simulations
```bash
# Additional simulation scripts in perspective_simulations/
ls perspective_simulations/
```

## Key Outputs

| Output | Description |
|--------|-------------|
| Diversity metrics | Shannon, Simpson, richness per environment |
| Ordination plots | NMDS, PCA, UMAP, t-SNE visualizations |
| ML benchmarks | Accuracy, F1, AUC for each classifier |
| Feature rankings | Consensus feature importance across methods |
| Risk indices | Composite contamination risk per environment |
| MDC thresholds | Minimum detectable contamination fractions |
| Entropy trajectories | Stability curves for microbial communities |

## Notes

- Modules 3 and 4 require output from Modules 1-2 (batch-corrected TPM matrices)
- ML models use nested cross-validation to avoid overfitting
- Contamination simulations are computationally intensive — expect longer runtimes on large datasets
