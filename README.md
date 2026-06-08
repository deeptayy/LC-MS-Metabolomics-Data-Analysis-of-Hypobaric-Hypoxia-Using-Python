# Hypobaric Hypoxia Metabolomics in Rats

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://www.python.org/)
[![Metabolomics](https://img.shields.io/badge/omics-metabolomics-green)]()

This is a bioinformatics project based on my M.Sc. thesis. The project shows how I organized metadata, wrote a simple Python analysis workflow, and reported metabolomics results cautiously without sharing unpublished raw data.

The public repository does not include raw metabolomics intensity values, complete statistical result tables, or confidential datasets. A small simulated demo dataset is included only so the code can be run and reviewed.

## Scientific Question

How does hypobaric hypoxia affect the metabolomic profile of rats compared with control animals?

## Study Overview

| Component | Public summary |
| --- | --- |
| Organism | Male Sprague Dawley rats |
| Groups | Control and hypobaric hypoxia |
| Sample size | 20 total samples, 10 per group |
| Assay | Untargeted LC-MS/MS metabolomics |
| Features detected in raw processed table | 2,433 metabolites/features |
| Public data included | Sample metadata summary and simulated demo data |
| Private data excluded | Raw intensities, complete result tables, detailed metabolite-level findings |
| Analysis focus | PCA, differential feature screening, volcano plot, exploratory PLS/VIP |

## Repository Structure

```text
.
├── README.md
├── LICENSE
├── requirements.txt
├── data/
│   ├── README.md
│   └── demo/
│       └── C_vs_H_demo/
│           ├── metadata.csv
│           └── raw_data.csv
├── docs/
│   ├── example_figures/
│   │   ├── PCA_demo.png
│   │   └── Volcano_plot_demo.png
│   ├── methodology.md
│   ├── public_summary.md
│   └── reproducibility.md
├── results/
│   └── README.md
└── scripts/
    └── metabolomics_pipeline.py
```

## Analysis Steps

The private thesis workflow included:

1. Loading the metabolomics feature table and metadata.
2. Matching sample IDs with group labels.
3. Filtering metabolites/features based on missing values.
4. Running PCA to check group-level variation.
5. Performing Welch's t-tests and calculating log2 fold change.
6. Applying Benjamini-Hochberg FDR correction.
7. Creating a volcano plot.
8. Running exploratory PLS/VIP analysis.
9. Interpreting candidate metabolites cautiously.

## Results Summary

The analysis suggested broad metabolomic differences between control and hypoxia groups. At a high level, the exploratory analysis showed:

- group-level separation in PCA;
- a subset of metabolites/features passing screening criteria;
- biological themes related to amino acid metabolism, lipid metabolism, stress response, and host-microbiome metabolic activity.

Exact metabolite names, intensity values, p-values, fold changes, VIP scores, AUC values, and complete result tables are not shared publicly because the work is unpublished or restricted.

These findings are hypothesis-generating. Candidate metabolites should not be treated as confirmed biomarkers without targeted LC-MS/MS validation, authentic standards, and further biological replication.

## Example Figures

The files in `docs/example_figures/` are generated from simulated demo data. They are placeholders that show the type of visualization used in the workflow, not private thesis results.


## Skills Demonstrated

- Organizing biological metadata for analysis.
- Reading and processing CSV files with Python.
- Using pandas, NumPy, SciPy, scikit-learn, statsmodels, and matplotlib.
- Applying missing-value filtering.
- Running PCA and differential screening.
- Creating exploratory visualizations such as PCA and volcano plots.


