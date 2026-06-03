# Methodology

## Experimental Context

Male Sprague Dawley rats were assigned to control and hypobaric hypoxia conditions. Gastrointestinal or fecal material was collected and analyzed using untargeted LC-MS/MS metabolomics. The public repository does not include the raw thesis intensity table or complete result files.

## Public Demo Data

The repository includes a small simulated dataset so the Python workflow can be run by recruiters, research labs, or graduate supervisors. The demo data use fake feature names such as `Feature_0001` and are not thesis results.

## Computational Workflow

### 1. Data Loading

The feature table is read with metabolites/features as rows and sample IDs as columns. Metadata are read as sample-level annotations. Sample IDs are stripped of whitespace before joining the intensity matrix to the phenotype table.

### 2. Missingness Filter

Features are retained only if they are observed in at least 50% of samples in both groups. This removes sparse features that may produce unstable fold changes or unreliable t-tests.

### 3. PCA

Principal component analysis is used as an unsupervised check for broad sample structure, possible outliers, and group-level variation.

### 4. Differential Feature Screening

Welch's t-test is applied feature-by-feature. P-values are adjusted with Benjamini-Hochberg FDR correction. Log2 fold change is calculated as:

```text
Hypoxia / Control
```

This screening is exploratory and should not be treated as final biomarker discovery.

### 5. Exploratory PLS/VIP

PLS regression is fit with binary group labels to approximate a beginner-level PLS-DA style workflow. VIP scores are calculated to rank features contributing to group separation.

Because the sample size is small, PLS/VIP results should be described cautiously and validated before making strong biological claims.

### 6. Candidate Feature Selection

Demo candidate features are selected using:

- nominal p-value < 0.05;
- absolute log2 fold change >= 1;
- VIP score > 1.

For the private thesis analysis, exact candidate names, p-values, fold changes, VIP scores, and result tables should remain private unless public sharing is approved.

### 7. ROC Analysis

The top-ranked demo feature is evaluated with an ROC curve. In publication-grade work, ROC results should be validated with held-out or independent samples.

## Biological Interpretation

The thesis context highlights possible hypoxia-linked changes in microbial metabolism, amino acid metabolism, lipid metabolism, stress response, gut barrier biology, and host-microbiome signaling.

The public repository reports these findings only at a high level. Untargeted metabolomics findings are hypothesis-generating and require targeted validation.
