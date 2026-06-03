# Reproducibility Guide

This public repository is reproducible with simulated demo data. The private thesis analysis requires restricted raw intensity values that are not included here.

## Environment

Install dependencies from the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

On Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

## Run The Public Demo

```bash
python scripts/metabolomics_pipeline.py --data-root data/demo --results-root results_demo
```

The demo data are simulated and use fake feature IDs. They are included only to show that the workflow runs.

## Expected Inputs

Each comparison folder should contain:

- `raw_data.csv`: features as rows and sample IDs as columns;
- `metadata.csv`: sample IDs and group labels.

The public demo is located at:

```text
data/demo/C_vs_H_demo/
```

For the private thesis analysis, restore the real intensity table locally in a private folder such as `data/private/`. Do not commit private raw data or complete private outputs.

## Expected Demo Outputs

Running the demo creates a folder such as:

```text
results_demo/C_vs_H_demo/
```

Typical demo outputs include:

- filtered data matrix;
- PCA scores and figure;
- t-test/FDR/log2FC table;
- volcano plot;
- exploratory PLS score table and figure;
- VIP score table and figure;
- demo biomarker table;
- demo pathway input table;
- ROC figure if a demo feature passes screening criteria;
- run summary text file.

These are simulated outputs, not thesis results.

## Validation Checklist

Before reporting private results, confirm:

- sample IDs match between raw data and metadata;
- group labels are correct;
- fold-change direction is `Hypoxia / Control`;
- missingness before and after filtering is reasonable;
- PCA is checked for outlier samples;
- nominal p-values and FDR-adjusted p-values are interpreted separately;
- PLS/VIP and ROC results are described as exploratory unless validated independently;
- possible contaminants, exogenous compounds, and uncertain annotations are flagged.
