# Data

This folder contains only public-safe files. The original thesis metabolomics intensity table is not included.

## Public Dataset Summary

| Item | Summary |
| --- | --- |
| Comparison | Control vs hypobaric hypoxia |
| Samples analyzed privately | 20 total, 10 per group |
| Features detected in raw processed table | 2,433 metabolites/features |
| Instrument/assay | Untargeted LC-MS/MS metabolomics |
| Public raw intensity values | Not shared |
| Complete private result tables | Not shared |

## Included Files

| File | Description |
| --- | --- |
| `demo/C_vs_H_demo/metadata.csv` | Simulated sample metadata with `SampleID` and `Group` columns. |
| `demo/C_vs_H_demo/raw_data.csv` | Simulated feature-intensity table with fake feature IDs. This is not thesis data. |

## Private Files Not Included

Do not commit these files unless public sharing is approved:

- unredacted raw metabolomics intensity tables;
- vendor/raw LC-MS/MS files;
- complete processed feature tables;
- full t-test, FDR, fold-change, VIP, biomarker, pathway, or ROC result tables from the thesis analysis;
- unpublished metabolite identities if they reveal research findings.

The demo files are included only to show that the beginner analysis workflow can run.
