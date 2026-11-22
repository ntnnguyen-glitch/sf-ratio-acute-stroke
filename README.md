# S/F-defined respiratory dysfunction in non-intubated ICU patients with acute stroke

R code for the study "Persistent early respiratory dysfunction defined by the SpO2/FiO2 ratio and mortality in non-intubated critically ill patients with acute stroke: a retrospective cohort study in two intensive care databases".

## Data

The analysis uses two databases hosted on PhysioNet. Access requires credentialing and signed data use agreements. No patient-level data are included in this repository.

- MIMIC-IV v3.1: https://doi.org/10.13026/kpb9-mt58
- eICU Collaborative Research Database v2.0: https://doi.org/10.13026/C2WM1R

Point the scripts to the compressed CSV files with two environment variables:

```
export MIMIC_DIR=/path/to/mimic-iv-3.1      # contains hosp/ and icu/
export EICU_DIR=/path/to/eicu-crd-2.0       # contains patient.csv.gz, lab.csv.gz, ...
```

## Pipeline

| Script | Content |
|---|---|
| `R/00_setup.R` | Paths and packages |
| `R/01_extract.R` | Extraction from the raw files into `work/extract.duckdb` |
| `R/02_cohort.R` | 24-h landmark, hourly SpO2/FiO2, exposure groups, covariates and outcomes (`work/analytic.rds`) |
| `R/03_analysis.R` | Multiple imputation, logistic and Fine-Gray models, splines, overlap weighting, E-values, sensitivity analyses and exploratory effect modification (`results/results.rds`) |
| `R/04_tables_figures.R` | Baseline tables, missing data and Figures 1-4 (`results/`) |

Run everything from the repository root:

```
Rscript run_all.R
```

Extraction reads several large files (chartevents, labevents, vitalPeriodic, respiratoryCharting) and needs about 16 GB of memory. Files in `work/` contain patient-level data and must not be shared.

## Definitions

- Cohort: adults at their first ICU stay for ischaemic stroke, intracerebral haemorrhage or subarachnoid haemorrhage, alive and in the ICU at 24 h and not invasively ventilated during the first 24 h.
- Exposure: each SpO2 reading (60-100%) divided by the most recent FiO2 within the previous 4 h, summarised as the hourly median. A charted FiO2 is used when available; otherwise oxygen flow is converted as 21% + 4% per L/min, and readings at a flow above 15 L/min without a charted FiO2 are not used; 21% when there is no oxygen record. Respiratory dysfunction is a lowest S/F <= 315 in the 0-12 h or 12-24 h window after ICU admission; it is persistent when present in both windows.
- Invasive ventilation (exclusion within 24 h and outcome after 24 h): MIMIC-IV procedure item 225792, endotracheal or tracheostomy tube as oxygen device (226732) or an invasive ventilator mode (223849); eICU-CRD respiratoryCare, treatment and nurse charting records and the APACHE intubation flag.
- Primary outcome: in-hospital death after the 24-h landmark, in both databases.
- Secondary outcomes: invasive ventilation within 7 days (death as a competing event), ICU length of stay, and 28-, 90- and 365-day mortality in MIMIC-IV.
- Missing covariates: multiple imputation by chained equations within each database (20 datasets), combined with Rubin's rules; complete-case analysis and median imputation as sensitivity analyses.

## Software

R 4.6.1 with DBI, duckdb, data.table, survival, cmprsk, splines, mice, WeightIt, cobalt, EValue, sandwich, ggplot2, ragg, scales and jsonlite.

## License

MIT
