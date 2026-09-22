# Stabilizing IPCW Conformal Lower Predictive Bounds for Right-Censored Survival Data by Weight Clipping

This repository contains the R code and outputs accompanying the manuscript:

**Nyangweso, H. N. and Wang, H.  
“Stabilizing IPCW Conformal Lower Predictive Bounds for Right-Censored Survival Data by Weight Clipping.”**

## Manuscript outputs

| Output | Analysis | Script |
|---|---|---|
| Table 1 | Simulation 1: internal ablation | `scripts/01_run-simulation1.R` |
| Figure 1 | Simulation 1 | `scripts/figure1_simulation1.R` |
| Table 2 | Simulation 2: censoring stress | `scripts/02_run-simulation2.R` |
| Figure 2 | Simulation 2 | `scripts/figure2_simulation2.R` |
| Table 3 | Simulation 3: external benchmark | `scripts/03b_run-simulation3-external-benchmark.R` |
| Figure 3 | Simulation 3 | `scripts/figure3_simulation3.R` |
| Table 4 | GBSG real-data application | `scripts/04_run-real-data-GBSG.R` |
| Figure 4 | GBSG application | `scripts/figure4_gbsg_application.R` |
| Table 5 | Simulation data-generating mechanisms | Simulation code and manuscript |
| Table 6 | Clipping sensitivity | `scripts/06_clip-sensitivity.R`, `scripts/07_make-clipping-sensitivity-table.R` |

## Repository structure

- `R/` — core calibration, weighting, simulation, and real-data functions.
- `scripts/` — simulation drivers, figure scripts, and sensitivity analyses.
- `results/` — simulation and real-data outputs used in the manuscript.
- `figures/` — manuscript figures in PDF, TIFF, and PNG formats.

## Manuscript figures

The four manuscript figures are generated using:

```r
source("scripts/figure1_simulation1.R")
source("scripts/figure2_simulation2.R")
source("scripts/figure3_simulation3.R")
source("scripts/figure4_gbsg_application.R")
