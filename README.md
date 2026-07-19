# Finite-Sample Stabilization of IPCW/AIPCW Conformal LPBs

Reproduction code for *Clipped Self-Normalized IPCW and AIPCW Conformal Lower
Predictive Bounds for Right-Censored Survival Data* (Nyangweso and Wang).
Running the scripts reproduces every result and figure in the paper.

## What each script produces

| Script | Reproduces |
|--------|------------|
| `scripts/00_smoke-test.R` | quick 3-replicate end-to-end check |
| `scripts/01_run-simulation1.R` | Table 2; Figures `simulation1_coverage`, `simulation1_stability` |
| `scripts/02_run-simulation2.R` | Table 3; Figures `simulation2_coverage`, `simulation2_maxweight` |
| `scripts/03_run-simulation3.R` | Table 4, internal ablation rows |
| `scripts/03b_run-simulation3-external-benchmark.R` | Table 4, external rows (DR-COSARC, KM, Oracle); Figure `simulation3_coverage` |
| `scripts/04_run-real-data-GBSG.R` | Table 5; Figure `gbsg_panels` |
| `scripts/06_clip-sensitivity.R` | clipping-level sensitivity (Appendix D); Figure `clip_sensitivity` |
| `scripts/05_make-figures.R` | builds all figures from `results/` |
| `scripts/run_all.R` | runs everything above in order |

## Quick start

From the **repository root** (the folder containing `R/` and `scripts/`):

```r
source("R/source-code.R")
smoke_test()                       # fast check
```

Reproduce everything (full run, R = 100 per cell; long):

```r
source("scripts/run_all.R")
```

Or run pieces individually, e.g. `source("scripts/01_run-simulation1.R")`.
All outputs are written to deterministic paths under `results/`, and
`05_make-figures.R` reads those and writes PNGs to `figures/`.

### Running from RStudio
Set the working directory to the repo root first
(Session -> Set Working Directory -> Choose Directory), then `source(...)`.
`file.exists("R/source-code.R")` should return `TRUE` before you run anything.

## Dependencies

Core simulations: `survival`. Figures: `ggplot2`. GBSG application additionally:
`TH.data`, `dplyr`, `readr`, `tibble`, `tidyr`. External benchmark (`03b`)
additionally: `tidyverse`, `R6`. Pinned versions are in `sessionInfo.txt`.

## Simulation 3 external benchmark (third-party code required)

`03b` calls the conformal-survival code of Sesia & Svetnik (2024), which is
**not** distributed here because their release carried no license. To run `03b`,
obtain their code and place `utils_survival.R`, `utils_censoring.R`,
`utils_conformal.R`, and `utils_decensoring.R` in
`external_methods_sesia/code/conf_surv/` (see that folder's `NOTICE.md`).
Every other script runs without it.

## Reproducibility

Replicates are seeded as `seed = 1000 * r`, so each scenario reproduces exactly
across machines and across the clipping-level sweep; the `lambda0` pilot uses a
separate, restored RNG stream. Target coverage is 0.90, the candidate grid is
`{0.01, ..., 0.49}`, and the clipping percentile defaults to the adaptive 90th
percentile (configurable via `options(herg_clip_q = ...)`; `06` uses this).

## Repository layout

```
R/                          calibration library and data-generating code
scripts/                    numbered entry points + run_all.R
external_methods_sesia/     placeholder for Sesia & Svetnik code (see NOTICE.md)
tools/                      one-off recovery utilities (NOT part of the clean rerun)
results/  figures/          generated outputs (created on first run)
```

The `tools/` scripts re-ran failed/overflowing replicates during development.
A clean run of `02` and `03b` should not need them; they are kept only for record.
