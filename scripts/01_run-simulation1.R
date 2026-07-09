## =====================================================================
## scripts/01_run-simulation1.R
## Simulation 1: internal ablation benchmark. Writes CSVs to results/.
## Run from the repository root:  source("scripts/01_run-simulation1.R")
## =====================================================================

source("R/source-code.R")
out_dir <- "results/simulation1"
res1 <- simulation1(R = 100, cores = 1)
save_simulation_result(res1, out_dir, "simulation1")
