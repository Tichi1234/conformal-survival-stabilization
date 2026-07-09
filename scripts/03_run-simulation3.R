## =====================================================================
## scripts/03_run-simulation3.R
## Simulation 3: internal ablation on the hardest cells. External competitors are in 03b.
## Run from the repository root:  source("scripts/03_run-simulation3.R")
## =====================================================================

source("R/source-code.R")
out_dir <- "results/simulation3"
res3 <- simulation3(R = 100, cores = 8)
save_simulation_result(res3, out_dir, "simulation3_internal_only")

## This driver runs the INTERNAL ablation (H-/Clip- IPCW and AIPCW) on the
## difficult Simulation 3 cells only. The external competitors (DR-COSARC
## fixed/adaptive, KM de-censoring, Oracle, and Sesia-Svetnik-DR) are produced
## by 03b_run-simulation3-external-benchmark.R, which calls the authors' code.
