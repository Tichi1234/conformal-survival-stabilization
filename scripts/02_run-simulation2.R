## =====================================================================
## scripts/02_run-simulation2.R
## Simulation 2: stress test (Settings A-C, mechanisms C1-C3). Writes CSVs to results/.
## Run from the repository root:  source("scripts/02_run-simulation2.R")
## =====================================================================

source("R/source-code.R")
out_dir <- "results/simulation2"
res2 <- simulation2(R = 100, cores = 8)
save_simulation_result(res2, out_dir, "simulation2_weibull_cox")

## Optional flexible-nuisance check; requires randomForestSRC.
## res2_rsf <- simulation2(R = 100, backend = "rsf", settings = "B", mechs = "C2",
##                         ns = 1200, censs = 0.60, cores = 1)
## save_simulation_result(res2_rsf, out_dir, "simulation2_rsf_selected")
