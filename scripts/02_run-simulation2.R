source("R/source-code.R")
out_dir <- file.path("results", format(Sys.time(), "%Y%m%d_%H%M%S"))
res2 <- simulation2(R = 100, cores = 4)
save_simulation_result(res2, out_dir, "simulation2_weibull_cox")

## Optional flexible-nuisance check; requires randomForestSRC.
## res2_rsf <- simulation2(R = 100, backend = "rsf", settings = "B", mechs = "C2",
##                         ns = 1200, censs = 0.60, cores = 1)
## save_simulation_result(res2_rsf, out_dir, "simulation2_rsf_selected")
