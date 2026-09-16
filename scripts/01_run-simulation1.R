source("R/source-code.R")
out_dir <- file.path("results", format(Sys.time(), "%Y%m%d_%H%M%S"))
res1 <- simulation1(R = 100, cores = 1)
save_simulation_result(res1, out_dir, "simulation1")
