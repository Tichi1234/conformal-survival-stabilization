source("R/source-code.R")
out_dir <- file.path("results", format(Sys.time(), "%Y%m%d_%H%M%S"))
res3 <- simulation3(R = 100, cores = 1)
save_simulation_result(res3, out_dir, "simulation3_internal_only")

## External competitors are intentionally not approximated here.
## Wire Sesia-Svetnik-DR to the authors'
## official code before reporting them as external benchmark results.
