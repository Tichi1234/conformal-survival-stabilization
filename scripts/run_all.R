## =====================================================================
## scripts/run_all.R
## Reproduce every result and figure in the paper, in order.
## Run from the repository root:  source("scripts/run_all.R")
##
## This is the FULL run (R = 100 per cell) and takes a long time. For a quick
## end-to-end check, run scripts/00_smoke-test.R instead, or lower R_MAIN inside
## each script.
## =====================================================================
message("== Simulation 1 =="); source("scripts/01_run-simulation1.R")
message("== Simulation 2 =="); source("scripts/02_run-simulation2.R")
message("== Simulation 3 (internal) =="); source("scripts/03_run-simulation3.R")
message("== Simulation 3 (external benchmark) =="); source("scripts/03b_run-simulation3-external-benchmark.R")
message("== GBSG real data =="); source("scripts/04_run-real-data-GBSG.R")
message("== Clipping sensitivity =="); source("scripts/06_clip-sensitivity.R")
message("== Figures =="); source("scripts/05_make-figures.R")
message("All experiments and figures complete. See results/ and figures/.")
