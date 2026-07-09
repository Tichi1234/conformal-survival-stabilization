## =====================================================================
## source-code.R -- master loader for the simulation code
## Run all scripts from the repository root using: source("R/source-code.R")
## =====================================================================

suppressPackageStartupMessages({
  library(survival)
})

source("R/simu.R")
source("R/Model-script.R")
source("R/calibration.R")
source("R/metrics.R")
source("R/run-sim.R")
