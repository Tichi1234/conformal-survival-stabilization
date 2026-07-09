## =====================================================================
## scripts/00_smoke-test.R
## Quick 3-replicate end-to-end check of the pipeline.
## Run from the repository root:  source("scripts/00_smoke-test.R")
## =====================================================================

source("R/source-code.R")
smoke <- smoke_test()
