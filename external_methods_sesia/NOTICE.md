# Third-party code required for Simulation 3 (NOT included)

The Simulation 3 external benchmark (DR-COSARC fixed/adaptive, KM de-censoring)
calls the conformal-survival utilities of

  Sesia, M. and Svetnik, V. (2024). Doubly robust conformalized survival
  analysis with right-censored data. arXiv:2412.09729.

These files are NOT distributed here, because the authors' release did not
include a license, and code without a license may not be redistributed.

TO RUN scripts/03b_run-simulation3-external-benchmark.R:
  1. Obtain the authors' code from their official source.
  2. Place these files in this folder (external_methods_sesia/code/conf_surv/):
       utils_survival.R  utils_censoring.R  utils_conformal.R  utils_decensoring.R
  3. Then run the script from the repository root.

Scripts 00, 01, 02, 03, 04, 05, and 06 do NOT need this code and run as-is.
