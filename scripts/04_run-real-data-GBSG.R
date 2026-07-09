## =====================================================================
## scripts/04_run-real-data-GBSG.R
## GBSG real-data application over 100 random splits. Writes CSVs to results/.
## Run from the repository root:  source("scripts/04_run-real-data-GBSG.R")
## =====================================================================

source("R/real-data-GBSG.R")

gbsg_res <- run_gbsg_real_data(
  R = 100,
  seed = 20260615,
  output_dir = "results/gbsg",
  time_scale = "years"
)

dir.create("results/gbsg", showWarnings = FALSE, recursive = TRUE)
write.csv(gbsg_res$summary, "results/gbsg/gbsg_realdata_summary_R100.csv", row.names = FALSE)
write.csv(gbsg_res$raw, "results/gbsg/gbsg_realdata_raw_R100.csv", row.names = FALSE)

cat("\nGBSG real-data analysis complete.\n")
cat("Summary also saved as results/gbsg/gbsg_realdata_summary_R100.csv\n")
cat("Raw results also saved as results/gbsg/gbsg_realdata_raw_R100.csv\n")
cat("Detailed output folder:\n")
cat(gbsg_res$output_dir, "\n")
