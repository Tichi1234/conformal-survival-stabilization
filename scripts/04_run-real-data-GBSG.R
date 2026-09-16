############################################################
## 04_run-real-data-GBSG.R
##
## Wrapper for the GBSG real-data application.
##
## GitHub placement:
##   scripts/04_run-real-data-GBSG.R
##
## Run from repository root:
##   source("scripts/04_run-real-data-GBSG.R")
############################################################

source("R/real-data-GBSG.R")

gbsg_res <- run_gbsg_real_data(
  R = 100,
  seed = 20260615,
  output_dir = "results/real_data_gbsg",
  time_scale = "years"
)

write.csv(gbsg_res$summary, "gbsg_realdata_summary_R100.csv", row.names = FALSE)
write.csv(gbsg_res$raw, "gbsg_realdata_raw_R100.csv", row.names = FALSE)

cat("\nGBSG real-data analysis complete.\n")
cat("Summary also saved as gbsg_realdata_summary_R100.csv\n")
cat("Raw results also saved as gbsg_realdata_raw_R100.csv\n")
cat("Detailed output folder:\n")
cat(gbsg_res$output_dir, "\n")
