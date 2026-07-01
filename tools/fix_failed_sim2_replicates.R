## fix_failed_sim2_replicates.R
## Run from the GitHub simulation folder, i.e. the folder containing README.md and R/.
## It detects missing/incomplete Simulation 2 scenario-replicates, reruns only those cells,
## appends the replacements, and writes fixed raw + summary CSV files.

source("R/source-code.R")

RAW_FILE <- "simulation2_full_raw_R100.csv"
TARGET_R <- 100
EXPECTED_METHODS <- c("H-IPCW", "Stab-IPCW", "Clip-IPCW", "H-AIPCW", "Clip-AIPCW")
BACKEND <- "weibull_cox"

if (!file.exists(RAW_FILE)) {
  stop("Cannot find ", RAW_FILE, ". Put this script in the folder where the raw Simulation 2 CSV is saved, or change RAW_FILE.")
}

raw <- read.csv(RAW_FILE, stringsAsFactors = FALSE)

## Expected Simulation 2 design
expected <- expand.grid(
  setting = c("A", "B", "C"),
  mech = c("C1", "C2", "C3"),
  n = c(300, 600, 1200),
  cens = c(0.40, 0.60, 0.80),
  rep = seq_len(TARGET_R),
  stringsAsFactors = FALSE
)

make_key <- function(d) {
  paste(d$setting, d$mech, d$n, sprintf("%.2f", as.numeric(d$cens)), d$rep, sep = "_")
}

counts <- aggregate(
  method ~ setting + mech + n + cens + rep,
  data = raw,
  FUN = function(x) length(unique(x))
)
names(counts)[names(counts) == "method"] <- "n_methods"

expected$key <- make_key(expected)
counts$key <- make_key(counts)
expected$n_methods <- counts$n_methods[match(expected$key, counts$key)]

missing <- expected[is.na(expected$n_methods) | expected$n_methods < length(EXPECTED_METHODS), ]
missing <- missing[order(missing$setting, missing$mech, missing$n, missing$cens, missing$rep), ]

cat("Number of missing/incomplete scenario-replicates:", nrow(missing), "\n")
if (nrow(missing) == 0) {
  cat("Nothing to fix. Writing a copy of the current files.\n")
  write.csv(raw, "simulation2_full_raw_R100_fixed.csv", row.names = FALSE)
  write.csv(aggregate_metrics(raw, ALPHA), "simulation2_full_summary_R100_fixed.csv", row.names = FALSE)
  quit(save = "no")
}

print(missing[, c("setting", "mech", "n", "cens", "rep", "n_methods")])

## Remove any incomplete rows for these replicate keys before appending replacements.
raw$key <- make_key(raw)
raw_clean <- raw[!(raw$key %in% missing$key), ]
raw_clean$key <- NULL

extras <- list()
for (i in seq_len(nrow(missing))) {
  s <- missing[i, ]
  seed_i <- 8800000 + i * 1000 + s$rep
  cat("\nRerunning replacement", i, "of", nrow(missing), ":",
      "setting=", s$setting,
      "mech=", s$mech,
      "n=", s$n,
      "cens=", s$cens,
      "rep=", s$rep,
      "seed=", seed_i, "\n")

  one <- run_one_replicate(
    n = s$n,
    setting = s$setting,
    mech = s$mech,
    cens = s$cens,
    backend = BACKEND,
    methods = EXPECTED_METHODS,
    seed = seed_i
  )

  one <- cbind(
    rep = s$rep,
    n = s$n,
    setting = s$setting,
    mech = s$mech,
    cens = s$cens,
    backend = BACKEND,
    one
  )
  extras[[i]] <- one
}

extra_df <- do.call(rbind, extras)
fixed <- rbind(raw_clean, extra_df)

## Sort for readability
fixed <- fixed[order(fixed$setting, fixed$mech, fixed$n, fixed$cens, fixed$rep, fixed$method), ]

## Verify all scenario-replicates are complete after replacement.
fixed_counts <- aggregate(
  method ~ setting + mech + n + cens + rep,
  data = fixed,
  FUN = function(x) length(unique(x))
)
names(fixed_counts)[names(fixed_counts) == "method"] <- "n_methods"
fixed_counts$key <- make_key(fixed_counts)
expected$n_methods_fixed <- fixed_counts$n_methods[match(expected$key, fixed_counts$key)]
still_bad <- expected[is.na(expected$n_methods_fixed) | expected$n_methods_fixed < length(EXPECTED_METHODS), ]

cat("\nAfter fixing, remaining incomplete scenario-replicates:", nrow(still_bad), "\n")
if (nrow(still_bad) > 0) {
  print(still_bad[, c("setting", "mech", "n", "cens", "rep", "n_methods_fixed")])
  warning("Some cells are still incomplete. Inspect the printed rows.")
}

write.csv(fixed, "simulation2_full_raw_R100_fixed.csv", row.names = FALSE)
write.csv(aggregate_metrics(fixed, ALPHA), "simulation2_full_summary_R100_fixed.csv", row.names = FALSE)

cat("\nWrote:\n")
cat("  simulation2_full_raw_R100_fixed.csv\n")
cat("  simulation2_full_summary_R100_fixed.csv\n")
cat("\nUse the *_fixed.csv files for final analysis.\n")
