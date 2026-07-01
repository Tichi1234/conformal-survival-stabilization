############################################################
## Repair Simulation 3 external benchmark to recover 1600
## valid scenario-replicates per method
##
## What this does:
## 1. Reads simulation3/results/simulation3_external_raw_R100.csv
## 2. Identifies failed/overflow scenario-replicates
## 3. Reruns only those scenario-replicates
## 4. If the same seed fails again, uses replacement seeds
## 5. Replaces bad rows and writes a repaired raw file
## 6. Creates a clean summary WITHOUT mean_lpb
##
## Run from the project folder:
## source("repair_simulation3_external_R100_to_1600.R")
############################################################

rm(list = ls())

## Run this script from the project root folder.
## This avoids Windows/Linux path problems.

############################################################
## 0. Controls
############################################################

ALPHA <- 0.10
NTEST <- 2000
TAUGRID <- seq(0.01, 0.49, by = 0.01)

SETTINGS <- c("B", "C")
MECHS <- c("C2", "C3")

OUR_METHODS <- c(
  "H-IPCW",
  "Clip-IPCW",
  "H-AIPCW",
  "Clip-AIPCW"
)

ALL_METHODS <- c(
  "H-IPCW",
  "Clip-IPCW",
  "H-AIPCW",
  "Clip-AIPCW",
  "DR-COSARC (fixed)",
  "DR-COSARC (adaptive)",
  "KM de-censoring",
  "Oracle"
)

KM_IMPUTE_R <- 10
MAX_REPAIR_ATTEMPTS <- 50

RAW_FILE <- "simulation3/results/simulation3_external_raw_R100.csv"
RAW_REPAIRED <- "simulation3/results/simulation3_external_raw_R100_repaired.csv"
SUMMARY_REPAIRED <- "simulation3/results/simulation3_external_summary_R100_repaired_no_mean_lpb.csv"
REPAIR_LOG <- "simulation3/results/simulation3_external_repair_log.csv"
ERR_OUT <- "simulation3/logs/simulation3_external_repair_errors.txt"

dir.create("simulation3/results", recursive = TRUE, showWarnings = FALSE)
dir.create("simulation3/logs", recursive = TRUE, showWarnings = FALSE)

if (file.exists(ERR_OUT)) file.remove(ERR_OUT)

############################################################
## 1. Packages
############################################################

pkgs <- c("survival", "tidyverse", "R6")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

suppressPackageStartupMessages({
  library(survival)
  library(tidyverse)
  library(R6)
})

############################################################
## 2. Source project code and Sesia code
############################################################

our_env <- new.env(parent = globalenv())
sesia_env <- new.env(parent = globalenv())

source_file_into_env <- function(fname, env, required = TRUE) {
  candidates <- c(fname, file.path("R", fname), file.path("scripts", fname))
  candidates <- candidates[file.exists(candidates)]

  if (length(candidates) == 0) {
    hits <- list.files(
      path = ".",
      pattern = paste0("^", gsub("\\.", "\\\\.", fname), "$"),
      recursive = TRUE,
      full.names = TRUE
    )
    hits <- hits[!grepl("external_methods", hits, ignore.case = TRUE)]
    candidates <- hits
  }

  if (length(candidates) == 0) {
    if (required) {
      stop(
        "Cannot find file: ", fname, "\n",
        "Copy ", fname, " into the project root, R/, or scripts/ folder.\n",
        "Current working directory: ", getwd()
      )
    }
    return(invisible(NULL))
  }

  hit <- candidates[order(nchar(candidates))][1]
  source(hit, local = env)
  cat("Sourced:", hit, "\n")
}

source_file_into_env("dgp.R", our_env)
source_file_into_env("nuisance.R", our_env)
source_file_into_env("calibration.R", our_env)
source_file_into_env("metrics.R", our_env, required = FALSE)

required_our_functions <- c(
  "tune_lambda0",
  "gen_data",
  "subset_data",
  "fit_nuisances",
  "run_internal_methods"
)

missing_our <- required_our_functions[!required_our_functions %in% ls(our_env)]
if (length(missing_our) > 0) {
  stop("Missing our functions: ", paste(missing_our, collapse = ", "))
}

sesia_dir <- "external_methods_sesia/code/conf_surv"

sesia_files <- c(
  "misc.R",
  "utils_data.R",
  "utils_survival.R",
  "utils_censoring.R",
  "utils_decensoring.R",
  "utils_conformal.R",
  "utils_plotting.R"
)

for (f in sesia_files) {
  source(file.path(sesia_dir, f), local = sesia_env)
  cat("Sourced Sesia:", f, "\n")
}

required_sesia_functions <- c(
  "SurvregModelWrapper",
  "CensoringModel",
  "predict_drcosarc",
  "predict_decensoring"
)

missing_sesia <- required_sesia_functions[!required_sesia_functions %in% ls(sesia_env)]
if (length(missing_sesia) > 0) {
  stop("Missing Sesia functions/classes: ", paste(missing_sesia, collapse = ", "))
}

cat("\nAll code loaded successfully.\n")

############################################################
## 3. Helper functions
############################################################

.lambda_cache <- new.env(parent = emptyenv())

get_lambda0_cached <- function(setting, mech, cens) {
  key <- paste(setting, mech, cens, sep = "_")

  if (is.null(.lambda_cache[[key]])) {
    cat("\nTuning lambda0 for:", key, "\n")
    .lambda_cache[[key]] <- our_env$tune_lambda0(
      setting = setting,
      mech = mech,
      target_cens = cens,
      n_pilot = 2e5,
      seed = 1
    )
  }

  .lambda_cache[[key]]
}

to_sesia_df <- function(d) {
  data.frame(
    d$X,
    time = d$Y,
    status = as.logical(d$Delta)
  )
}

true_oracle_lpb <- function(X, setting, alpha = 0.10) {
  if (setting == "B") {
    mu <- 0.8 * sin(X$X1) +
      0.5 * X$X2^2 -
      0.7 * X$X3 * X$X4

    return(as.numeric(exp(mu + 0.60 * qnorm(alpha))))
  }

  if (setting == "C") {
    mu <- 0.5 * X$X1 -
      0.4 * X$X2 +
      0.6 * X$X3

    sigma <- 0.5 + 0.3 * abs(X$X1)

    return(as.numeric(exp(mu + sigma * qnorm(alpha))))
  }

  stop("Unknown setting: ", setting)
}

safe_lpb <- function(label, n_test, expr) {
  out <- tryCatch(
    {
      z <- eval.parent(substitute(expr))
      as.numeric(z)
    },
    error = function(e) {
      msg <- paste0("ERROR in ", label, ": ", conditionMessage(e))
      message(msg)
      cat(msg, file = ERR_OUT, append = TRUE, sep = "\n")
      rep(NA_real_, n_test)
    }
  )

  if (length(out) != n_test) {
    msg <- paste0(
      "ERROR in ", label, ": returned length ",
      length(out), " but expected ", n_test
    )
    message(msg)
    cat(msg, file = ERR_OUT, append = TRUE, sep = "\n")
    out <- rep(NA_real_, n_test)
  }

  out
}

is_bad_lpb <- function(x) {
  any(!is.finite(x)) || any(abs(x) > 1e6, na.rm = TRUE)
}

evaluate_lpb <- function(lpb, test, method, alpha = 0.10, tau = NA_real_, diag = NULL) {
  ok <- is.finite(lpb)

  if (sum(ok) == 0) {
    coverage <- NA_real_
    shortfall <- NA_real_
    abs_cal_err <- NA_real_
    mean_lpb <- NA_real_
    med_lpb <- NA_real_
    worst_slice <- NA_real_
  } else {
    coverage <- mean(test$T[ok] >= lpb[ok])
    shortfall <- max(0, (1 - alpha) - coverage)
    abs_cal_err <- abs(coverage - (1 - alpha))
    mean_lpb <- mean(lpb[ok])
    med_lpb <- median(lpb[ok])

    slice_values <- c()

    if ("X3" %in% names(test$X)) {
      slice_values <- c(
        slice_values,
        tapply(test$T[ok] >= lpb[ok], test$X$X3[ok], mean)
      )
    }

    if ("X4" %in% names(test$X)) {
      slice_values <- c(
        slice_values,
        tapply(test$T[ok] >= lpb[ok], test$X$X4[ok], mean)
      )
    }

    worst_slice <- if (length(slice_values) > 0) {
      min(slice_values, na.rm = TRUE)
    } else {
      NA_real_
    }
  }

  if (is.null(diag)) {
    diag <- list(
      max_raw_wt = NA_real_,
      max_norm_wt = NA_real_,
      ess = NA_real_,
      cv = NA_real_,
      clipfrac = NA_real_
    )
  }

  data.frame(
    method = method,
    tau = tau,
    coverage = coverage,
    shortfall = shortfall,
    abs_cal_err = abs_cal_err,
    mean_lpb = mean_lpb,
    med_lpb = med_lpb,
    worst_slice = worst_slice,
    max_raw_wt = diag$max_raw_wt,
    max_norm_wt = diag$max_norm_wt,
    ess = diag$ess,
    wt_cv = diag$cv,
    clipfrac = diag$clipfrac,
    stringsAsFactors = FALSE
  )
}

mean_or_na <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

median_or_na <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}

sd_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(NA_real_)
  sd(x)
}

mcse_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(NA_real_)
  sd(x) / sqrt(length(x))
}

summarise_results_no_mean_lpb <- function(df, alpha = 0.10) {
  df %>%
    group_by(method) %>%
    summarise(
      R = sum(!is.na(coverage)),
      mean_cov = mean_or_na(coverage),
      mcse_cov = mcse_or_na(coverage),
      med_cov = median_or_na(coverage),
      pac_success = mean_or_na(coverage >= (1 - alpha)),
      mean_shortfall = mean_or_na(shortfall),
      mean_abs_cal_err = mean_or_na(abs_cal_err),
      med_lpb = median_or_na(med_lpb),
      worst_slice = mean_or_na(worst_slice),
      sd_tau = sd_or_na(tau),
      sd_cov = sd_or_na(coverage),
      mean_max_wt = mean_or_na(max_raw_wt),
      mean_ess = mean_or_na(ess),
      mean_clipfrac = mean_or_na(clipfrac),
      .groups = "drop"
    ) %>%
    mutate(method = factor(method, levels = ALL_METHODS)) %>%
    arrange(method)
}

base_seed <- function(setting, mech, n, cens, rep_id) {
  100000 +
    rep_id +
    1000 * match(setting, SETTINGS) +
    100 * match(mech, MECHS) +
    n +
    round(100 * cens)
}

############################################################
## 4. One repaired replicate
############################################################

run_one_replicate <- function(setting, mech, n, cens, rep_id, repair_attempt = 0) {

  seed0 <- base_seed(setting, mech, n, cens, rep_id)

  ## attempt 0 uses original seed.
  ## attempt >0 uses replacement seed, same scenario.
  seed <- seed0 + repair_attempt * 10000000

  set.seed(seed)

  lambda0 <- get_lambda0_cached(setting, mech, cens)

  obs <- our_env$gen_data(n, setting, mech, lambda0)

  idx_train <- sample.int(n, floor(0.60 * n))
  idx_cal <- setdiff(seq_len(n), idx_train)

  train <- our_env$subset_data(obs, idx_train)
  cal <- our_env$subset_data(obs, idx_cal)
  test <- our_env$gen_data(NTEST, setting, mech, lambda0)

  actual_cens_obs <- mean(obs$Delta == 0)
  actual_cens_cal <- mean(cal$Delta == 0)
  actual_cens_test <- mean(test$Delta == 0)

  ##########################################################
  ## Our methods
  ##########################################################

  nu <- our_env$fit_nuisances(train, backend = "weibull_cox")

  internal <- our_env$run_internal_methods(
    nu = nu,
    cal = cal,
    tau_grid = TAUGRID,
    X_test = test$X,
    alpha = ALPHA,
    methods = OUR_METHODS,
    aug_grid_size = 50
  )

  ## Reject replicate if any internal LPB overflows.
  for (m in names(internal)) {
    if (is_bad_lpb(internal[[m]]$lpb)) {
      stop("Internal method produced nonfinite/extreme LPB: ", m)
    }
  }

  rows_ours <- bind_rows(lapply(names(internal), function(m) {
    evaluate_lpb(
      lpb = internal[[m]]$lpb,
      test = test,
      method = m,
      alpha = ALPHA,
      tau = internal[[m]]$tau,
      diag = internal[[m]]$diag
    )
  }))

  ##########################################################
  ## Sesia methods
  ##########################################################

  data_train <- to_sesia_df(train)
  data_cal <- to_sesia_df(cal)
  data_test <- to_sesia_df(test)

  surv_model <- sesia_env$SurvregModelWrapper$new(dist = "lognormal")
  surv_model$fit(survival::Surv(time, status) ~ ., data = data_train)

  cens_base <- sesia_env$SurvregModelWrapper$new(dist = "lognormal")
  cens_model <- sesia_env$CensoringModel$new(model = cens_base)
  cens_model$fit(data = data_train)

  km_fit <- survival::survfit(
    survival::Surv(time, status) ~ 1,
    data = data_train
  )

  pred_fixed <- safe_lpb(
    "DR-COSARC (fixed)",
    NTEST,
    sesia_env$predict_drcosarc(
      data.test = data_test,
      surv_model = surv_model,
      cens_imputator = cens_model,
      data.cal = data_cal,
      alpha = ALPHA,
      cutoffs = "candes-fixed",
      doubly_robust = TRUE
    )
  )

  pred_adaptive <- safe_lpb(
    "DR-COSARC (adaptive)",
    NTEST,
    sesia_env$predict_drcosarc(
      data.test = data_test,
      surv_model = surv_model,
      cens_imputator = cens_model,
      data.cal = data_cal,
      alpha = ALPHA,
      cutoffs = "adaptive",
      finite_sample_correction = FALSE,
      doubly_robust = TRUE
    )
  )

  pred_km <- safe_lpb(
    "KM de-censoring",
    NTEST,
    sesia_env$predict_decensoring(
      data.test = data_test,
      surv_model = surv_model,
      km_fit = km_fit,
      data.cal = data_cal,
      alpha = ALPHA,
      R = KM_IMPUTE_R
    )
  )

  pred_oracle <- true_oracle_lpb(
    X = test$X,
    setting = setting,
    alpha = ALPHA
  )

  sesia_preds <- list(
    "DR-COSARC (fixed)" = pred_fixed,
    "DR-COSARC (adaptive)" = pred_adaptive,
    "KM de-censoring" = pred_km,
    "Oracle" = pred_oracle
  )

  for (m in names(sesia_preds)) {
    if (is_bad_lpb(sesia_preds[[m]])) {
      stop("Sesia/oracle method produced nonfinite/extreme LPB: ", m)
    }
  }

  rows_sesia <- bind_rows(
    evaluate_lpb(pred_fixed, test, "DR-COSARC (fixed)", ALPHA),
    evaluate_lpb(pred_adaptive, test, "DR-COSARC (adaptive)", ALPHA),
    evaluate_lpb(pred_km, test, "KM de-censoring", ALPHA),
    evaluate_lpb(pred_oracle, test, "Oracle", ALPHA)
  )

  bind_rows(rows_ours, rows_sesia) %>%
    mutate(
      setting = setting,
      mech = mech,
      n = n,
      cens = cens,
      rep = rep_id,
      seed = seed,
      repaired = TRUE,
      repair_attempt = repair_attempt,
      actual_cens_obs = actual_cens_obs,
      actual_cens_cal = actual_cens_cal,
      actual_cens_test = actual_cens_test,
      .before = method
    )
}

############################################################
## 5. Identify bad keys from existing raw file
############################################################

if (!file.exists(RAW_FILE)) {
  stop("Cannot find raw file: ", RAW_FILE)
}

raw <- read_csv(RAW_FILE, show_col_types = FALSE)

key_cols <- c("setting", "mech", "n", "cens", "rep")

## Bad rows: failed rows, overflow rows, or scenario-replicates with fewer than 8 valid methods.
bad_keys_1 <- raw %>%
  filter(
    !is.finite(coverage) |
      !is.finite(med_lpb) |
      abs(med_lpb) > 1e6 |
      !is.finite(mean_lpb) |
      abs(mean_lpb) > 1e50
  ) %>%
  distinct(across(all_of(key_cols)))

bad_keys_2 <- raw %>%
  group_by(across(all_of(key_cols))) %>%
  summarise(
    n_methods = n_distinct(method[!is.na(coverage)]),
    .groups = "drop"
  ) %>%
  filter(n_methods < length(ALL_METHODS)) %>%
  select(all_of(key_cols))

bad_keys <- bind_rows(bad_keys_1, bad_keys_2) %>%
  distinct(across(all_of(key_cols))) %>%
  arrange(setting, mech, n, cens, rep)

cat("\nBad scenario-replicates identified:\n")
print(bad_keys)

if (nrow(bad_keys) == 0) {
  cat("\nNo bad keys found. The raw file appears complete.\n")
  repaired_raw_no_mean <- raw %>% select(-mean_lpb)
  summary_repaired <- summarise_results_no_mean_lpb(raw, alpha = ALPHA)
  write_csv(repaired_raw_no_mean, RAW_REPAIRED)
  write_csv(summary_repaired, SUMMARY_REPAIRED)
  quit(save = "no")
}

############################################################
## 6. Repair only the bad keys
############################################################

repair_rows <- vector("list", nrow(bad_keys))
repair_log <- vector("list", nrow(bad_keys))

for (i in seq_len(nrow(bad_keys))) {
  bk <- bad_keys[i, ]

  cat("\n============================================================\n")
  cat("Repairing bad key", i, "of", nrow(bad_keys), "\n")
  print(bk)
  cat("============================================================\n")

  repaired_i <- NULL
  success_attempt <- NA_integer_
  success_seed <- NA_real_
  success <- FALSE
  last_error <- NA_character_

  for (attempt in 0:MAX_REPAIR_ATTEMPTS) {
    cat("Attempt:", attempt, "\n")

    tmp <- tryCatch(
      {
        run_one_replicate(
          setting = bk$setting,
          mech = bk$mech,
          n = bk$n,
          cens = bk$cens,
          rep_id = bk$rep,
          repair_attempt = attempt
        )
      },
      error = function(e) {
        last_error <<- conditionMessage(e)
        message("Repair attempt failed: ", last_error)
        NULL
      }
    )

    if (!is.null(tmp)) {
      ## Check that all 8 methods are valid and no median LPB overflow.
      valid_tmp <- tmp %>%
        filter(
          method %in% ALL_METHODS,
          is.finite(coverage),
          is.finite(med_lpb),
          abs(med_lpb) <= 1e6
        )

      if (n_distinct(valid_tmp$method) == length(ALL_METHODS)) {
        repaired_i <- tmp
        success_attempt <- attempt
        success_seed <- unique(tmp$seed)[1]
        success <- TRUE
        cat("Repair succeeded at attempt", attempt, "\n")
        break
      }
    }
  }

  if (!success) {
    stop(
      "Could not repair key after ", MAX_REPAIR_ATTEMPTS, " attempts:\n",
      paste(capture.output(print(bk)), collapse = "\n"),
      "\nLast error: ", last_error
    )
  }

  repair_rows[[i]] <- repaired_i

  repair_log[[i]] <- bk %>%
    mutate(
      repaired = success,
      repair_attempt = success_attempt,
      repair_seed = success_seed
    )
}

repair_rows <- bind_rows(repair_rows)
repair_log <- bind_rows(repair_log)

############################################################
## 7. Replace bad rows, write repaired raw and summary
############################################################

raw_good <- raw %>%
  anti_join(bad_keys, by = key_cols)

raw_repaired <- bind_rows(raw_good, repair_rows) %>%
  arrange(setting, mech, n, cens, rep, factor(method, levels = ALL_METHODS))

## Check final completeness.
final_counts <- raw_repaired %>%
  group_by(method) %>%
  summarise(R = sum(!is.na(coverage)), .groups = "drop") %>%
  mutate(method = factor(method, levels = ALL_METHODS)) %>%
  arrange(method)

cat("\nFinal R counts by method:\n")
print(final_counts)

## Drop mean_lpb for the final raw analysis file.
raw_repaired_no_mean <- raw_repaired %>%
  select(-mean_lpb)

summary_repaired <- summarise_results_no_mean_lpb(raw_repaired, alpha = ALPHA)

write_csv(raw_repaired_no_mean, RAW_REPAIRED)
write_csv(summary_repaired, SUMMARY_REPAIRED)
write_csv(repair_log, REPAIR_LOG)

cat("\n============================================================\n")
cat("Repair complete.\n")
cat("Bad keys repaired:", nrow(bad_keys), "\n")
cat("Repaired raw file:", RAW_REPAIRED, "\n")
cat("Repaired summary file:", SUMMARY_REPAIRED, "\n")
cat("Repair log:", REPAIR_LOG, "\n")
cat("============================================================\n\n")

print(summary_repaired)
