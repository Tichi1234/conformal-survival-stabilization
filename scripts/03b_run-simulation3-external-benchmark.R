## =====================================================================
## scripts/03b_run-simulation3-external-benchmark.R
##
## Simulation 3 EXTERNAL benchmark, rebuilt cleanly on top of the R/ library
## and the authors' (Sesia & Svetnik) conf_surv utilities.
##
## For every difficult cell and replicate it evaluates, on identical
## train/calibration/test splits and seeds:
##   * internal ablation   : H-IPCW, Clip-IPCW, H-AIPCW, Clip-AIPCW
##                           (R/ library; carries weight diagnostics)
##   * external competitors : DR-COSARC (fixed), DR-COSARC (adaptive),
##                           KM de-censoring, Oracle  (Sesia code; no weights)
##
## Canonical external settings match Sesia's experiment_1.R:
##   DR-COSARC (fixed)    -> predict_drcosarc(cutoffs="candes-fixed", doubly_robust=TRUE)
##   DR-COSARC (adaptive) -> predict_drcosarc(cutoffs="adaptive",     doubly_robust=TRUE,
##                                            finite_sample_correction=FALSE)
##   KM de-censoring      -> predict_decensoring(R = 10)
##   Nuisances            -> SurvregModelWrapper(dist="lognormal"); KM on TRAIN fold.
##
## Requirements: survival, tidyverse, R6, and Sesia's conf_surv utilities at
##   external_methods_sesia/code/conf_surv/
## Output: simulation3/results/simulation3_external_{raw,summary}_R<R>.csv
##
## Run from the repository root:
##   source("scripts/03b_run-simulation3-external-benchmark.R")
## =====================================================================

suppressMessages({ library(survival); library(dplyr); library(tidyr); library(readr); library(purrr); library(tibble); library(ggplot2);  library(R6) })

## ---- 1. Internal library --------------------------------------------------------
source("R/source-code.R")

## ---- 2. Sesia conf_surv utilities (isolated env) --------------------------------
SESIA_DIR <- "external_methods_sesia/code/conf_surv"
sesia <- new.env()
for (f in c("utils_survival.R", "utils_censoring.R",
            "utils_conformal.R", "utils_decensoring.R")) {
  sys.source(file.path(SESIA_DIR, f), envir = sesia)
}

## ---- 3. Configuration -----------------------------------------------------------
R_MAIN    <- 100
KM_R      <- 10
SURV_DIST <- "lognormal"
INTERNAL  <- c("H-IPCW", "Clip-IPCW", "H-AIPCW", "Clip-AIPCW")

CELLS <- expand.grid(
  setting = c("B", "C"), mech = c("C2", "C3"),
  n = c(600, 1200), cens = c(0.60, 0.80),
  stringsAsFactors = FALSE
)

OUT_DIR <- "results/simulation3"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- 4. Adapter: R/ list -> Sesia data frame ------------------------------------
to_sesia_df <- function(d) data.frame(d$X, time = d$Y, status = as.logical(d$Delta))

## ---- 5. Oracle: true conditional alpha-quantile of T | X (matches simu.R) --------
true_oracle_lpb <- function(X, setting, alpha = ALPHA) {
  if (setting == "A") {
    mu <- 0.5*X$X1 - 0.5*X$X2 + 0.8*X$X3 - 0.6*X$X4
    return(as.numeric(exp(mu + 0.60 * log(-log(1 - alpha)))))
  }
  if (setting == "B") {
    mu <- 0.8*sin(X$X1) + 0.5*X$X2^2 - 0.7*X$X3*X$X4
    return(as.numeric(exp(mu + 0.60 * qnorm(alpha))))
  }
  if (setting == "C") {
    mu    <- 0.5*X$X1 - 0.4*X$X2 + 0.6*X$X3
    sigma <- 0.5 + 0.3*abs(X$X1)
    return(as.numeric(exp(mu + sigma * qnorm(alpha))))
  }
  stop("Unknown setting: ", setting)
}

## ---- 6. Evaluator for a bare lower-bound vector (same columns as eval_method) ----
eval_bound <- function(lpb, test, method, alpha = ALPHA, ev_lp = NULL) {
  lpb <- as.numeric(lpb)
  cov    <- mean(test$T >= lpb)
  short  <- max(0, (1 - alpha) - cov)
  abscal <- abs(cov - (1 - alpha))
  groups <- list(X3 = test$X$X3, X4 = test$X$X4)
  if (!is.null(ev_lp)) {
    brks <- quantile(ev_lp, probs = seq(0, 1, .25), na.rm = TRUE)
    if (length(unique(brks)) > 1) groups$LPq <- cut(ev_lp, breaks = brks, include.lowest = TRUE)
  }
  cov_g <- unlist(lapply(groups, function(g) tapply(test$T >= lpb, g, mean)))
  worst <- if (length(cov_g)) min(cov_g, na.rm = TRUE) else NA_real_
  data.frame(method = method, tau = NA_real_, coverage = cov, shortfall = short,
             abs_cal_err = abscal, mean_lpb = mean(lpb), med_lpb = median(lpb),
             mean_loglpb = mean(log(pmax(lpb, 1e-8))), worst_slice = worst,
             max_raw_wt = NA_real_, max_norm_wt = NA_real_, ess = NA_real_,
             wt_cv = NA_real_, clipfrac = NA_real_)
}

safe_bound <- function(label, n_test, expr) {
  out <- tryCatch(as.numeric(expr),
                  error = function(e) { message(sprintf("  [%s] %s", label, conditionMessage(e)));
                                        rep(NA_real_, n_test) })
  if (length(out) != n_test) out <- rep(NA_real_, n_test)
  out
}

## ---- 7. One replicate of one cell ----------------------------------------------
run_external_replicate <- function(n, setting, mech, cens,
                                   backend = "weibull_cox", seed = NULL, ntest = NTEST) {
  lam <- get_lambda0(setting, mech, cens)
  if (!is.null(seed)) set.seed(seed)

  obs   <- gen_data(n, setting, mech, lam)
  itr   <- sample.int(n, floor(0.6 * n))
  train <- subset_data(obs, itr)
  cal   <- subset_data(obs, setdiff(seq_len(n), itr))
  test  <- gen_data(ntest, setting, mech, lam)

  nu    <- fit_nuisances(train, backend = backend)
  ev_lp <- tryCatch(as.numeric(nu$qhat(test$X, 0.5)), error = function(e) NULL)

  ## internal ablation (identical splits/seed to Simulation 3)
  res <- run_internal_methods(nu, cal, TAUGRID, test$X, alpha = ALPHA,
                              methods = INTERNAL, aug_grid_size = 50)
  rows_int <- do.call(rbind, lapply(names(res), function(m)
    cbind(method = m, eval_method(res[[m]], test, nu, alpha = ALPHA, event_lp = ev_lp))))

  ## external competitors (Sesia)
  tr_df <- to_sesia_df(train); ca_df <- to_sesia_df(cal); te_df <- to_sesia_df(test)

  surv_model <- sesia$SurvregModelWrapper$new(dist = SURV_DIST)
  surv_model$fit(Surv(time, status) ~ ., data = tr_df)
  cens_model <- sesia$CensoringModel$new(model = sesia$SurvregModelWrapper$new(dist = SURV_DIST))
  cens_model$fit(data = tr_df)
  km_fit <- survival::survfit(Surv(time, status) ~ 1, data = tr_df)

  pred_fixed <- safe_bound("DR-COSARC (fixed)", ntest,
    sesia$predict_drcosarc(te_df, surv_model, cens_model, ca_df, ALPHA,
                           doubly_robust = TRUE, cutoffs = "candes-fixed"))
  pred_adapt <- safe_bound("DR-COSARC (adaptive)", ntest,
    sesia$predict_drcosarc(te_df, surv_model, cens_model, ca_df, ALPHA,
                           doubly_robust = TRUE, cutoffs = "adaptive",
                           finite_sample_correction = FALSE))
  pred_km <- safe_bound("KM de-censoring", ntest,
    sesia$predict_decensoring(te_df, surv_model, km_fit, ca_df, ALPHA, R = KM_R))
  pred_oracle <- true_oracle_lpb(test$X, setting, ALPHA)

  rows_ext <- rbind(
    eval_bound(pred_fixed,  test, "DR-COSARC (fixed)",    ALPHA, ev_lp),
    eval_bound(pred_adapt,  test, "DR-COSARC (adaptive)", ALPHA, ev_lp),
    eval_bound(pred_km,     test, "KM de-censoring",      ALPHA, ev_lp),
    eval_bound(pred_oracle, test, "Oracle",               ALPHA, ev_lp))

  cbind(setting = setting, mech = mech, n = n, cens = cens, rbind(rows_int, rows_ext))
}

## ---- 8. Driver ------------------------------------------------------------------
run_external_benchmark <- function(R = R_MAIN, ntest = NTEST, verbose = TRUE) {
  all_rows <- list()
  for (i in seq_len(nrow(CELLS))) {
    cell <- CELLS[i, ]
    if (verbose) message(sprintf("[cell %d/%d] setting=%s mech=%s n=%d cens=%.2f",
                                 i, nrow(CELLS), cell$setting, cell$mech, cell$n, cell$cens))
    for (r in seq_len(R)) {
      row <- tryCatch(
        cbind(rep = r, run_external_replicate(cell$n, cell$setting, cell$mech, cell$cens,
                                              seed = 1000 * r, ntest = ntest)),
        error = function(e) { message(sprintf("  replicate %d failed: %s", r, conditionMessage(e))); NULL })
      if (!is.null(row)) all_rows[[length(all_rows) + 1]] <- row
    }
  }
  do.call(rbind, all_rows)
}

## ---- 9. Run + save --------------------------------------------------------------
if (!exists("SOURCE_ONLY") || !isTRUE(SOURCE_ONLY)) {
  raw <- run_external_benchmark(R = R_MAIN)
  summary_tbl <- aggregate_metrics_by_scenario(raw, ALPHA)
  write.csv(raw, file.path(OUT_DIR, sprintf("simulation3_external_raw_R%d.csv", R_MAIN)), row.names = FALSE)
  write.csv(summary_tbl, file.path(OUT_DIR, sprintf("simulation3_external_summary_R%d.csv", R_MAIN)), row.names = FALSE)
  message("Saved raw and summary CSVs to ", OUT_DIR)
}
