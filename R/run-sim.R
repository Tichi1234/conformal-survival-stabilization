## =====================================================================
## run_simulations.R -- Simulation drivers (Section 3.5)
##
## Protocol (matches revised Section 3): for each replicate, generate an
## observed sample of size n, split into train/calibration, fit nuisances on
## train, calibrate on the calibration fold, and evaluate on an INDEPENDENT
## synthetic test set (true T known). R replicates per scenario.
##
## Usage:
##   source("R/run_simulations.R")
##   smoke_test()                 # quick check
##   res1 <- simulation1()        # internal ablation benchmark
##   res2 <- simulation2()        # stress test (heavier; consider R, cores)
## =====================================================================

ALPHA   <- 0.10
TAUGRID <- seq(0.01, 0.49, by = 0.01)
NTEST   <- 2000

## cache tuned lambda0 across replicates
.lam_cache <- new.env()
get_lambda0 <- function(setting, mech, cens) {
  key <- paste(setting, mech, cens, sep = "_")
  if (is.null(.lam_cache[[key]]))
    .lam_cache[[key]] <- tune_lambda0(setting, mech, cens)
  .lam_cache[[key]]
}

## ---- one replicate --------------------------------------------------------
run_one_replicate <- function(n, setting, mech, cens, backend = "weibull_cox",
                              p_train = 0.6, methods = NULL,
                              aug_grid_size = 50, seed = NULL) {
  lam <- get_lambda0(setting, mech, cens)
  if (!is.null(seed)) set.seed(seed)

  obs  <- gen_data(n, setting, mech, lam)
  itr  <- sample.int(n, floor(p_train * n))
  train <- subset_data(obs, itr)
  cal   <- subset_data(obs, setdiff(seq_len(n), itr))
  test  <- gen_data(NTEST, setting, mech, lam)

  nu <- fit_nuisances(train, backend = backend)
  if (is.null(methods))
    methods <- c("Naive-Y","CC","HT-IPCW","H-IPCW","Stab-IPCW",
                 "Clip-IPCW","H-AIPCW","Clip-AIPCW")

  ## event linear predictor on test for the grouped-coverage strata
  ev_lp <- tryCatch(as.numeric(nu$qhat(test$X, 0.5)), error = function(e) NULL)

  res <- run_internal_methods(nu, cal, TAUGRID, test$X, alpha = ALPHA,
                              methods = methods, aug_grid_size = aug_grid_size)
  rows <- lapply(names(res), function(m) {
    e <- eval_method(res[[m]], test, nu, alpha = ALPHA, event_lp = ev_lp)
    cbind(method = m, e)
  })
  do.call(rbind, rows)
}

## ---- scenario runner (R replicates) --------------------------------------
run_scenario <- function(n, setting, mech, cens, R = 100,
                         backend = "weibull_cox", cores = 1, ...) {
  one <- function(r) tryCatch(
    cbind(rep = r, n = n, setting = setting, mech = mech, cens = cens,
          backend = backend,
          run_one_replicate(n, setting, mech, cens, backend = backend,
                            seed = 1000 * r, ...)),
    error = function(e) NULL)
  if (cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    L <- parallel::mclapply(seq_len(R), one, mc.cores = cores)
  } else {
    L <- lapply(seq_len(R), one)
  }
  do.call(rbind, L)
}

## ---- Simulation 1: internal ablation benchmark ---------------------------
## Setting A, mechanism C1, n in {600,1200}, censoring 20/40/60.
simulation1 <- function(R = 100, cores = 1) {
  grid <- expand.grid(n = c(600, 1200), cens = c(0.20, 0.40, 0.60),
                      setting = "A", mech = "C1", stringsAsFactors = FALSE)
  out <- do.call(rbind, Map(function(n, cens, s, m)
    run_scenario(n, s, m, cens, R = R, cores = cores),
    grid$n, grid$cens, grid$setting, grid$mech))
  list(raw = out,
       summary = aggregate_metrics(out, ALPHA),
       scenario_summary = aggregate_metrics_by_scenario(out, ALPHA))
}

## ---- Simulation 2: stress test -------------------------------------------
## Settings A-C, mechanisms C1-C3, n in {300,600,1200}, censoring 40/60/80.
## Set backend = "rsf" for the flexible-nuisance comparison on selected cells.
simulation2 <- function(R = 100, cores = 8, backend = "weibull_cox",
                        settings = c("A","B","C"), mechs = c("C1","C2","C3"),
                        ns = c(300,600,1200), censs = c(0.40,0.60,0.80),
                        methods = c("H-IPCW","Stab-IPCW","Clip-IPCW",
                                    "H-AIPCW","Clip-AIPCW")) {
  grid <- expand.grid(n = ns, cens = censs, setting = settings, mech = mechs,
                      stringsAsFactors = FALSE)
  out <- do.call(rbind, Map(function(n, cens, s, m)
    run_scenario(n, s, m, cens, R = R, cores = cores, backend = backend,
                 methods = methods),
    grid$n, grid$cens, grid$setting, grid$mech))
  list(raw = out,
       summary = aggregate_metrics(out, ALPHA),
       scenario_summary = aggregate_metrics_by_scenario(out, ALPHA))
}

## ---- Simulation 3: external benchmark ------------------------------------
## Wire up fit_davidov() / fit_sesia_svetnik() before running; until then this
## evaluates the proposed methods on the selected difficult cells.
simulation3 <- function(R = 100, cores = 1) {
  grid <- expand.grid(n = c(600, 1200), cens = c(0.60, 0.80),
                      setting = c("B","C"), mech = c("C2","C3"),
                      stringsAsFactors = FALSE)
  out <- do.call(rbind, Map(function(n, cens, s, m)
    run_scenario(n, s, m, cens, R = R, cores = cores,
                 methods = c("H-IPCW","Clip-IPCW","H-AIPCW","Clip-AIPCW")),
    grid$n, grid$cens, grid$setting, grid$mech))
  list(raw = out,
       summary = aggregate_metrics(out, ALPHA),
       scenario_summary = aggregate_metrics_by_scenario(out, ALPHA))
}

## ---- utilities for saving outputs ---------------------------------------
save_simulation_result <- function(res, out_dir, prefix) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  write.csv(res$raw, file.path(out_dir, paste0(prefix, "_raw.csv")), row.names = FALSE)
  write.csv(res$summary, file.path(out_dir, paste0(prefix, "_overall_summary.csv")), row.names = FALSE)
  if (!is.null(res$scenario_summary))
    write.csv(res$scenario_summary, file.path(out_dir, paste0(prefix, "_scenario_summary.csv")), row.names = FALSE)
  invisible(out_dir)
}

## ---- smoke test -----------------------------------------------------------
smoke_test <- function() {
  cat("Running 3-replicate smoke test (Setting A, C1, 40% censoring, n=600)...\n")
  df <- run_scenario(600, "A", "C1", 0.40, R = 3, cores = 1)
  print(aggregate_metrics(df, ALPHA)[,
        c("method","mean_cov","pac_success","mean_abscal",
          "med_lpb","mean_maxwt","mean_ess")])
  invisible(df)
}
