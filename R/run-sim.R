## =====================================================================
## run_sim.R -- Simulation drivers
##
## Protocol: for each replicate, generate an observed sample of size n,
## split into train/calibration, fit nuisances on train,
## calibrate on the calibration fold, and evaluate on an INDEPENDENT
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
  if (!is.null(seed)) {set.seed(seed)
  }
  
  ## generate observed dataset
  obs  <- gen_data(n, setting, mech, lam)
  
  ## train/calibration split
  itr  <- sample.int(n, floor(p_train * n))
  train <- subset_data(obs, itr)
  cal   <- subset_data(obs, setdiff(seq_len(n), itr))
  
  ## independent synthetic test set
  ## true T is known here
  test  <- gen_data(NTEST, setting, mech, lam)

  ## fit nuisance models using TRAINING fold only
  nu <- fit_nuisances(train, backend = backend)
  
  ## default internal methods
  if (is.null(methods))
    methods <- c("Naive-Y","CC","HT-IPCW","H-IPCW","Stab-IPCW",
                 "Clip-IPCW","H-AIPCW","Clip-AIPCW")
  
  ## select clipping threshold from Training fold
  
  needs_clip <- any(
    c("Clip-IPCW", "Clip-AIPCW") %in% methods
  )
  
  if (needs_clip) {
    
    clip_sel <- select_clip_percentile(
      nu = nu,
      train = train,
      q = getOption("herg_clip_q", 0.90)
    )
    
    cclip <- clip_sel$cclip
    
  } else {
    
    cclip <- NULL
  }

  ## event-model score used only for grouped-coverage diagnostics
  ev_lp <- tryCatch(as.numeric(nu$qhat(test$X, 0.5)), error = function(e) NULL)

  ## calibrate on calibration fold
  
  res <- run_internal_methods(
    nu = nu,
    cal = cal,
    tau_grid = TAUGRID,
    X_test = test$X,
    alpha = ALPHA,
    methods = methods,
    aug_grid_size = aug_grid_size,
    cclip = cclip
  )
  
  ## evaluate on independent test set
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
    error = function(e) {
     message(
       sprintf(
       "Replicate %d failed: %s",
       r,
       conditionMessage(e)
       )
     )
     NULL
    }
  )
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
simulation2 <- function(R = 100, cores = 1, backend = "weibull_cox",
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

## =====================================================================
## Simulation 3: external benchmark against Sesia-Svetnik DR-COSARC
##
## Selected difficult cells:
##   settings B-C
##   mechanisms C2-C3
##   n = 600, 1200
##   censoring = 60%, 80%
##
## All methods use exactly the same generated sample,
## train/calibration split, and independent test sample.
## =====================================================================


## ---- one Simulation 3 replicate -------------------------------------

run_one_replicate_sim3 <- function(
  n,
  setting,
  mech,
  cens,
  backend = "weibull_cox",
  p_train = 0.6,
  aug_grid_size = 50,
  seed = NULL
) {

  lam <- get_lambda0(
    setting,
    mech,
    cens
  )

  if (!is.null(seed)) {
    set.seed(seed)
  }


  ## ---------------------------------------------------------------
  ## Generate observed sample
  ## ---------------------------------------------------------------

  obs <- gen_data(
    n,
    setting,
    mech,
    lam
  )


  ## ---------------------------------------------------------------
  ## Common train/calibration split
  ## ---------------------------------------------------------------

  itr <- sample.int(
    n,
    floor(p_train * n)
  )

  train <- subset_data(
    obs,
    itr
  )

  cal <- subset_data(
    obs,
    setdiff(
      seq_len(n),
      itr
    )
  )


  ## ---------------------------------------------------------------
  ## Common independent test sample
  ## ---------------------------------------------------------------

  test <- gen_data(
    NTEST,
    setting,
    mech,
    lam
  )


  ## ---------------------------------------------------------------
  ## Fit our nuisance models on TRAINING fold
  ## ---------------------------------------------------------------

  nu <- fit_nuisances(
    train,
    backend = backend
  )


  ## ---------------------------------------------------------------
  ## Select clipping threshold on TRAINING fold
  ## ---------------------------------------------------------------

  clip_sel <- select_clip_percentile(
    nu = nu,
    train = train,
    q = getOption(
      "herg_clip_q",
      0.90
    )
  )

  cclip <- clip_sel$cclip


  ## ---------------------------------------------------------------
  ## Internal methods
  ## ---------------------------------------------------------------

  internal_methods <- c(
    "H-IPCW",
    "Clip-IPCW",
    "H-AIPCW",
    "Clip-AIPCW"
  )

  res_internal <- run_internal_methods(
    nu = nu,
    cal = cal,
    tau_grid = TAUGRID,
    X_test = test$X,
    alpha = ALPHA,
    methods = internal_methods,
    aug_grid_size = aug_grid_size,
    cclip = cclip
  )

  ## ---------------------------------------------------------------
  ## Sesia-Svetnik adaptive DR-COSARC
  ##
  ## IMPORTANT:
  ## Sesia failures are isolated from the internal methods.
  ## A failed DR-COSARC fit is retained as an NA row with
  ## fit_ok = FALSE so numerical failure rates can be reported.
  ## ---------------------------------------------------------------

  sesia_attempt <- tryCatch(

    list(
      ok = TRUE,

      result = fit_sesia_svetnik(
        train = train,
        cal = cal,
        X_test = test$X,
        alpha = ALPHA,
        finite_sample_correction = FALSE
      ),

      error = NA_character_
    ),

    error = function(e) {

      list(
        ok = FALSE,
        result = NULL,
        error = conditionMessage(e)
      )
    }
  )


  ## ---------------------------------------------------------------
  ## Event-model score for grouped diagnostics
  ## ---------------------------------------------------------------

  ev_lp <- tryCatch(
    as.numeric(
      nu$qhat(
        test$X,
        0.5
      )
    ),
    error = function(e) NULL
  )


  ## ---------------------------------------------------------------
  ## Evaluate four internal methods
  ## ---------------------------------------------------------------

  rows <- lapply(
    names(res_internal),
    function(m) {

      e <- eval_method(
        res_internal[[m]],
        test,
        nu,
        alpha = ALPHA,
        event_lp = ev_lp
      )

      cbind(
        method = m,
        fit_ok = TRUE,
        error_msg = NA_character_,
        e
      )
    }
  )


  ## ---------------------------------------------------------------
  ## Evaluate DR-COSARC if successful
  ## ---------------------------------------------------------------

  if (sesia_attempt$ok) {

    e_sesia <- eval_method(
      sesia_attempt$result,
      test,
      nu,
      alpha = ALPHA,
      event_lp = ev_lp
    )

    rows[["DR-COSARC"]] <- cbind(
      method = "DR-COSARC",
      fit_ok = TRUE,
      error_msg = NA_character_,
      e_sesia
    )

  } else {

    message(
      paste0(
        "DR-COSARC failed: ",
        sesia_attempt$error
      )
    )


    ## Same metric structure as eval_method(), but NA because
    ## the external method did not return a valid prediction.
    e_sesia_fail <- data.frame(
      tau = NA_real_,
      coverage = NA_real_,
      shortfall = NA_real_,
      abs_cal_err = NA_real_,
      mean_lpb = NA_real_,
      med_lpb = NA_real_,
      mean_loglpb = NA_real_,
      worst_slice = NA_real_,
      max_raw_wt = NA_real_,
      max_norm_wt = NA_real_,
      ess = NA_real_,
      wt_cv = NA_real_,
      clipfrac = NA_real_,
      clip_level = NA_real_,
      nu_hat_c = NA_real_,
      runtime = NA_real_
    )

    rows[["DR-COSARC"]] <- cbind(
      method = "DR-COSARC",
      fit_ok = FALSE,
      error_msg = sesia_attempt$error,
      e_sesia_fail
    )
  }


  do.call(
    rbind,
    rows
  )
}

 ## ---- Simulation 3 scenario runner -----------------------------------

run_scenario_sim3 <- function(
  n,
  setting,
  mech,
  cens,
  R = 100,
  backend = "weibull_cox",
  cores = 1
) {

  one <- function(r) {

    tryCatch(

      cbind(
        rep = r,
        n = n,
        setting = setting,
        mech = mech,
        cens = cens,
        backend = backend,

        run_one_replicate_sim3(
          n = n,
          setting = setting,
          mech = mech,
          cens = cens,
          backend = backend,
          seed = 1000 * r
        )
      ),

      error = function(e) {

        message(
          sprintf(
            "Simulation 3 replicate %d failed: %s",
            r,
            conditionMessage(e)
          )
        )

        NULL
      }
    )
  }


  if (
    cores > 1 &&
    requireNamespace(
      "parallel",
      quietly = TRUE
    )
  ) {

    L <- parallel::mclapply(
      seq_len(R),
      one,
      mc.cores = cores
    )

  } else {

    L <- lapply(
      seq_len(R),
      one
    )
  }


  do.call(
    rbind,
    L
  )
}


## ---- Simulation 3 full benchmark ------------------------------------

simulation3 <- function(
  R = 100,
  cores = 1
) {

  grid <- expand.grid(
    n = c(
      600,
      1200
    ),

    cens = c(
      0.60,
      0.80
    ),

    setting = c(
      "B",
      "C"
    ),

    mech = c(
      "C2",
      "C3"
    ),

    stringsAsFactors = FALSE
  )


  out <- do.call(
    rbind,

    Map(
      function(
        n,
        cens,
        s,
        m
      ) {

        run_scenario_sim3(
          n = n,
          setting = s,
          mech = m,
          cens = cens,
          R = R,
          cores = cores
        )
      },

      grid$n,
      grid$cens,
      grid$setting,
      grid$mech
    )
  )


  list(
    raw = out,

    summary = aggregate_metrics(
      out,
      ALPHA
    ),

    scenario_summary =
      aggregate_metrics_by_scenario(
        out,
        ALPHA
      )
  )
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
  cat
  paste0(("Running 3-replicate smoke test,
          (Setting A, C1, 40% censoring, n=600)...\n")
  )
  df <- run_scenario(600, "A", "C1", 0.40, R = 3, cores = 1)
  sm <- aggregate_metrics(df, ALPHA)
  print(sm[, c(
    "method",
    "mean_cov",
    "cov_q05",
    "cov_q10",
    "pac_success_nominal",
    "pac_success_tol",
    "mean_abscal",
    "med_lpb",
    "mean_maxwt",
    "mean_ess",
    "mean_clipfrac",
    "mean_clip_level",
    "mean_nu_hat_c"
  )]
  )
  invisible(df)
}
