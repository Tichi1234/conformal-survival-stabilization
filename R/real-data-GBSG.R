############################################################
## real-data-GBSG.R
## Real-data application for clipped self-normalized IPCW/AIPCW
## conformal lower predictive bounds under right censoring.
##
## GitHub placement:
##   R/real-data-GBSG.R
##
## Run from repository root:
##   source("R/real-data-GBSG.R")
##   gbsg_res <- run_gbsg_real_data(R = 100)
############################################################

## -------------------------------------------------------------------------
## 0. Packages
## -------------------------------------------------------------------------

gbsg_required_pkgs <- c("survival", "TH.data", "dplyr", "readr", "tibble", "tidyr", "ggplot2")

gbsg_install_missing <- function(pkgs = gbsg_required_pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) {
    message("Installing missing packages: ", paste(miss, collapse = ", "))
    install.packages(miss)
  }
}

gbsg_load_packages <- function() {
  gbsg_install_missing()
  suppressPackageStartupMessages({
    library(survival)
    library(TH.data)
    library(dplyr)
    library(readr)
    library(tibble)
    library(tidyr)
    library(ggplot2)
  })
}

## -------------------------------------------------------------------------
## 1. User controls
## -------------------------------------------------------------------------

GBSG_ALPHA <- 0.10
GBSG_TARGET_COV <- 1 - GBSG_ALPHA
GBSG_TAU_GRID <- seq(0.01, 0.49, by = 0.01)

GBSG_METHODS <- c(
  "Naive-Y", "CC", "HT-IPCW", "H-IPCW",
  "Stab-IPCW", "Clip-IPCW", "H-AIPCW", "Clip-AIPCW"
)

GBSG_CLIP_RULE <- "q90"
GBSG_AUG_GRID_SIZE <- 60
GBSG_G_FLOOR <- 1e-4
GBSG_MIN_POS_TIME <- 1e-8

GBSG_TRAIN_PROP <- 0.50
GBSG_CAL_PROP <- 0.25
GBSG_TEST_PROP <- 0.25

GBSG_LPB_CAP_QUANTILE <- 0.995
GBSG_LPB_CAP_MULTIPLIER <- 1.25

GBSG_COLORS <- c(
  "Naive-Y"    = "#7F7F7F",
  "CC"         = "#1F77B4",
  "HT-IPCW"    = "#D62728",
  "H-IPCW"     = "#FF7F0E",
  "Stab-IPCW"  = "#9467BD",
  "Clip-IPCW"  = "#2CA02C",
  "H-AIPCW"    = "#8C564B",
  "Clip-AIPCW" = "#17BECF"
)

## -------------------------------------------------------------------------
## 2. Data loading and preprocessing
## -------------------------------------------------------------------------

load_gbsg2_data <- function(time_scale = c("years", "days")) {
  time_scale <- match.arg(time_scale)

  if (!requireNamespace("TH.data", quietly = TRUE)) {
    stop("Package 'TH.data' is required. Install it using install.packages('TH.data').")
  }

  data("GBSG2", package = "TH.data")
  dat <- get("GBSG2")

  required <- c(
    "time", "cens", "horTh", "age", "menostat",
    "tsize", "tgrade", "pnodes", "progrec", "estrec"
  )

  missing_cols <- setdiff(required, names(dat))
  if (length(missing_cols) > 0) {
    stop("GBSG2 is missing expected columns: ", paste(missing_cols, collapse = ", "))
  }

  dat <- dat[, required]

  dat <- dat %>%
    mutate(
      Y = if (time_scale == "years") time / 365.25 else time,
      Delta = as.integer(cens == 1),
      age = as.numeric(age),
      tsize = as.numeric(tsize),
      pnodes = as.numeric(pnodes),
      progrec = as.numeric(progrec),
      estrec = as.numeric(estrec),
      horTh = factor(horTh),
      menostat = factor(menostat),
      tgrade = factor(tgrade)
    ) %>%
    filter(
      is.finite(Y),
      Y > 0,
      !is.na(Delta),
      complete.cases(age, tsize, pnodes, progrec, estrec, horTh, menostat, tgrade)
    )

  dat <- as.data.frame(dat)

  attr(dat, "time_scale") <- time_scale
  dat
}

gbsg_covariates <- function() {
  c("age", "menostat", "tsize", "tgrade", "pnodes", "progrec", "estrec", "horTh")
}

gbsg_surv_formula <- function(covariates = gbsg_covariates(), event = "Delta") {
  as.formula(paste0("Surv(Y, ", event, ") ~ ", paste(covariates, collapse = " + ")))
}

gbsg_split_data <- function(dat,
                            train_prop = GBSG_TRAIN_PROP,
                            cal_prop = GBSG_CAL_PROP,
                            seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(dat)
  id <- sample(seq_len(n))

  n_train <- floor(train_prop * n)
  n_cal <- floor(cal_prop * n)

  train_id <- id[seq_len(n_train)]
  cal_id <- id[(n_train + 1):(n_train + n_cal)]
  test_id <- id[(n_train + n_cal + 1):n]

  list(
    train = dat[train_id, , drop = FALSE],
    cal   = dat[cal_id, , drop = FALSE],
    test  = dat[test_id, , drop = FALSE]
  )
}

## -------------------------------------------------------------------------
## 3. Nuisance models
## -------------------------------------------------------------------------

gbsg_fit_event_weibull <- function(train, covariates = gbsg_covariates()) {
  if (sum(train$Delta == 1, na.rm = TRUE) < 10) return(NULL)

  f <- gbsg_surv_formula(covariates, event = "Delta")

  fit <- tryCatch(
    suppressWarnings(
      survreg(
        f,
        data = train,
        dist = "weibull",
        control = survreg.control(maxiter = 500)
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)
  if (!is.finite(fit$scale) || fit$scale <= 0 || fit$scale > 20) return(NULL)
  if (any(!is.finite(stats::coef(fit)))) return(NULL)

  fit
}

gbsg_fit_censoring_cox <- function(train, covariates = gbsg_covariates()) {
  if (sum(train$Delta == 0, na.rm = TRUE) < 10) return(NULL)

  f <- gbsg_surv_formula(covariates, event = "I(1 - Delta)")

  fit <- tryCatch(
    suppressWarnings(
      coxph(
        f,
        data = train,
        x = TRUE,
        ties = "breslow",
        control = coxph.control(iter.max = 100)
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)
  if (any(!is.finite(stats::coef(fit)))) return(NULL)

  fit
}

gbsg_fit_censoring_km0 <- function(train) {
  if (sum(train$Delta == 0, na.rm = TRUE) < 2) return(NULL)
  tryCatch(survfit(Surv(Y, 1 - Delta) ~ 1, data = train), error = function(e) NULL)
}

gbsg_basehaz_df <- function(fitC) {
  bh <- tryCatch(basehaz(fitC, centered = FALSE), error = function(e) NULL)
  if (is.null(bh) || nrow(bh) == 0) return(data.frame(time = 0, hazard = 0))
  bh <- bh[is.finite(bh$time) & is.finite(bh$hazard), , drop = FALSE]
  bh[order(bh$time), , drop = FALSE]
}

gbsg_step_H0 <- function(fitC, time) {
  bh <- gbsg_basehaz_df(fitC)
  approx(
    x = c(0, bh$time),
    y = c(0, bh$hazard),
    xout = time,
    method = "constant",
    rule = 2,
    f = 0
  )$y
}

gbsg_predict_lp_cox_zero <- function(fitC, newdata) {
  tryCatch(
    as.numeric(predict(fitC, newdata = newdata, type = "lp", reference = "zero")),
    error = function(e) rep(NA_real_, nrow(newdata))
  )
}

gbsg_predict_G_cox <- function(fitC, newdata, time) {
  if (is.null(fitC)) return(rep(NA_real_, nrow(newdata)))
  time <- rep_len(time, nrow(newdata))
  H0 <- gbsg_step_H0(fitC, time)
  lp <- gbsg_predict_lp_cox_zero(fitC, newdata)
  G <- exp(-H0 * exp(lp))
  G[!is.finite(G)] <- NA_real_
  pmin(pmax(G, GBSG_G_FLOOR), 1)
}

gbsg_predict_G0_km <- function(fitKM0, time) {
  if (is.null(fitKM0)) return(rep(NA_real_, length(time)))
  tt <- fitKM0$time
  ss <- fitKM0$surv
  if (length(tt) == 0 || length(ss) == 0) return(rep(1, length(time)))
  G0 <- approx(
    x = c(0, tt),
    y = c(1, ss),
    xout = time,
    method = "constant",
    rule = 2,
    f = 0
  )$y
  pmin(pmax(G0, GBSG_G_FLOOR), 1)
}

gbsg_lpb_cap <- function(train) {
  y <- train$Y[is.finite(train$Y) & train$Y > 0]
  if (length(y) == 0) return(Inf)
  cap <- as.numeric(stats::quantile(y, GBSG_LPB_CAP_QUANTILE, na.rm = TRUE)) *
    GBSG_LPB_CAP_MULTIPLIER
  if (!is.finite(cap) || cap <= GBSG_MIN_POS_TIME) Inf else cap
}

gbsg_predict_q_weibull <- function(fitT, newdata, tau, lpb_cap = Inf) {
  if (is.null(fitT) || !is.finite(tau)) return(rep(NA_real_, nrow(newdata)))

  q <- tryCatch(
    as.numeric(predict(fitT, newdata = newdata, type = "quantile", p = tau)),
    error = function(e) rep(NA_real_, nrow(newdata))
  )

  q[!is.finite(q)] <- NA_real_
  q <- pmax(q, GBSG_MIN_POS_TIME)
  if (is.finite(lpb_cap)) q <- pmin(q, lpb_cap)
  q
}

## -------------------------------------------------------------------------
## 4. Weights and diagnostics
## -------------------------------------------------------------------------

gbsg_clip_prob <- function(rule) {
  switch(rule, q80 = 0.80, q90 = 0.90, q95 = 0.95, q99 = 0.99, none = NA_real_,
         stop("Unknown clipping rule: ", rule))
}

gbsg_get_clip_c <- function(cal, fitC, rule = GBSG_CLIP_RULE) {
  if (rule == "none") return(Inf)
  if (is.null(fitC)) return(Inf)
  GY <- gbsg_predict_G_cox(fitC, cal, cal$Y)
  omega <- 1 / GY
  omega_pos <- omega[cal$Delta == 1 & is.finite(omega) & omega > 0]
  if (length(omega_pos) == 0) return(Inf)
  as.numeric(stats::quantile(omega_pos, probs = gbsg_clip_prob(rule), na.rm = TRUE))
}

gbsg_raw_ipcw_weights <- function(dat, fitC, clip_c = Inf) {
  GY <- gbsg_predict_G_cox(fitC, dat, dat$Y)
  omega <- 1 / GY
  if (is.finite(clip_c)) omega <- pmin(omega, clip_c)
  omega[!is.finite(omega)] <- NA_real_
  dat$Delta * omega
}

gbsg_raw_stabilized_weights <- function(dat, fitC, fitKM0, clip_c = Inf) {
  GY <- gbsg_predict_G_cox(fitC, dat, dat$Y)
  G0Y <- gbsg_predict_G0_km(fitKM0, dat$Y)
  omega <- G0Y / GY
  if (is.finite(clip_c)) omega <- pmin(omega, clip_c)
  omega[!is.finite(omega)] <- NA_real_
  dat$Delta * omega
}

gbsg_hajek_normalize <- function(raw_w) {
  denom <- mean(raw_w, na.rm = TRUE)
  if (!is.finite(denom) || denom <= 0) return(rep(NA_real_, length(raw_w)))
  raw_w / denom
}

gbsg_method_uses_clip <- function(method) method %in% c("Clip-IPCW", "Clip-AIPCW")
gbsg_method_uses_censoring <- function(method) !method %in% c("Naive-Y", "CC")
gbsg_method_uses_aipcw <- function(method) method %in% c("H-AIPCW", "Clip-AIPCW")

gbsg_get_weights <- function(dat, fitC, fitKM0, method, clip_c = Inf, normalized = TRUE) {
  if (method == "Naive-Y") return(rep(1, nrow(dat)))
  if (method == "CC") return(dat$Delta)

  if (method == "HT-IPCW") return(gbsg_raw_ipcw_weights(dat, fitC, clip_c = Inf))

  if (method %in% c("H-IPCW", "H-AIPCW")) {
    raw_w <- gbsg_raw_ipcw_weights(dat, fitC, clip_c = Inf)
    return(if (normalized) gbsg_hajek_normalize(raw_w) else raw_w)
  }

  if (method == "Stab-IPCW") {
    raw_w <- gbsg_raw_stabilized_weights(dat, fitC, fitKM0, clip_c = Inf)
    return(if (normalized) gbsg_hajek_normalize(raw_w) else raw_w)
  }

  if (method %in% c("Clip-IPCW", "Clip-AIPCW")) {
    raw_w <- gbsg_raw_ipcw_weights(dat, fitC, clip_c = clip_c)
    return(if (normalized) gbsg_hajek_normalize(raw_w) else raw_w)
  }

  stop("Unknown method: ", method)
}

gbsg_weight_diagnostics <- function(cal, fitC, fitKM0, method, clip_c = Inf) {
  raw_w <- gbsg_get_weights(cal, fitC, fitKM0, method, clip_c = clip_c, normalized = FALSE)
  norm_w <- gbsg_get_weights(cal, fitC, fitKM0, method, clip_c = clip_c, normalized = TRUE)

  raw0 <- ifelse(is.finite(raw_w), raw_w, 0)
  norm0 <- ifelse(is.finite(norm_w), norm_w, 0)

  ess <- if (sum(raw0^2) > 0) sum(raw0)^2 / sum(raw0^2) else NA_real_
  cv_raw <- if (mean(raw0) > 0) stats::sd(raw0) / mean(raw0) else NA_real_

  clip_frac <- NA_real_
  if (gbsg_method_uses_clip(method) && is.finite(clip_c)) {
    GY <- gbsg_predict_G_cox(fitC, cal, cal$Y)
    omega <- 1 / GY
    idx <- cal$Delta == 1
    clip_frac <- if (sum(idx) > 0) mean(omega[idx] > clip_c, na.rm = TRUE) else NA_real_
  }

  data.frame(
    max_raw_weight = max(raw0, na.rm = TRUE),
    max_normalized_weight = max(norm0, na.rm = TRUE),
    weight_variance = stats::var(raw0, na.rm = TRUE),
    weight_cv = cv_raw,
    effective_sample_size = ess,
    clipped_fraction = clip_frac
  )
}

## -------------------------------------------------------------------------
## 5. Calibration moments and AIPCW augmentation
## -------------------------------------------------------------------------

gbsg_thin_grid <- function(u, max_grid = GBSG_AUG_GRID_SIZE) {
  u <- sort(unique(u[is.finite(u) & u > 0]))
  if (length(u) <= max_grid) return(u)
  idx <- unique(round(seq(1, length(u), length.out = max_grid)))
  sort(u[idx])
}

gbsg_naive_y_moment <- function(tau, cal, fitT, lpb_cap = Inf) {
  q <- gbsg_predict_q_weibull(fitT, cal, tau, lpb_cap = lpb_cap)
  Icov <- as.numeric(cal$Y >= q)
  mean(Icov - GBSG_TARGET_COV, na.rm = TRUE)
}

gbsg_cc_moment <- function(tau, cal, fitT, lpb_cap = Inf) {
  q <- gbsg_predict_q_weibull(fitT, cal, tau, lpb_cap = lpb_cap)
  idx <- cal$Delta == 1
  if (sum(idx, na.rm = TRUE) < 2) return(NA_real_)
  Icov <- as.numeric(cal$Y[idx] >= q[idx])
  mean(Icov, na.rm = TRUE) - GBSG_TARGET_COV
}

gbsg_ht_ipcw_moment <- function(tau, cal, fitT, fitC, lpb_cap = Inf) {
  q <- gbsg_predict_q_weibull(fitT, cal, tau, lpb_cap = lpb_cap)
  Icov <- as.numeric(cal$Y >= q)
  raw_w <- gbsg_raw_ipcw_weights(cal, fitC, clip_c = Inf)
  mean(raw_w * Icov, na.rm = TRUE) - GBSG_TARGET_COV
}

gbsg_weighted_centered_moment <- function(tau, cal, fitT, weights, lpb_cap = Inf) {
  q <- gbsg_predict_q_weibull(fitT, cal, tau, lpb_cap = lpb_cap)
  Icov <- as.numeric(cal$Y >= q)
  if (all(!is.finite(weights)) || all(!is.finite(Icov))) return(NA_real_)
  mean(weights * (Icov - GBSG_TARGET_COV), na.rm = TRUE)
}

gbsg_predict_S_weibull <- function(fitT, newdata, time) {
  if (is.null(fitT)) return(rep(NA_real_, nrow(newdata)))
  lp <- tryCatch(as.numeric(predict(fitT, newdata = newdata, type = "lp")),
                 error = function(e) rep(NA_real_, nrow(newdata)))
  time <- rep_len(time, length(lp))
  shape <- 1 / fitT$scale
  scale <- exp(lp)
  S <- exp(-(pmax(time, GBSG_MIN_POS_TIME) / scale)^shape)
  S[!is.finite(S)] <- NA_real_
  pmin(pmax(S, GBSG_MIN_POS_TIME), 1)
}

gbsg_make_aipcw_cache <- function(cal, fitT, fitC, clip_c = Inf, lpb_cap = Inf) {
  if (is.null(fitT) || is.null(fitC)) return(NULL)

  ncal <- nrow(cal)
  bh <- gbsg_basehaz_df(fitC)
  ugrid_full <- sort(unique(c(cal$Y[cal$Delta == 0],
                              bh$time[bh$time <= max(cal$Y, na.rm = TRUE)])))
  ugrid <- gbsg_thin_grid(ugrid_full, GBSG_AUG_GRID_SIZE)
  ugrid <- ugrid[is.finite(ugrid) & ugrid > 0]
  if (length(ugrid) == 0) return(NULL)

  K <- length(ugrid)
  lower_grid <- c(0, ugrid[-K])

  lpC <- gbsg_predict_lp_cox_zero(fitC, cal)
  lpT <- tryCatch(as.numeric(predict(fitT, newdata = cal, type = "lp")),
                  error = function(e) rep(NA_real_, ncal))
  if (all(!is.finite(lpC)) || all(!is.finite(lpT))) return(NULL)

  exp_lpC <- exp(lpC)
  shapeT <- 1 / fitT$scale
  scaleT <- exp(lpT)

  H0u <- gbsg_step_H0(fitC, ugrid)
  H0lower <- gbsg_step_H0(fitC, lower_grid)
  dH0 <- pmax(H0u - H0lower, 0)

  at_risk <- outer(cal$Y, lower_grid, FUN = ">=") * 1

  dN <- matrix(0, nrow = ncal, ncol = K)
  for (k in seq_len(K)) {
    dN[, k] <- as.numeric(cal$Delta == 0 & cal$Y > lower_grid[k] & cal$Y <= ugrid[k])
  }

  dM <- dN - at_risk * outer(exp_lpC, dH0)

  GU <- exp(-outer(exp_lpC, H0u))
  GU <- pmin(pmax(GU, GBSG_G_FLOOR), 1)

  omega_u <- 1 / GU
  if (is.finite(clip_c)) omega_u <- pmin(omega_u, clip_c)
  omega_u[!is.finite(omega_u)] <- NA_real_

  Umat <- matrix(ugrid, nrow = ncal, ncol = K, byrow = TRUE)
  ScaleMat <- matrix(scaleT, nrow = ncal, ncol = K, byrow = FALSE)

  ST_u <- exp(-(pmax(Umat, GBSG_MIN_POS_TIME) / ScaleMat)^shapeT)
  ST_u <- pmin(pmax(ST_u, GBSG_MIN_POS_TIME), 1)

  list(cal = cal, fitT = fitT, ugrid = ugrid, shapeT = shapeT,
       scaleT = scaleT, ST_u = ST_u, omega_u = omega_u,
       dM = dM, ncal = ncal, lpb_cap = lpb_cap)
}

gbsg_augmentation_term <- function(tau, cache) {
  if (is.null(cache)) return(NA_real_)

  cal <- cache$cal
  ncal <- cache$ncal
  K <- length(cache$ugrid)

  q <- gbsg_predict_q_weibull(cache$fitT, cal, tau, lpb_cap = cache$lpb_cap)
  if (all(!is.finite(q))) return(NA_real_)
  q <- pmax(q, GBSG_MIN_POS_TIME)

  Qmat <- matrix(q, nrow = ncal, ncol = K, byrow = FALSE)
  Umat <- matrix(cache$ugrid, nrow = ncal, ncol = K, byrow = TRUE)
  ScaleMat <- matrix(cache$scaleT, nrow = ncal, ncol = K, byrow = FALSE)

  ST_max <- exp(-(pmax(Qmat, Umat) / ScaleMat)^cache$shapeT)
  ST_max <- pmin(pmax(ST_max, GBSG_MIN_POS_TIME), 1)

  eta <- pmin(pmax(ST_max / cache$ST_u, 0), 1)

  val <- sum(cache$omega_u * (eta - GBSG_TARGET_COV) * cache$dM, na.rm = TRUE) / ncal
  if (!is.finite(val)) NA_real_ else val
}

gbsg_select_tau <- function(cal, fitT, fitC, fitKM0, method,
                            clip_rule = GBSG_CLIP_RULE,
                            lpb_cap = Inf) {
  if (is.null(fitT)) return(list(tau_hat = NA_real_, phi_at_tau = NA_real_, clip_c = NA_real_))
  if (gbsg_method_uses_censoring(method) && is.null(fitC)) {
    return(list(tau_hat = NA_real_, phi_at_tau = NA_real_, clip_c = NA_real_))
  }

  clip_c <- if (gbsg_method_uses_clip(method)) gbsg_get_clip_c(cal, fitC, rule = clip_rule) else Inf

  wH <- NULL
  if (method %in% c("H-IPCW", "Stab-IPCW", "Clip-IPCW", "H-AIPCW", "Clip-AIPCW")) {
    wH <- gbsg_get_weights(cal, fitC, fitKM0, method, clip_c = clip_c, normalized = TRUE)
  }

  cache <- NULL
  if (gbsg_method_uses_aipcw(method)) {
    cache <- gbsg_make_aipcw_cache(
      cal, fitT, fitC,
      clip_c = if (method == "Clip-AIPCW") clip_c else Inf,
      lpb_cap = lpb_cap
    )
  }

  phi <- rep(NA_real_, length(GBSG_TAU_GRID))

  for (j in seq_along(GBSG_TAU_GRID)) {
    tau <- GBSG_TAU_GRID[j]

    if (method == "Naive-Y") {
      phi[j] <- gbsg_naive_y_moment(tau, cal, fitT, lpb_cap)
    } else if (method == "CC") {
      phi[j] <- gbsg_cc_moment(tau, cal, fitT, lpb_cap)
    } else if (method == "HT-IPCW") {
      phi[j] <- gbsg_ht_ipcw_moment(tau, cal, fitT, fitC, lpb_cap)
    } else if (method %in% c("H-IPCW", "Stab-IPCW", "Clip-IPCW")) {
      phi[j] <- gbsg_weighted_centered_moment(tau, cal, fitT, wH, lpb_cap)
    } else if (method %in% c("H-AIPCW", "Clip-AIPCW")) {
      phi_ipcw <- gbsg_weighted_centered_moment(tau, cal, fitT, wH, lpb_cap)
      phi_aug <- gbsg_augmentation_term(tau, cache)
      phi[j] <- if (is.finite(phi_ipcw) && is.finite(phi_aug)) phi_ipcw + phi_aug else NA_real_
    } else {
      stop("Unknown method: ", method)
    }
  }

  phi2 <- phi
  phi2[!is.finite(phi2)] <- -Inf
  feasible <- which(cummin(phi2) >= 0)
  tau_hat <- if (length(feasible) == 0) min(GBSG_TAU_GRID) else max(GBSG_TAU_GRID[feasible])
  idx <- which.min(abs(GBSG_TAU_GRID - tau_hat))

  list(tau_hat = tau_hat, phi_at_tau = phi[idx], clip_c = clip_c)
}

## -------------------------------------------------------------------------
## 6. Real-data held-out evaluation
## -------------------------------------------------------------------------

gbsg_eval_ipcw_coverage <- function(test, lpb, fitC_eval) {
  GY <- gbsg_predict_G_cox(fitC_eval, test, test$Y)
  w <- test$Delta / GY
  Icov <- as.numeric(test$Y >= lpb)

  cov_ht <- mean(w * Icov, na.rm = TRUE)

  cov_hajek <- if (sum(w, na.rm = TRUE) > 0) {
    sum(w * Icov, na.rm = TRUE) / sum(w, na.rm = TRUE)
  } else {
    NA_real_
  }

  obs_cov_y <- mean(test$Y >= lpb, na.rm = TRUE)

  c(coverage_ipcw_ht = cov_ht,
    coverage_ipcw_hajek = cov_hajek,
    observed_y_coverage = obs_cov_y)
}

gbsg_worst_slice_coverage <- function(test, lpb, fitC_eval) {
  GY <- gbsg_predict_G_cox(fitC_eval, test, test$Y)
  w <- test$Delta / GY
  Icov <- as.numeric(test$Y >= lpb)

  group_list <- list(
    age_quartile = cut(test$age,
                       breaks = unique(stats::quantile(test$age, probs = seq(0, 1, 0.25), na.rm = TRUE)),
                       include.lowest = TRUE),
    nodes_group = cut(test$pnodes,
                      breaks = unique(stats::quantile(test$pnodes, probs = seq(0, 1, 0.25), na.rm = TRUE)),
                      include.lowest = TRUE),
    horTh = test$horTh,
    tgrade = test$tgrade,
    menostat = test$menostat
  )

  covs <- c()

  for (g in group_list) {
    levs <- unique(g)
    levs <- levs[!is.na(levs)]

    for (lev in levs) {
      idx <- which(g == lev)
      if (length(idx) >= 10 && sum(w[idx], na.rm = TRUE) > 0) {
        covs <- c(covs, sum(w[idx] * Icov[idx], na.rm = TRUE) / sum(w[idx], na.rm = TRUE))
      }
    }
  }

  if (length(covs) == 0) NA_real_ else min(covs, na.rm = TRUE)
}

gbsg_run_single_method <- function(train, cal, test,
                                   fitT, fitC, fitKM0, fitC_eval,
                                   method,
                                   clip_rule = GBSG_CLIP_RULE,
                                   lpb_cap = Inf) {
  start_time <- proc.time()[3]

  tau_res <- gbsg_select_tau(cal, fitT, fitC, fitKM0, method,
                             clip_rule = clip_rule, lpb_cap = lpb_cap)

  lpb <- gbsg_predict_q_weibull(fitT, test, tau_res$tau_hat, lpb_cap = lpb_cap)

  covs <- gbsg_eval_ipcw_coverage(test, lpb, fitC_eval)
  worst <- gbsg_worst_slice_coverage(test, lpb, fitC_eval)
  wd <- gbsg_weight_diagnostics(cal, fitC, fitKM0, method, clip_c = tau_res$clip_c)

  coverage_primary <- unname(covs["coverage_ipcw_hajek"])

  data.frame(
    method = method,
    tau_hat = tau_res$tau_hat,
    calibration_moment = tau_res$phi_at_tau,
    clip_c = tau_res$clip_c,
    coverage_ipcw_hajek = coverage_primary,
    coverage_ipcw_ht = unname(covs["coverage_ipcw_ht"]),
    observed_y_coverage = unname(covs["observed_y_coverage"]),
    coverage_shortfall = max(0, GBSG_TARGET_COV - coverage_primary),
    abs_calibration_error = abs(coverage_primary - GBSG_TARGET_COV),
    pac_success = as.integer(coverage_primary >= GBSG_TARGET_COV),
    mean_lpb = mean(lpb, na.rm = TRUE),
    median_lpb = stats::median(lpb, na.rm = TRUE),
    mean_log_lpb = mean(log(pmax(lpb, GBSG_MIN_POS_TIME)), na.rm = TRUE),
    sd_lpb = stats::sd(lpb, na.rm = TRUE),
    worst_slice_coverage = worst,
    wd,
    computation_time_sec = unname(proc.time()[3] - start_time)
  )
}

## -------------------------------------------------------------------------
## 7. Main real-data runner
## -------------------------------------------------------------------------

run_gbsg_real_data <- function(R = 100,
                               seed = 20260615,
                               output_dir = "results/real_data_gbsg",
                               time_scale = "years",
                               methods = GBSG_METHODS) {
  gbsg_load_packages()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out_dir <- file.path(output_dir, paste0("run_", stamp))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

  cat("\nGBSG real-data application\n")
  cat("Output directory:", out_dir, "\n\n")

  dat <- load_gbsg2_data(time_scale = time_scale)
  covariates <- gbsg_covariates()

  data_summary <- data.frame(
    n = nrow(dat),
    events = sum(dat$Delta == 1),
    censored = sum(dat$Delta == 0),
    event_fraction = mean(dat$Delta == 1),
    censoring_fraction = mean(dat$Delta == 0),
    median_followup = median(dat$Y),
    time_scale = time_scale
  )

  readr::write_csv(data_summary, file.path(out_dir, "gbsg_data_summary.csv"))

  all_rows <- list()

  for (r in seq_len(R)) {
    cat("GBSG split", r, "of", R, "\n")

    sp <- gbsg_split_data(dat, seed = seed + r)

    train <- sp$train
    cal <- sp$cal
    test <- sp$test

    fitT <- gbsg_fit_event_weibull(train, covariates)
    fitC <- gbsg_fit_censoring_cox(train, covariates)
    fitKM0 <- gbsg_fit_censoring_km0(train)

    ## Evaluation model uses train + calibration data only.
    fitC_eval <- gbsg_fit_censoring_cox(rbind(train, cal), covariates)

    lpb_cap <- gbsg_lpb_cap(train)

    split_rows <- lapply(methods, function(m) {
      tryCatch(
        gbsg_run_single_method(
          train = train,
          cal = cal,
          test = test,
          fitT = fitT,
          fitC = fitC,
          fitKM0 = fitKM0,
          fitC_eval = fitC_eval,
          method = m,
          clip_rule = GBSG_CLIP_RULE,
          lpb_cap = lpb_cap
        ),
        error = function(e) {
          message("Method failed in split ", r, ": ", m, " | ", e$message)
          data.frame(
            method = m,
            tau_hat = NA_real_,
            calibration_moment = NA_real_,
            clip_c = NA_real_,
            coverage_ipcw_hajek = NA_real_,
            coverage_ipcw_ht = NA_real_,
            observed_y_coverage = NA_real_,
            coverage_shortfall = NA_real_,
            abs_calibration_error = NA_real_,
            pac_success = NA_integer_,
            mean_lpb = NA_real_,
            median_lpb = NA_real_,
            mean_log_lpb = NA_real_,
            sd_lpb = NA_real_,
            worst_slice_coverage = NA_real_,
            max_raw_weight = NA_real_,
            max_normalized_weight = NA_real_,
            weight_variance = NA_real_,
            weight_cv = NA_real_,
            effective_sample_size = NA_real_,
            clipped_fraction = NA_real_,
            computation_time_sec = NA_real_
          )
        }
      )
    })

    split_rows <- dplyr::bind_rows(split_rows)

    split_rows$split <- r
    split_rows$n_total <- nrow(dat)
    split_rows$n_train <- nrow(train)
    split_rows$n_cal <- nrow(cal)
    split_rows$n_test <- nrow(test)
    split_rows$events_train <- sum(train$Delta == 1)
    split_rows$events_cal <- sum(cal$Delta == 1)
    split_rows$events_test <- sum(test$Delta == 1)
    split_rows$censoring_test <- mean(test$Delta == 0)

    all_rows[[length(all_rows) + 1]] <- split_rows

    if (r %% 10 == 0) {
      partial <- dplyr::bind_rows(all_rows)
      readr::write_csv(partial, file.path(out_dir, "gbsg_partial_raw_results.csv"))
    }
  }

  raw <- dplyr::bind_rows(all_rows)

  summary <- raw %>%
    group_by(method) %>%
    summarise(
      R = sum(is.finite(coverage_ipcw_hajek)),
      mean_cov = mean(coverage_ipcw_hajek, na.rm = TRUE),
      mcse_cov = stats::sd(coverage_ipcw_hajek, na.rm = TRUE) / sqrt(R),
      med_cov = stats::median(coverage_ipcw_hajek, na.rm = TRUE),
      pac_success = mean(pac_success, na.rm = TRUE),
      mean_short = mean(coverage_shortfall, na.rm = TRUE),
      mean_abscal = mean(abs_calibration_error, na.rm = TRUE),
      med_lpb = stats::median(median_lpb, na.rm = TRUE),
      mean_lpb = mean(mean_lpb, na.rm = TRUE),
      worst_slice = mean(worst_slice_coverage, na.rm = TRUE),
      sd_tau = stats::sd(tau_hat, na.rm = TRUE),
      sd_cov = stats::sd(coverage_ipcw_hajek, na.rm = TRUE),
      mean_maxwt = mean(max_raw_weight, na.rm = TRUE),
      mean_ess = mean(effective_sample_size, na.rm = TRUE),
      mean_clipfrac = mean(clipped_fraction, na.rm = TRUE),
      mean_time_sec = mean(computation_time_sec, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(method = factor(method, levels = GBSG_METHODS)) %>%
    arrange(method)

  readr::write_csv(raw, file.path(out_dir, "gbsg_realdata_raw_R100.csv"))
  readr::write_csv(summary, file.path(out_dir, "gbsg_realdata_summary_R100.csv"))

  gbsg_make_figures(raw, summary, out_dir)

  sink(file.path(out_dir, "session_info.txt"))
  cat("GBSG real-data application completed at:\n")
  print(Sys.time())
  cat("\nData summary:\n")
  print(data_summary)
  cat("\nControls:\n")
  print(list(
    R = R,
    alpha = GBSG_ALPHA,
    target_coverage = GBSG_TARGET_COV,
    tau_grid = GBSG_TAU_GRID,
    train_prop = GBSG_TRAIN_PROP,
    cal_prop = GBSG_CAL_PROP,
    test_prop = GBSG_TEST_PROP,
    clip_rule = GBSG_CLIP_RULE,
    time_scale = time_scale,
    methods = methods
  ))
  cat("\nSummary:\n")
  print(summary)
  cat("\nSession info:\n")
  print(sessionInfo())
  sink()

  cat("\nDone. Main files written to:\n")
  cat(out_dir, "\n\n")
  print(summary)

  invisible(list(raw = raw, summary = summary, data_summary = data_summary, output_dir = out_dir))
}

## -------------------------------------------------------------------------
## 8. Figures
## -------------------------------------------------------------------------

gbsg_theme_pub <- function() {
  theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1),
      strip.background = element_rect(fill = "grey90", colour = "grey40"),
      panel.grid.minor = element_blank()
    )
}

gbsg_save_fig <- function(p, out_dir, name, width = 8.5, height = 5.8) {
  fig_dir <- file.path(out_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  ggsave(file.path(fig_dir, paste0(name, ".png")), plot = p,
         width = width, height = height, dpi = 300)

  ggsave(file.path(fig_dir, paste0(name, ".pdf")), plot = p,
         width = width, height = height)
}

gbsg_make_figures <- function(raw, summary, out_dir) {
  summary <- summary %>% mutate(method = factor(as.character(method), levels = GBSG_METHODS))
  raw <- raw %>% mutate(method = factor(as.character(method), levels = GBSG_METHODS))

  p1 <- ggplot(summary, aes(x = med_lpb, y = mean_cov, colour = method, size = mean_ess)) +
    geom_hline(yintercept = GBSG_TARGET_COV, linetype = "dashed", linewidth = 0.8) +
    geom_point(alpha = 0.90) +
    geom_text(aes(label = method), size = 3, hjust = -0.05, vjust = -0.4, show.legend = FALSE) +
    scale_colour_manual(values = GBSG_COLORS) +
    scale_size_continuous(range = c(3, 10)) +
    labs(
      title = "GBSG real-data application: coverage-informativeness trade-off",
      x = "Median lower predictive bound",
      y = "IPCW-estimated held-out coverage",
      size = "Mean ESS"
    ) +
    gbsg_theme_pub()

  gbsg_save_fig(p1, out_dir, "GBSG_Figure_1_coverage_informativeness")

  p2 <- ggplot(summary, aes(x = method, y = mean_maxwt, fill = method)) +
    geom_col(colour = "black", linewidth = 0.3, width = 0.75) +
    scale_fill_manual(values = GBSG_COLORS) +
    scale_y_log10() +
    labs(
      title = "GBSG real-data application: maximum inverse-censoring weights",
      x = "Method",
      y = "Mean maximum raw weight (log scale)"
    ) +
    gbsg_theme_pub() +
    theme(legend.position = "none")

  gbsg_save_fig(p2, out_dir, "GBSG_Figure_2_max_weight_logscale")

  dashboard <- summary %>%
    mutate(log10_mean_maxwt = log10(pmax(mean_maxwt, 1e-8))) %>%
    select(method, mean_cov, mean_short, mean_abscal, med_lpb,
           mean_ess, log10_mean_maxwt, mean_time_sec) %>%
    tidyr::pivot_longer(cols = -method, names_to = "metric", values_to = "value") %>%
    mutate(metric = dplyr::recode(
      metric,
      mean_cov = "Coverage",
      mean_short = "Shortfall",
      mean_abscal = "Abs. cal. error",
      med_lpb = "Median LPB",
      mean_ess = "ESS",
      log10_mean_maxwt = "log10 max weight",
      mean_time_sec = "Time"
    ))

  p3 <- ggplot(dashboard, aes(x = method, y = value, fill = method)) +
    geom_col(colour = "black", linewidth = 0.2, width = 0.75) +
    facet_wrap(~ metric, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = GBSG_COLORS) +
    labs(title = "GBSG real-data application: summary metrics", x = "Method", y = NULL) +
    gbsg_theme_pub() +
    theme(legend.position = "none")

  gbsg_save_fig(p3, out_dir, "GBSG_Figure_3_metric_dashboard", width = 11, height = 7)

  raw_long <- raw %>%
    mutate(log10_max_raw_weight = log10(pmax(max_raw_weight, 1e-8))) %>%
    select(method, coverage_ipcw_hajek, coverage_shortfall,
           abs_calibration_error, median_lpb, effective_sample_size,
           log10_max_raw_weight) %>%
    tidyr::pivot_longer(cols = -method, names_to = "metric", values_to = "value") %>%
    mutate(metric = dplyr::recode(
      metric,
      coverage_ipcw_hajek = "Coverage",
      coverage_shortfall = "Shortfall",
      abs_calibration_error = "Abs. cal. error",
      median_lpb = "Median LPB",
      effective_sample_size = "ESS",
      log10_max_raw_weight = "log10 max weight"
    ))

  p4 <- ggplot(raw_long, aes(x = method, y = value, fill = method)) +
    geom_boxplot(outlier.alpha = 0.10, width = 0.65) +
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = GBSG_COLORS) +
    labs(title = "GBSG real-data application: split-level metric distributions", x = "Method", y = NULL) +
    gbsg_theme_pub() +
    theme(legend.position = "none")

  gbsg_save_fig(p4, out_dir, "GBSG_Figure_4_metric_boxplots", width = 11, height = 7)
}

############################################################
## End of file
############################################################
