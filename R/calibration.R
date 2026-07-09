## =====================================================================
## calibration.R -- Calibration methods (Section 2.3-2.5)
##
## Sign convention (fixes the inconsistent selection rules in the old draft):
## every method forms a CENTERED coverage moment phi(tau) = theta(tau)-(1-alpha)
## and selects via the common monotone-envelope rule
##     tau_hat = sup{ tau in T : inf_{tau' <= tau} phi(tau') >= 0 }.
## =====================================================================

## ---- common selection rule ------------------------------------------------
select_tau <- function(phi, tau_grid) {
  ord <- order(tau_grid)
  phi_o <- phi[ord]; tau_o <- tau_grid[ord]
  run_inf <- cummin(phi_o)                        # inf_{tau' <= tau}
  ok <- which(run_inf >= 0)
  if (!length(ok)) tau_o[1] else tau_o[max(ok)]   # fallback = most conservative
}

## weight diagnostics (Eq. 26 ESS; max/normalized weights; CV; clip fraction)
.weight_diag <- function(w_raw, clipfrac = NA_real_) {
  sw <- sum(w_raw); sw2 <- sum(w_raw^2)
  ess  <- if (sw2 > 0) sw^2 / sw2 else NA_real_
  wn   <- if (sw > 0) w_raw / mean(w_raw) else w_raw
  list(max_raw_wt  = max(w_raw),
       max_norm_wt = if (sw > 0) max(wn) else NA_real_,
       ess         = ess,
       cv          = if (mean(w_raw) > 0) sd(w_raw) / mean(w_raw) else NA_real_,
       clipfrac    = clipfrac)
}

## ---- shared calibration ingredients --------------------------------------
## Returns the coverage-indicator matrix Ind[i, tau] = 1{Y_i >= qhat(tau|X_i)},
## raw inverse-censoring weights, marginal-KM stabilized weights, clip level.
calib_ingredients <- function(nu, cal, tau_grid) {
  Q   <- nu$qhat(cal$X, tau_grid)                 # n_cal x |T|
  Ind <- (cal$Y >= Q) * 1
  GY  <- pmax(nu$Gvec(cal$X, cal$Y), 1e-8)        # G(Y_i | X_i)
  w_raw <- cal$Delta / GY                         # raw ICW
  ## marginal KM of the censoring distribution G0(t) for stabilized weights
  km <- survfit(Surv(cal$Y, 1 - cal$Delta) ~ 1)
  G0fun <- approxfun(km$time, km$surv, method = "constant",
                     yleft = 1, rule = 2)
  w_stab <- cal$Delta * G0fun(cal$Y) / GY
  ## adaptive 90th-percentile clip among uncensored
  unc <- cal$Delta == 1
  ## clip percentile is configurable via option (default 0.90); 1.0 => no clipping
  .q <- getOption("herg_clip_q", 0.90)
  cclip <- as.numeric(quantile(1 / GY[unc], .q, names = FALSE))
  clipfrac <- mean((1 / GY[unc]) > cclip)
  list(Ind = Ind, GY = GY, w_raw = w_raw, w_stab = w_stab,
       cclip = cclip, clipfrac = clipfrac)
}

## ---- AIPCW augmentation term Pi(tau) (raw or clipped) --------------------
## Backend-agnostic: hazard increments from -d log Ghat; eta uses the identity
## eta_i(tau,u) = min(1-tau, S_i(u)) / S_i(u).  clip = NULL -> raw 1/G factor.
aipcw_augmentation <- function(nu, cal, tau_grid, alpha,
                               clip = NULL, aug_grid_size = 50) {
  ctimes <- sort(unique(cal$Y[cal$Delta == 0]))
  if (!length(ctimes)) return(rep(0, length(tau_grid)))
  if (length(ctimes) > aug_grid_size)
    ctimes <- as.numeric(quantile(ctimes,
                probs = seq(0, 1, length.out = aug_grid_size), names = FALSE))
  U  <- sort(unique(ctimes)); Lg <- length(U)
  ncal <- cal$n

  SU <- pmax(nu$Smat(cal$X, U), 1e-8)             # n_cal x Lg
  GU <- pmax(nu$Gmat(cal$X, U), 1e-8)
  logGU <- log(GU)
  dLam <- cbind(-(logGU[, 1] - 0),
                -(logGU[, -1, drop = FALSE] - logGU[, -Lg, drop = FALSE]))
  dLam[dLam < 0] <- 0                             # hazard increments >= 0
  Omega <- if (is.null(clip)) 1 / GU else pmin(1 / GU, clip)
  atrisk <- outer(cal$Y, U, ">=") * 1             # Y_i >= u_k

  ## jump term at own censoring time
  SY  <- pmax(nu$Svec(cal$X, cal$Y), 1e-8)
  GY  <- pmax(nu$Gvec(cal$X, cal$Y), 1e-8)
  OmY <- if (is.null(clip)) 1 / GY else pmin(1 / GY, clip)
  censored <- (cal$Delta == 0) * 1

  Pi <- numeric(length(tau_grid))
  for (j in seq_along(tau_grid)) {
    s <- 1 - tau_grid[j]
    Eta  <- pmin(s, SU) / SU                      # n_cal x Lg
    comp <- rowSums(atrisk * (Eta - (1 - alpha)) * Omega * dLam)
    etaY <- pmin(s, SY) / SY
    jump <- censored * (etaY - (1 - alpha)) * OmY
    Pi[j] <- mean(jump - comp)
  }
  Pi
}

## ---- the eight internal methods -------------------------------------------
## Each returns list(tau, lpb (on X_test), diag).
run_internal_methods <- function(nu, cal, tau_grid, X_test, alpha = 0.10,
                                  methods = c("Naive-Y","CC","HT-IPCW","H-IPCW",
                                              "Stab-IPCW","Clip-IPCW",
                                              "H-AIPCW","Clip-AIPCW"),
                                  aug_grid_size = 50) {
  ing <- calib_ingredients(nu, cal, tau_grid)
  Ind <- ing$Ind; n <- cal$n
  one_minus_a <- 1 - alpha
  out <- list()

  finish <- function(phi, w_raw_diag, clipfrac = NA_real_) {
    th <- select_tau(phi, tau_grid)
    list(tau = th,
         lpb = as.numeric(nu$qhat(X_test, th)),
         diag = .weight_diag(w_raw_diag, clipfrac))
  }

  if ("Naive-Y" %in% methods) {
    phi <- colMeans(Ind) - one_minus_a
    out[["Naive-Y"]] <- finish(phi, rep(1, n))
  }
  if ("CC" %in% methods) {
    phi <- colSums(cal$Delta * Ind) / sum(cal$Delta) - one_minus_a
    out[["CC"]] <- finish(phi, cal$Delta)
  }
  if ("HT-IPCW" %in% methods) {
    phi <- colMeans(ing$w_raw * Ind) - one_minus_a    # unnormalized HT
    out[["HT-IPCW"]] <- finish(phi, ing$w_raw)
  }
  if ("H-IPCW" %in% methods) {
    wt <- ing$w_raw / mean(ing$w_raw)
    phi <- colMeans(wt * (Ind - one_minus_a))
    out[["H-IPCW"]] <- finish(phi, ing$w_raw)
  }
  if ("Stab-IPCW" %in% methods) {
    wt <- ing$w_stab / mean(ing$w_stab)
    phi <- colMeans(wt * (Ind - one_minus_a))
    out[["Stab-IPCW"]] <- finish(phi, ing$w_stab)
  }
  ## clipped raw weights shared by Clip-IPCW / Clip-AIPCW
  w_clip <- cal$Delta * pmin(1 / ing$GY, ing$cclip)
  if ("Clip-IPCW" %in% methods) {
    wt <- w_clip / mean(w_clip)
    phi <- colMeans(wt * (Ind - one_minus_a))
    out[["Clip-IPCW"]] <- finish(phi, w_clip, ing$clipfrac)
  }
  if ("H-AIPCW" %in% methods) {
    wt <- ing$w_raw / mean(ing$w_raw)
    phiH <- colMeans(wt * (Ind - one_minus_a))
    Pi <- aipcw_augmentation(nu, cal, tau_grid, alpha,
                             clip = NULL, aug_grid_size = aug_grid_size)
    out[["H-AIPCW"]] <- finish(phiH + Pi, ing$w_raw)
  }
  if ("Clip-AIPCW" %in% methods) {
    wt <- w_clip / mean(w_clip)
    phiC <- colMeans(wt * (Ind - one_minus_a))
    Pi <- aipcw_augmentation(nu, cal, tau_grid, alpha,
                             clip = ing$cclip, aug_grid_size = aug_grid_size)
    out[["Clip-AIPCW"]] <- finish(phiC + Pi, w_clip, ing$clipfrac)
  }
  out
}

## =====================================================================
## External competitors (Section 2.8 / Simulation 3) -- STUBS.
##
## These are intentionally NOT reimplemented. Reproducing a published method
## from memory risks an unfair comparison. Plug in the authors' released code:
##   Davidov et al. (2025, ICLR)  : Python (Romano group). Run on exported
##                                   splits, write an LPB vector back as CSV.
##   Sesia & Svetnik (2025, ICML) : R, github.com/msesia/conformal_survival.
##
## Each stub must return list(tau = NA, lpb = <numeric vector on X_test>).
## Until wired up they return NULL and the driver skips them.
## =====================================================================
fit_davidov <- function(train, cal, X_test, alpha = 0.10,
                        variant = c("focused", "fused")) {
  variant <- match.arg(variant)
  message("[competitor] Davidov-", variant,
          ": stub -- call authors' Python implementation here.")
  NULL
}
fit_sesia_svetnik <- function(train, cal, X_test, alpha = 0.10) {
  message("[competitor] Sesia-Svetnik-DR: stub -- call authors' R code here.")
  NULL
}
