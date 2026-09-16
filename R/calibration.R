## =====================================================================
## calibration.R -- Calibration methods (Section 2.3-3.2)
##All methods use the centered coverage moment
##phi(tau)=theta(tau)-(1-alpha)
##and the common monotone_envelope selection rule
##     tau_hat = sup{ tau in T : inf_{tau' <= tau} phi(tau') >= 0 }.
##IMPORTANT
##The clipping threshold c is selected using the TRAINING fold and then 
##passed into run_internal_methods(). It is NOT selected using the 
##calibration fold.
## =====================================================================

## ---- common monotone-envelope selection rule ------------------------
##======================================================================
select_tau <- function(phi, tau_grid) {
  ord <- order(tau_grid)
  phi_o <- phi[ord]; tau_o <- tau_grid[ord]
  run_inf <- cummin(phi_o)                        # inf_{tau' <= tau}
  ok <- which(run_inf >= 0)
## If no candidate satisfies the calibration constrain,
##return the most conservative candidate on the grid.
if (!length(ok)){
  tau_o[1]
}else{
  tau_o[max(ok)]
  }
}

##=====================================================================
## ---- Training-fold percentile clipping rule -------------------------
##======================================================================
select_clip_percentile <- function(nu, train,
    q = getOption("herg_clip_q", 0.90)) {
  
  ## Estimated censoring survival at observed training times
  GY_train <- pmax(nu$Gvec(train$X, train$Y),1e-8
  )
  
  unc <- train$Delta == 1
  
  if (!any(unc)) {
    stop("No uncensored observations in the training fold.")
  }
  
  ## Positive inverse-censoring factors among uncensored subjects
  pos_icw <- 1 / GY_train[unc]
  
  cclip <- as.numeric(
    quantile(
      pos_icw,
      probs = q,
      names = FALSE,
      na.rm = TRUE
    )
  )
  
  list(
    cclip = cclip,
    q = q
  )
}

##===========================================================================
##---Alternative:bound-minimizing clipping rule -----------------------------
##===========================================================================
select_clip_bound <- function(nu, train, tau_grid, n_cal, delta = 0.05, q_grid = c(
      0.50, 0.60, 0.70, 0.80, 0.85, 0.90, 0.95, 0.975, 0.99, 1.00)) {
  
  GY_train <- pmax(
    nu$Gvec(train$X, train$Y), 1e-8
  )
  
  unc <- train$Delta == 1
  
  if (!any(unc)) {
    stop("No uncensored observations in the training fold.")
  }
  
  pos_icw <- 1 / GY_train[unc]
  
  ## Candidate clipping constants generated from training-fold
  ## positive-weight quantiles
  c_grid <- unique(
    as.numeric(
      quantile(
        pos_icw,
        probs = q_grid,
        names = FALSE,
        na.rm = TRUE
      )
    )
  )
  
  tab <- do.call(
    rbind,
    lapply(c_grid, function(cval) {
      
      ## IMPORTANT:
      ## nu_hat is averaged over ALL training observations.
      ## Censored subjects contribute zero through Delta.
      w_clip <- train$Delta *
        pmin(
          1 / GY_train,
          cval
        )
      
      nu_hat <- mean(
        w_clip,
        na.rm = TRUE
      )
      
      if (!is.finite(nu_hat) || nu_hat <= 0) {
        
        return(
          data.frame(
            c = cval,
            nu_hat = nu_hat,
            b_hat = NA_real_,
            eps_hat = Inf,
            criterion = Inf
          )
        )
      }
      
      ## Training-fold proxy for clipping discrepancy
      b_hat <- 1 - nu_hat
      
      ## Empirical analogue of the stochastic term
      eps_hat <-
        (4 * cval / nu_hat) *
        sqrt(
          log(
            2 * (length(tau_grid) + 1) / delta
          ) /
            (2 * n_cal)
        )
      
      data.frame(
        c = cval,
        nu_hat = nu_hat,
        b_hat = b_hat,
        eps_hat = eps_hat,
        criterion = b_hat + eps_hat
      )
    })
  )
  
  valid <- which(
    is.finite(tab$criterion)
  )
  
  if (!length(valid)) {
    stop(
      "No valid candidate clipping level was found by the bound-minimizing rule."
    )
  }
  
  best_idx <- valid[
    which.min(tab$criterion[valid])
  ]
  
  best <- tab[
    best_idx,
    ,
    drop = FALSE
  ]
  
  list(
    cclip = best$c,
    nu_hat_train = best$nu_hat,
    b_hat = best$b_hat,
    eps_hat = best$eps_hat,
    criterion = best$criterion,
    diagnostics = tab
  )
}

##=============================================================================
## ----------- weight diagnostics ---------------------------------------------
##=============================================================================

.weight_diag <- function(w_raw, clipfrac = NA_real_, clip_level =NA_real_, 
                         nu_hat_c =NA_real_, runtime = NA_real_) {
  w_raw <- as.numeric(w_raw)
  
  sw <- sum(w_raw, na.rm = TRUE)
  sw2 <- sum(w_raw^2, na.rm =TRUE)
  ess  <- if (is.finite(sw2) && sw2 > 0) { sw^2 / sw2 }else {NA_real_}
  mean_w <- mean(w_raw, na.rm = TRUE)
  wn <- if (is.finite(mean_w) && mean_w > 0) {w_raw / mean_w} else { 
    rep( NA_real_, length(w_raw))}
  list(max_raw_wt  = if (length(w_raw) && any(is.finite(w_raw))) {
    max(w_raw, na.rm = TRUE)} else {NA_real_},
       max_norm_wt = if (any(is.finite(wn))) {
         max(wn,na.rm = TRUE)} else {NA_real_},
       ess = ess,
    cv =if (is.finite(mean_w) && mean_w > 0) {
        sd(w_raw, na.rm = TRUE) / mean_w} else {NA_real_},
    
    clipfrac = clipfrac,
    clip_level = clip_level,
    nu_hat_c = nu_hat_c,
    runtime = runtime
  )
}

## =====================================================================
## ---------Shared calibration ingredients------------------------------
## =====================================================================

## Returns:
##
## Ind[i,tau] = 1{Y_i >= qhat(tau | X_i)}
## raw IPCW weights
## conventionally stabilized weights
## externally supplied clipping threshold
## calibration-fold clipping fraction

calib_ingredients <- function(nu, cal, tau_grid, cclip = NULL) {
  
## Candidate lower predictive bounds
  Q <- nu$qhat(cal$X, tau_grid)
  
  ## Coverage-indicator matrix
  Ind <- (cal$Y >= Q) * 1
  
  ## Estimated censoring survival at observed times
  GY <- pmax(nu$Gvec(cal$X, cal$Y),1e-8)
  
  ## Raw IPCW weights
  w_raw <- cal$Delta / GY
  
  ## ---------------------------------------------------------------
  ## Marginal KM stabilization
  ## ---------------------------------------------------------------
  
  km <- survfit(
    Surv(
      cal$Y,
      1 - cal$Delta
    ) ~ 1)
  
  G0fun <- approxfun(
    km$time, km$surv,
    method = "constant",
    yleft = 1,
    rule = 2
  )
  
  w_stab <- cal$Delta *
    G0fun(cal$Y) / GY
  
  ## ---------------------------------------------------------------
  ## Calibration-fold clipping diagnostics
  ## cclip itself was selected on TRAINING data.
  ## ---------------------------------------------------------------
  
  if (!is.null(cclip)) {
    
    unc <- cal$Delta == 1
    
    if (any(unc)) {
      
      pos_icw <- 1 / GY[unc]
      
      clipfrac <- mean(
        pos_icw > cclip,
        na.rm = TRUE
      )
      
    } else {
      
      clipfrac <- NA_real_
    }
    
  } else {
    
    clipfrac <- NA_real_
  }
  
  list(Ind = Ind, GY = GY, w_raw = w_raw, w_stab = w_stab, cclip = cclip,
    clipfrac = clipfrac)
}

## =====================================================================
## ------------AIPCW augmentation term Pi(tau) -------------------------
## =====================================================================
## Backend-agnostic.
## eta_i(tau,u) = min(1-tau, S_i(u)) / S_i(u)
##
## clip = NULL: use raw 1/G factor
##
## clip = c: use min(1/G, c)

aipcw_augmentation <- function(nu, cal, tau_grid, alpha, clip = NULL,
    aug_grid_size = 50) {
  
  ctimes <- sort(
    unique(
      cal$Y[
        cal$Delta == 0]
    )
  )
  
  if (!length(ctimes)) {
    return(
      rep(0, length(tau_grid))
    )
  }
  
  if (length(ctimes) > aug_grid_size) {
    
    ctimes <- as.numeric(
      quantile(ctimes, probs = seq(0, 1, length.out = aug_grid_size),
        names = FALSE)
    )
  }
  
  U <- sort(unique(ctimes))
  
  Lg <- length(U)
  
  ## Event survival
  SU <- pmax(nu$Smat(cal$X, U), 1e-8)
  
  ## Censoring survival
  GU <- pmax(nu$Gmat(cal$X, U), 1e-8)
  
  ## Estimated cumulative censoring-hazard increments
  logGU <- log(GU)
  
  if (Lg == 1) {dLam <- matrix(-logGU[, 1], ncol = 1)}
    else {dLam <- cbind(
      -logGU[, 1],-(
        logGU[, -1, drop = FALSE] -
          logGU[, -Lg, drop = FALSE])
    )
  }
  
  ## Numerical safeguard
  dLam[!is.finite(dLam)] <- 0
  
  dLam[ dLam < 0] <- 0
  
  ## Raw or clipped inverse-censoring factor
  Omega <- if (is.null(clip)) {1 / GU}
   else {pmin( 1 / GU, clip)
  }
  
  ## At-risk process
  atrisk <- outer(cal$Y, U, ">=") * 1
  
  ## ---------------------------------------------------------------
  ## ------Jump term at each subject's observed censoring time -----
  ## ---------------------------------------------------------------
  
  SY <- pmax(nu$Svec(cal$X, cal$Y), 1e-8)
  
  GY <- pmax(nu$Gvec(cal$X, cal$Y), 1e-8)
  
  OmY <- if (is.null(clip)) {
    1 / GY}
   else {pmin(1 / GY, clip)}
  
  censored <- (cal$Delta == 0) * 1
  
  Pi <- numeric(length(tau_grid))
  
  for (j in seq_along(tau_grid)) {
    
    s <- 1 - tau_grid[j]
    
    Eta <- pmin(s, SU) / SU
    
    comp <- rowSums(atrisk *(Eta -(1 - alpha)) *Omega *dLam
    )
    
    etaY <- pmin(s, SY) / SY
    
    jump <- censored *(etaY - (1 - alpha)) * OmY
    
    Pi[j] <- mean(jump - comp)
  }
  
  Pi
}

## =====================================================================
## ----------Internal calibration methods ------------------------------
## =====================================================================

## Each method returns:
## list(
##   tau  = selected calibration level,
##   lpb  = lower predictive bound on X_test,
##   diag = weight diagnostics
## )

run_internal_methods <- function(nu, cal, tau_grid, X_test, alpha = 0.10,
    methods = c("Naive-Y", "CC", "HT-IPCW", "H-IPCW", "Stab-IPCW", 
                "Clip-IPCW", "H-AIPCW", "Clip-AIPCW"),
    aug_grid_size = 50,
    cclip = NULL) {
  
  ## ---------------------------------------------------------------
  ## Check whether a clipping threshold is required
  ## ---------------------------------------------------------------
  
  needs_clip <- any(
    c("Clip-IPCW", "Clip-AIPCW") %in% methods
  )
  
  if (
    needs_clip && is.null(cclip)
  ) {
    stop(
      paste0(
        "Clip-IPCW or Clip-AIPCW was requested, ",
        "but no clipping threshold cclip was supplied. ",
        "Select cclip from the training fold before calibration."
      )
    )
  }
  
  ## Shared ingredients
  ing <- calib_ingredients(nu = nu, cal = cal, tau_grid = tau_grid, 
                    cclip = cclip
  )
  
  Ind <- ing$Ind
  
  n <- cal$n
  
  one_minus_a <- 1 - alpha
  
  out <- list()
  
  
  ## ---------------------------------------------------------------
  ## Helper for constructing final method object
  ## ---------------------------------------------------------------
  
  finish <- function(phi, w_raw_diag, clipfrac = NA_real_, 
    clip_level = NA_real_, nu_hat_c = NA_real_, runtime = NA_real_) {
    
    th <- select_tau(phi, tau_grid
    )
    
    list(tau = th,
      
      lpb = as.numeric(
        nu$qhat(X_test, th)
      ),
      
      diag = .weight_diag(w_raw_diag, clipfrac = clipfrac,
        clip_level = clip_level, nu_hat_c = nu_hat_c,
        runtime = runtime)
    )
  }
  
  
  ## =================================================================
  ## Naive-Y
  ## =================================================================
  
  if ("Naive-Y" %in% methods) {
    
    phi <- colMeans(Ind) - one_minus_a
    
    out[["Naive-Y"]] <- finish(phi, rep(1, n))
  }
  
  
  ## =================================================================
  ## Complete case
  ## =================================================================
  
  if ("CC" %in% methods) {
    
    n_unc <- sum(cal$Delta)
    
    if (n_unc > 0) {
      
      theta_cc <- colSums(cal$Delta * Ind) /  n_unc
      
      phi <- theta_cc - one_minus_a
      
    } else {
      phi <- rep(-Inf, length(tau_grid))
    }
    
    out[["CC"]] <- finish(phi, cal$Delta)
  }
  
  
  ## =================================================================
  ## Horvitz-Thompson IPCW
  ## =================================================================
  
  if ("HT-IPCW" %in% methods) {
    
    phi <- colMeans(ing$w_raw * Ind) - one_minus_a
    
    out[["HT-IPCW"]] <- finish( phi, ing$w_raw)
  }
  
  
  ## =================================================================
  ## Hajek IPCW
  ## =================================================================
  
  if ("H-IPCW" %in% methods) {
    
    mean_w <- mean(ing$w_raw)
    
    if (
      is.finite(mean_w) && mean_w > 0) {
      
      wt <- ing$w_raw / mean_w
      
      phi <- colMeans(
        wt *
          (
            Ind -
              one_minus_a
          )
      )
      
    } else {
      
      phi <- rep(
        -Inf,
        length(tau_grid)
      )
    }
    
    out[["H-IPCW"]] <- finish(
      phi,
      ing$w_raw
    )
  }
  
  
  ## =================================================================
  ## Conventionally stabilized IPCW
  ## =================================================================
  
  if ("Stab-IPCW" %in% methods) {
    
    mean_w <- mean(
      ing$w_stab
    )
    
    if (
      is.finite(mean_w) &&
      mean_w > 0
    ) {
      
      wt <- ing$w_stab /
        mean_w
      
      phi <- colMeans(
        wt *
          (
            Ind -
              one_minus_a
          )
      )
      
    } else {
      
      phi <- rep(
        -Inf,
        length(tau_grid)
      )
    }
    
    out[["Stab-IPCW"]] <- finish(
      phi,
      ing$w_stab
    )
  }
  
  
  ## =================================================================
  ## Construct clipped raw weights once for both clipped methods
  ## =================================================================
  
  if (needs_clip) {
    
    w_clip <- cal$Delta *
      pmin(1 / ing$GY, ing$cclip)
    
    ## Realized calibration-fold mean clipped weight
    nu_hat_c_cal <- mean(
      w_clip, na.rm = TRUE)
  }
  
  
  ## =================================================================
  ## Clip-IPCW
  ## =================================================================
  
  if ("Clip-IPCW" %in% methods) {
    
    mean_w <- mean(
      w_clip
    )
    
    if (
      is.finite(mean_w) &&
      mean_w > 0
    ) {
      
      wt <- w_clip /
        mean_w
      
      phi <- colMeans(
        wt *
          (
            Ind -
              one_minus_a
          )
      )
      
    } else {
      
      phi <- rep(
        -Inf,
        length(tau_grid)
      )
    }
    
    out[["Clip-IPCW"]] <- finish(phi, w_clip, clipfrac = ing$clipfrac,
      clip_level = ing$cclip, nu_hat_c = nu_hat_c_cal)
  }
  
  
  ## =================================================================
  ## Hajek-AIPCW
  ## =================================================================
  
  if ("H-AIPCW" %in% methods) {
    
    mean_w <- mean(
      ing$w_raw
    )
    
    if (
      is.finite(mean_w) &&
      mean_w > 0
    ) {
      
      wt <- ing$w_raw /
        mean_w
      
      phiH <- colMeans(
        wt *
          (
            Ind -
              one_minus_a
          )
      )
      
      Pi <- aipcw_augmentation(nu = nu, cal = cal, tau_grid = tau_grid,
        alpha = alpha, clip = NULL, aug_grid_size = aug_grid_size)
      
      phi <- phiH + Pi
      
    } else {
      
      phi <- rep(
        -Inf,
        length(tau_grid)
      )
    }
    
    out[["H-AIPCW"]] <- finish(
      phi,
      ing$w_raw
    )
  }
  
  
  ## =================================================================
  ## Clip-AIPCW
  ## =================================================================
  
  if ("Clip-AIPCW" %in% methods) {
    
    mean_w <- mean(
      w_clip
    )
    
    if (
      is.finite(mean_w) &&
      mean_w > 0
    ) {
      
      wt <- w_clip /
        mean_w
      
      phiC <- colMeans(
        wt *
          (
            Ind -
              one_minus_a
          )
      )
      
      Pi <- aipcw_augmentation(nu = nu, cal = cal, tau_grid = tau_grid,
        alpha = alpha, clip = ing$cclip, aug_grid_size = aug_grid_size
      )
      
      phi <- phiC + Pi
      
    } else {
      
      phi <- rep(-Inf, length(tau_grid))
    }
    
    out[["Clip-AIPCW"]] <- finish(phi, w_clip, clipfrac = ing$clipfrac,
      clip_level = ing$cclip, nu_hat_c = nu_hat_c_cal)
  }
  out
}


## =====================================================================
## -------------------- External competitors ----------------------------
## =====================================================================
## These remain deliberately as stubs. Published external methods should
## continue to be called through the authors' released implementations.
## =====================================================================

fit_sesia_svetnik <- function(train, cal, X_test, alpha = 0.10) {
  
  message(
    paste0(
      "[competitor] Sesia-Svetnik-DR: stub -- ",
      "call authors' R implementation here.")
  )
  NULL
}


