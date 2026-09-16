## =====================================================================
## nuisance.R -- Nuisance estimation 
##

## Two backends, both returning the SAME interface object `nu` with:
##   nu$qhat(Xnew, tau)   -> matrix [nrow(Xnew) x length(tau)]  candidate LPBs
##   nu$Smat(Xnew, tgrid) -> matrix [nrow(Xnew) x length(tgrid)] event survival
##   nu$Gmat(Xnew, tgrid) -> matrix [nrow(Xnew) x length(tgrid)] censoring survival
##   nu$Svec(Xnew, t)     -> vector, S(t_i | x_i)   (t length 1 or nrow)
##   nu$Gvec(Xnew, t)     -> vector, G(t_i | x_i)
##
## RSF is used ONLY as a flexible nuisance estimator, never as a calibration
## method. Switching `backend` swaps the nuisances under every method.
## =====================================================================

suppressPackageStartupMessages(library(survival))

.feat <- function(X) as.matrix(X[, c("X1", "X2", "X3", "X4")])

## ---- Safe Weibull survreg fit --------------------------------------------

.fit_survreg_safe <- function(
  formula,
  data,
  dist = "weibull",
  prediction_multiplier = 1e6
) {

  ## ------------------------------------------------------------
  ## Validate a survreg fit.
  ##
  ## A fit is rejected not only when coefficients/scale are
  ## non-finite, but also when fitted event times are numerically
  ## implausible relative to the observed time scale.
  ## ------------------------------------------------------------
  fit_is_valid <- function(fit) {

    if (is.null(fit)) {
      return(FALSE)
    }

    if (
      !all(is.finite(coef(fit))) ||
      !is.finite(fit$scale) ||
      fit$scale <= 0
    ) {
      return(FALSE)
    }

    pred <- tryCatch(
      predict(
        fit,
        newdata = data,
        type = "response"
      ),
      error = function(e) rep(NA_real_, nrow(data))
    )

    if (
      length(pred) == 0 ||
      any(!is.finite(pred)) ||
      any(pred <= 0)
    ) {
      return(FALSE)
    }

    ## Very generous data-adaptive numerical ceiling.
    ## This is intended only to catch catastrophic fits such as
    ## the 10^136 predictions observed in Simulation 1, replicate 92.
    y_ref <- max(
      data$Y[is.finite(data$Y) & data$Y > 0],
      na.rm = TRUE
    )

    if (!is.finite(y_ref) || y_ref <= 0) {
      y_ref <- 1
    }

    if (max(pred, na.rm = TRUE) > prediction_multiplier * y_ref) {
      return(FALSE)
    }

    TRUE
  }


  ## ------------------------------------------------------------
  ## Attempt 1
  ## ------------------------------------------------------------
  warn1 <- FALSE

  fit1 <- withCallingHandlers(
    survreg(
      formula,
      data = data,
      dist = dist,
      control = survreg.control(
        maxiter = 200
      )
    ),
    warning = function(w) {

      if (grepl(
        "Ran out of iterations|did not converge",
        conditionMessage(w),
        ignore.case = TRUE
      )) {

        warn1 <<- TRUE
        invokeRestart("muffleWarning")
      }
    }
  )

  good1 <- !warn1 && fit_is_valid(fit1)

  if (good1) {
    return(
      list(
        fit = fit1,
        retried = FALSE,
        fallback = FALSE
      )
    )
  }


  ## ------------------------------------------------------------
  ## Attempt 2
  ## ------------------------------------------------------------
  warn2 <- FALSE

  fit2 <- withCallingHandlers(
    survreg(
      formula,
      data = data,
      dist = dist,
      control = survreg.control(
        maxiter = 1000,
        rel.tolerance = 1e-7
      )
    ),
    warning = function(w) {

      if (grepl(
        "Ran out of iterations|did not converge",
        conditionMessage(w),
        ignore.case = TRUE
      )) {

        warn2 <<- TRUE
        invokeRestart("muffleWarning")
      }
    }
  )

  good2 <- !warn2 && fit_is_valid(fit2)

  if (good2) {
    return(
      list(
        fit = fit2,
        retried = TRUE,
        fallback = FALSE
      )
    )
  }


  ## ------------------------------------------------------------
  ## Attempt 3: empirical starting values
  ## ------------------------------------------------------------

  p <- length(
    attr(
      terms(formula),
      "term.labels"
    )
  )

  logY_event <- log(
    pmax(
      data$Y[data$Delta == 1],
      1e-8
    )
  )

  if (length(logY_event) >= 2) {

    init_intercept <- median(
      logY_event,
      na.rm = TRUE
    )

    init_scale <- sd(
      logY_event,
      na.rm = TRUE
    )

  } else {

    logY_all <- log(
      pmax(
        data$Y,
        1e-8
      )
    )

    init_intercept <- median(
      logY_all,
      na.rm = TRUE
    )

    init_scale <- sd(
      logY_all,
      na.rm = TRUE
    )
  }

  if (
    !is.finite(init_scale) ||
    init_scale <= 1e-4
  ) {
    init_scale <- 1
  }

  init_full <- c(
    init_intercept,
    rep(0, p),
    log(init_scale)
  )

  warn3 <- FALSE

  fit3 <- withCallingHandlers(
    survival::survreg(
      formula,
      data = data,
      dist = dist,
      init = init_full,
      control = survival::survreg.control(
        maxiter = 2000,
        rel.tolerance = 1e-7
      )
    ),
    warning = function(w) {

      if (grepl(
        "Ran out of iterations|did not converge",
        conditionMessage(w),
        ignore.case = TRUE
      )) {

        warn3 <<- TRUE
        invokeRestart("muffleWarning")
      }
    }
  )

  good3 <- !warn3 && fit_is_valid(fit3)

  if (!good3) {
    stop(
      paste(
        "Weibull survreg nuisance model failed",
        "numerical validity checks after all three fitting attempts."
      )
    )
  }

  list(
    fit = fit3,
    retried = TRUE,
    fallback = TRUE
  )
}

## ---- Weibull AFT event + Cox censoring -----------------------------------
fit_nuisances_weibull_cox <- function(train) {
  df <- data.frame(Y = train$Y, Delta = train$Delta, train$X)

## event-time Weibull AFT
fe_obj <- .fit_survreg_safe(Surv(Y, Delta) ~ X1 + X2 + X3 + X4,
  data = df, dist = "weibull")

fe <- fe_obj$fit
sig <- fe$scale

## censoring Cox (censoring = event of interest)
  fc <- coxph(Surv(Y, 1 - Delta) ~ X1 + X2 + X3 + X4, data = df)
  bc <- coef(fc)
  bh <- basehaz(fc, centered = FALSE)            # H0(t) at covariate = 0
  H0 <- approxfun(bh$time, bh$hazard, method = "constant",
                  yleft = 0, rule = 2)

  lp_event <- function(Xnew) as.numeric(predict(fe, newdata = as.data.frame(Xnew),
                                                type = "lp"))
  lp_cens  <- function(Xnew) as.numeric(.feat(Xnew) %*% bc)

  qhat <- function(Xnew, tau) {
    q <- predict(fe, newdata = as.data.frame(Xnew), type = "quantile", p = tau)
    matrix(q, nrow = nrow(as.data.frame(Xnew)), ncol = length(tau))
  }
  Smat <- function(Xnew, tgrid) {
    lp <- lp_event(Xnew)
    tgrid <- pmax(tgrid, 1e-8)
    outer(lp, tgrid, function(l, t) exp(-exp((log(t) - l) / sig)))
  }
  Gmat <- function(Xnew, tgrid) {
    e <- exp(lp_cens(Xnew))
    H <- H0(tgrid)                                # length(tgrid)
    exp(-outer(e, H))                             # nrow x length(tgrid)
  }
  Svec <- function(Xnew, t) {
    lp <- lp_event(Xnew); t <- pmax(rep_len(t, length(lp)), 1e-8)
    exp(-exp((log(t) - lp) / sig))
  }
  Gvec <- function(Xnew, t) {
    e <- exp(lp_cens(Xnew)); t <- rep_len(t, length(e))
    exp(-H0(t) * e)
  }
  list(backend = "weibull_cox",
       qhat = qhat, Smat = Smat, Gmat = Gmat, Svec = Svec, Gvec = Gvec, 
       survreg_retried = fe_obj$retried, survreg_fallback = fe_obj$fallback)
}

## ---- Random survival forest backend --------------------------------------
## step-survival evaluator: S_i(t) = surv at largest model time <= t
.step_surv <- function(surv, tm, tgrid) {
  idx <- findInterval(tgrid, tm)                  # 0 means t < min(tm) -> S=1
  out <- matrix(1, nrow = nrow(surv), ncol = length(tgrid))
  pos <- which(idx >= 1)
  if (length(pos)) out[, pos] <- surv[, idx[pos], drop = FALSE]
  out
}

fit_nuisances_rsf <- function(train, ntree = 500, nodesize = 15) {
  if (!requireNamespace("randomForestSRC", quietly = TRUE))
    stop("Install 'randomForestSRC' to use the RSF backend.")
  df <- data.frame(Y = train$Y, Delta = train$Delta,
                   Cind = 1 - train$Delta, train$X)

  fe <- randomForestSRC::rfsrc(Surv(Y, Delta) ~ X1 + X2 + X3 + X4,
                               data = df, ntree = ntree, nodesize = nodesize)
  fc <- randomForestSRC::rfsrc(Surv(Y, Cind) ~ X1 + X2 + X3 + X4,
                               data = df, ntree = ntree, nodesize = nodesize)
  te <- fe$time.interest
  tc <- fc$time.interest

  surv_e <- function(Xnew) predict(fe, newdata = as.data.frame(Xnew))$survival
  surv_c <- function(Xnew) predict(fc, newdata = as.data.frame(Xnew))$survival

  qhat <- function(Xnew, tau) {
    S <- surv_e(Xnew)                             # n x length(te)
    sapply(tau, function(tt) {
      thr <- 1 - tt
      apply(S, 1, function(s) {
        j <- which(s <= thr)[1]
        if (is.na(j)) max(te) else te[j]
      })
    })
  }
  Smat <- function(Xnew, tgrid) .step_surv(surv_e(Xnew), te, tgrid)
  Gmat <- function(Xnew, tgrid) .step_surv(surv_c(Xnew), tc, tgrid)
  Svec <- function(Xnew, t) {
    t <- rep_len(t, nrow(as.data.frame(Xnew)))
    S <- surv_e(Xnew); diag(.step_surv(S, te, t))
  }
  Gvec <- function(Xnew, t) {
    t <- rep_len(t, nrow(as.data.frame(Xnew)))
    S <- surv_c(Xnew); diag(.step_surv(S, tc, t))
  }
  list(backend = "rsf",
       qhat = qhat, Smat = Smat, Gmat = Gmat, Svec = Svec, Gvec = Gvec)
}

fit_nuisances <- function(train, backend = c("weibull_cox", "rsf"), ...) {
  backend <- match.arg(backend)
  if (backend == "weibull_cox") fit_nuisances_weibull_cox(train)
  else fit_nuisances_rsf(train, ...)
}
