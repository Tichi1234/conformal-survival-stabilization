## =====================================================================
## sesia-adapter.R
##
## Adapter to Sesia & Svetnik (2025) DR-COSARC implementation.
## The authors' source code is kept unchanged under external/sesia/.
## =====================================================================

SESIA_DIR <- file.path(
  "external",
  "sesia",
  "conformal_survival-main",
  "code",
  "conf_surv"
)

if (!dir.exists(SESIA_DIR)) {
  stop(
    paste(
      "Sesia source directory not found:",
      SESIA_DIR
    )
  )
}

## Authors' utilities
source(file.path(SESIA_DIR, "utils_survival.R"))
source(file.path(SESIA_DIR, "utils_censoring.R"))
source(file.path(SESIA_DIR, "utils_decensoring.R"))
source(file.path(SESIA_DIR, "utils_conformal.R"))


## ---------------------------------------------------------------------
## Convert our simulation objects to the Sesia data format
## ---------------------------------------------------------------------

.to_sesia_observed <- function(dat) {

  data.frame(
    time   = as.numeric(dat$Y),
    status = as.integer(dat$Delta),
    X1     = as.numeric(dat$X$X1),
    X2     = as.numeric(dat$X$X2),
    X3     = as.numeric(dat$X$X3),
    X4     = as.numeric(dat$X$X4)
  )
}


.to_sesia_test <- function(X) {

  data.frame(
    X1 = as.numeric(X$X1),
    X2 = as.numeric(X$X2),
    X3 = as.numeric(X$X3),
    X4 = as.numeric(X$X4)
  )
}


## ---------------------------------------------------------------------
## Sesia-Svetnik adaptive DR-COSARC
## ---------------------------------------------------------------------

fit_sesia_svetnik <- function(
  train,
  cal,
  X_test,
  alpha = 0.10,
  finite_sample_correction = FALSE
) {

  data.train <- .to_sesia_observed(train)
  data.cal   <- .to_sesia_observed(cal)
  data.test  <- .to_sesia_test(X_test)


  ## Event model
  ##
  ## Use Weibull survreg so that the event-model family is aligned with
  ## our default simulation nuisance model.
  surv_model <- SurvregModelWrapper$new(
    dist = "weibull"
  )

  surv_model$fit(
    survival::Surv(time, status) ~ .,
    data = data.train
  )


  ## Censoring model
  ##
  ## Cox model, wrapped using the authors' original CensoringModel class.
  cens_base_model <- CoxphModelWrapper$new()

  cens_imputator <- CensoringModel$new(
    model = cens_base_model
  )

  cens_imputator$fit(
    data = data.train
  )


  ## Adaptive doubly robust COSARC
  lpb <- predict_drcosarc(
    data.test = data.test,
    surv_model = surv_model,
    cens_imputator = cens_imputator,
    data.cal = data.cal,
    alpha = alpha,
    cutoffs = "adaptive",
    finite_sample_correction = finite_sample_correction,
    doubly_robust = TRUE
  )

  lpb <- as.numeric(lpb)


  ## Basic numerical checks
  if (length(lpb) != nrow(data.test)) {

    stop(
      sprintf(
        "Sesia-Svetnik returned %d predictions; expected %d.",
        length(lpb),
        nrow(data.test)
      )
    )
  }

  if (any(!is.finite(lpb))) {

    stop(
      "Sesia-Svetnik returned non-finite lower predictive bounds."
    )
  }


  ## Return object compatible with our evaluation framework
  list(
    tau = NA_real_,

    lpb = lpb,

    diag = list(
      max_raw_wt = NA_real_,
      max_norm_wt = NA_real_,
      ess = NA_real_,
      cv = NA_real_,
      clipfrac = NA_real_,
      clip_level = NA_real_,
      nu_hat_c = NA_real_,
      runtime = NA_real_
    )
  )
}

