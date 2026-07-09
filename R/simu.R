## =====================================================================
## dgp.R -- Data-generating mechanisms (Section 3.2 of the paper)
##
## Event-time settings A (Weibull AFT, correctly specified), B (nonlinear
## misspecified), C (heteroscedastic). Censoring mechanisms C1 (PH, aligned
## with the Cox working model), C2 (nonlinear, misspecifies Cox), C3
## (near-positivity stress via a high-hazard subgroup). lambda0 is tuned per
## scenario to a target censoring proportion.
## =====================================================================

## Covariates: X1,X2 ~ N(0,1); X3,X4 ~ Bernoulli(0.5)
gen_covariates <- function(n) {
  data.frame(
    X1 = rnorm(n), X2 = rnorm(n),
    X3 = rbinom(n, 1, 0.5), X4 = rbinom(n, 1, 0.5)
  )
}

## log event time given covariates and event setting
gen_logT <- function(X, setting = c("A", "B", "C")) {
  setting <- match.arg(setting)
  n <- nrow(X)
  if (setting == "A") {
    mu  <- 0.5 * X$X1 - 0.5 * X$X2 + 0.8 * X$X3 - 0.6 * X$X4
    ## standard extreme-value (Gumbel-min) innovation: W = log(E), E ~ Exp(1)
    eps <- log(rexp(n))
    logT <- mu + 0.60 * eps
  } else if (setting == "B") {
    mu  <- 0.8 * sin(X$X1) + 0.5 * X$X2^2 - 0.7 * X$X3 * X$X4
    eps <- rnorm(n)
    logT <- mu + 0.60 * eps
  } else {
    mu    <- 0.5 * X$X1 - 0.4 * X$X2 + 0.6 * X$X3
    sigma <- 0.5 + 0.3 * abs(X$X1)
    eps   <- rnorm(n)
    logT  <- mu + sigma * eps
  }
  logT
}

## censoring linear predictor m(X): rate = lambda0 * exp(m(X))
gen_cens_lp <- function(X, mech = c("C1", "C2", "C3")) {
  mech <- match.arg(mech)
  if (mech == "C1") {
    0.5 * X$X1 - 0.5 * X$X2
  } else if (mech == "C2") {
    0.8 * sin(X$X1) + 0.5 * X$X2^2 - 0.7 * X$X3 * X$X4
  } else {
    0.5 * X$X1 - 0.5 * X$X2 + 1.2 * (X$X1 > 1.25)
  }
}

## Tune lambda0 so that the marginal censoring proportion ~= target_cens.
## P(censored | T, X) = 1 - exp(-lambda0 * exp(m) * T), solved on a pilot sample.
tune_lambda0 <- function(setting, mech, target_cens,
                         n_pilot = 2e5, seed = 1) {
  ## Tune lambda0 without disturbing the replicate-level RNG stream.
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  X  <- gen_covariates(n_pilot)
  Tt <- exp(gen_logT(X, setting))
  m  <- gen_cens_lp(X, mech)
  emt <- exp(m) * Tt
  f <- function(loglam) mean(1 - exp(-exp(loglam) * emt)) - target_cens
  ur <- uniroot(f, lower = -25, upper = 25)
  exp(ur$root)
}

## Full observed dataset (latent T and C retained for synthetic evaluation)
gen_data <- function(n, setting, mech, lambda0) {
  X  <- gen_covariates(n)
  Tt <- exp(gen_logT(X, setting))
  m  <- gen_cens_lp(X, mech)
  C  <- rexp(n, rate = lambda0 * exp(m))
  Y  <- pmin(Tt, C)
  Delta <- as.integer(Tt <= C)
  list(X = X, T = Tt, C = C, Y = Y, Delta = Delta, n = n)
}

## Convenience: subset a generated dataset by row index
subset_data <- function(d, idx) {
  list(X = d$X[idx, , drop = FALSE], T = d$T[idx], C = d$C[idx],
       Y = d$Y[idx], Delta = d$Delta[idx], n = length(idx))
}
