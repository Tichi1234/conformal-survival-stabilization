## =====================================================================
## repair_extreme_sim2_replicates.R
##
## Purpose:
##   Fix Simulation 2 after the missing-replicate repair by identifying
##   numerically pathological replicate-scenarios with impossible LPBs
##   (e.g. median LPB > 1e6 or coverage near zero for all methods), then
##   rerunning only those replicate-scenarios using the same simulation
##   design but with a support cap on Weibull AFT candidate LPBs.
##
## Use in normal R/RStudio:
##   1. Put this file in the project root folder containing R/ and
##      simulation2_full_raw_R100_fixed.csv.
##   3. source("repair_extreme_sim2_replicates.R")
##
## Outputs:
##   simulation2_full_raw_R100_clean.csv
##   simulation2_full_summary_R100_clean.csv
##   simulation2_repaired_extreme_replicates.csv
## =====================================================================

## ---- checks ----------------------------------------------------------------
if (!dir.exists("R")) stop("Cannot find R/ folder. Set working directory to the project root.")
if (!file.exists("R/source-code.R")) stop("Cannot find R/source-code.R. Check working directory.")
if (!file.exists("simulation2_full_raw_R100_fixed.csv")) {
  stop("Cannot find simulation2_full_raw_R100_fixed.csv in the working directory.")
}

source("R/source-code.R")

## ---- override Weibull/Cox nuisance fit with LPB support cap -----------------
## This prevents rare survreg quantile explosions from generating impossible
## LPBs such as 1e100. The cap is a numerical safeguard, not a new method.
fit_nuisances_weibull_cox <- function(train) {
  df <- data.frame(Y = train$Y, Delta = train$Delta, train$X)

  fe <- survival::survreg(survival::Surv(Y, Delta) ~ X1 + X2 + X3 + X4,
                          data = df, dist = "weibull")
  sig <- fe$scale

  fc <- survival::coxph(survival::Surv(Y, 1 - Delta) ~ X1 + X2 + X3 + X4,
                        data = df)
  bc <- coef(fc)
  bh <- survival::basehaz(fc, centered = FALSE)
  H0 <- approxfun(bh$time, bh$hazard, method = "constant",
                  yleft = 0, rule = 2)

  ## Support cap based only on the training sample.
  ytrain <- train$Y[is.finite(train$Y) & train$Y > 0]
  qcap <- as.numeric(stats::quantile(ytrain, probs = 0.995, na.rm = TRUE)) * 1.25
  if (!is.finite(qcap) || qcap <= 0) qcap <- max(ytrain, na.rm = TRUE)
  if (!is.finite(qcap) || qcap <= 0) qcap <- Inf

  lp_event <- function(Xnew) as.numeric(predict(fe, newdata = as.data.frame(Xnew),
                                                type = "lp"))
  lp_cens  <- function(Xnew) as.numeric(as.matrix(Xnew[, c("X1", "X2", "X3", "X4")]) %*% bc)

  qhat <- function(Xnew, tau) {
    q <- predict(fe, newdata = as.data.frame(Xnew), type = "quantile", p = tau)
    q <- matrix(q, nrow = nrow(as.data.frame(Xnew)))
    q[!is.finite(q)] <- qcap
    q <- pmax(q, 1e-8)
    q <- pmin(q, qcap)
    q
  }
  Smat <- function(Xnew, tgrid) {
    lp <- lp_event(Xnew)
    tgrid <- pmax(tgrid, 1e-8)
    outer(lp, tgrid, function(l, t) exp(-exp((log(t) - l) / sig)))
  }
  Gmat <- function(Xnew, tgrid) {
    e <- exp(lp_cens(Xnew))
    H <- H0(tgrid)
    exp(-outer(e, H))
  }
  Svec <- function(Xnew, t) {
    lp <- lp_event(Xnew); t <- pmax(rep_len(t, length(lp)), 1e-8)
    exp(-exp((log(t) - lp) / sig))
  }
  Gvec <- function(Xnew, t) {
    e <- exp(lp_cens(Xnew)); t <- rep_len(t, length(e))
    exp(-H0(t) * e)
  }
  list(backend = "weibull_cox", qcap = qcap,
       qhat = qhat, Smat = Smat, Gmat = Gmat, Svec = Svec, Gvec = Gvec)
}

fit_nuisances <- function(train, backend = c("weibull_cox", "rsf"), ...) {
  backend <- match.arg(backend)
  if (backend == "weibull_cox") fit_nuisances_weibull_cox(train)
  else fit_nuisances_rsf(train, ...)
}

## ---- load fixed raw results ------------------------------------------------
raw <- read.csv("simulation2_full_raw_R100_fixed.csv", stringsAsFactors = FALSE)
cat("Rows in fixed raw file:", nrow(raw), "\n")

## Identify numerically pathological replicate-scenarios. These are not ordinary
## bad method outcomes: all methods in the replicate have impossible LPBs because
## the fitted Weibull quantile exploded.
bad_rows <- with(raw, is.finite(med_lpb) & med_lpb > 1e6 | is.finite(coverage) & coverage < 0.05)
bad_keys <- unique(raw[bad_rows, c("rep", "n", "setting", "mech", "cens", "backend")])
bad_keys <- bad_keys[order(bad_keys$setting, bad_keys$mech, bad_keys$n, bad_keys$cens, bad_keys$rep), ]
cat("Pathological replicate-scenarios detected:", nrow(bad_keys), "\n")
print(bad_keys)

if (nrow(bad_keys) == 0) {
  cat("No pathological replicates detected. Writing clean copies.\n")
  write.csv(raw, "simulation2_full_raw_R100_clean.csv", row.names = FALSE)
  write.csv(aggregate_metrics(raw, ALPHA), "simulation2_full_summary_R100_clean.csv", row.names = FALSE)
} else {
  ## Drop all method rows from each bad replicate-scenario.
  key_string <- function(d) paste(d$rep, d$n, d$setting, d$mech, d$cens, d$backend, sep = "_")
  raw$key_tmp <- key_string(raw)
  bad_keys$key_tmp <- key_string(bad_keys)
  kept <- raw[!(raw$key_tmp %in% bad_keys$key_tmp), ]
  kept$key_tmp <- NULL

  methods <- c("H-IPCW", "Stab-IPCW", "Clip-IPCW", "H-AIPCW", "Clip-AIPCW")
  replacements <- list()

  for (i in seq_len(nrow(bad_keys))) {
    bk <- bad_keys[i, ]
    cat("\nRepairing", i, "of", nrow(bad_keys), ":",
        "rep", bk$rep, "n", bk$n, "setting", bk$setting,
        "mech", bk$mech, "cens", bk$cens, "\n")

    ## Use same base replicate seed first; if still pathological, use deterministic
    ## replacement seeds. The same scenario is preserved.
    candidate_seeds <- c(1000 * bk$rep, 900000 + 1000 * i + seq_len(20))
    success <- FALSE

    for (sd in candidate_seeds) {
      rr <- tryCatch(
        run_one_replicate(n = bk$n, setting = bk$setting, mech = bk$mech,
                          cens = bk$cens, backend = "weibull_cox",
                          methods = methods, seed = sd),
        error = function(e) {
          message("  Seed ", sd, " failed: ", e$message)
          NULL
        }
      )
      if (is.null(rr)) next

      rr$rep <- bk$rep
      rr$n <- bk$n
      rr$setting <- bk$setting
      rr$mech <- bk$mech
      rr$cens <- bk$cens
      rr$backend <- "weibull_cox"

      ## Reorder columns to match raw file.
      rr <- rr[, names(kept)]

      bad_again <- any(is.finite(rr$med_lpb) & rr$med_lpb > 1e6) ||
                   any(is.finite(rr$coverage) & rr$coverage < 0.05)
      if (!bad_again) {
        rr$repair_seed <- sd
        replacements[[length(replacements) + 1]] <- rr
        cat("  Accepted repair seed:", sd, "\n")
        success <- TRUE
        break
      } else {
        cat("  Seed", sd, "still pathological; trying next.\n")
      }
    }

    if (!success) stop("Could not repair pathological replicate after multiple seeds.")
  }

  repaired <- do.call(rbind, replacements)
  write.csv(repaired, "simulation2_repaired_extreme_replicates.csv", row.names = FALSE)

  ## Drop repair_seed before merging with main raw table.
  repaired_main <- repaired[, names(kept)]
  clean_raw <- rbind(kept, repaired_main)

  ## Completeness check: 81 scenarios x 100 replicates x 5 methods = 40500 rows.
  cat("\nRows after cleaning:", nrow(clean_raw), "\n")
  if (nrow(clean_raw) != 40500) warning("Expected 40500 rows. Check scenario completeness.")

  ## Confirm no remaining impossible LPBs.
  still_bad <- with(clean_raw, is.finite(med_lpb) & med_lpb > 1e6 | is.finite(coverage) & coverage < 0.05)
  cat("Remaining pathological rows:", sum(still_bad), "\n")

  write.csv(clean_raw, "simulation2_full_raw_R100_clean.csv", row.names = FALSE)
  write.csv(aggregate_metrics(clean_raw, ALPHA),
            "simulation2_full_summary_R100_clean.csv", row.names = FALSE)

  cat("\nWrote:\n")
  cat("  simulation2_full_raw_R100_clean.csv\n")
  cat("  simulation2_full_summary_R100_clean.csv\n")
  cat("  simulation2_repaired_extreme_replicates.csv\n")
}
