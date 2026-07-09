## =====================================================================
## metrics.R -- Performance metrics (Section 3.7-3.8)
## =====================================================================

## Per-replicate, per-method evaluation on the synthetic test set (true T known).
eval_method <- function(res, test, nu, alpha = 0.10, event_lp = NULL) {
  lpb <- res$lpb
  cov <- mean(test$T >= lpb)
  short <- max(0, (1 - alpha) - cov)
  abscal <- abs(cov - (1 - alpha))

  ## group-conditional and worst-slice coverage (Section 3.9)
  groups <- list(X3 = test$X$X3, X4 = test$X$X4)
  if (!is.null(event_lp)) {
    q <- cut(event_lp, breaks = quantile(event_lp, probs = seq(0, 1, .25)),
             include.lowest = TRUE)
    groups$LPq <- q
  }
  cov_g <- unlist(lapply(groups, function(g)
    tapply(test$T >= lpb, g, mean)))
  worst <- if (length(cov_g)) min(cov_g, na.rm = TRUE) else NA_real_

  data.frame(
    tau        = res$tau,
    coverage   = cov,
    shortfall  = short,
    abs_cal_err = abscal,
    mean_lpb   = mean(lpb),
    med_lpb    = median(lpb),
    mean_loglpb = mean(log(pmax(lpb, 1e-8))),
    worst_slice = worst,
    max_raw_wt  = res$diag$max_raw_wt,
    max_norm_wt = res$diag$max_norm_wt,
    ess         = res$diag$ess,
    wt_cv       = res$diag$cv,
    clipfrac    = res$diag$clipfrac
  )
}

## Aggregate a long per-replicate data.frame to method-level summaries with
## Monte Carlo SEs, PAC-success, and across-replicate variability of tau/coverage.
aggregate_metrics <- function(df, alpha = 0.10) {
  mcse <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  by_method <- split(df, df$method)
  do.call(rbind, lapply(names(by_method), function(m) {
    d <- by_method[[m]]
    data.frame(
      method       = m,
      R            = nrow(d),
      mean_cov     = mean(d$coverage),
      mcse_cov     = mcse(d$coverage),
      med_cov      = median(d$coverage),
      pac_success  = mean(d$coverage >= (1 - alpha)),     # PAC-success rate
      mean_short   = mean(d$shortfall),
      mean_abscal  = mean(d$abs_cal_err),
      med_lpb      = median(d$med_lpb),
      worst_slice  = mean(d$worst_slice, na.rm = TRUE),
      sd_tau       = sd(d$tau),                            # threshold variability
      sd_cov       = sd(d$coverage),                       # coverage variability
      mean_maxwt   = mean(d$max_raw_wt, na.rm = TRUE),
      mean_ess     = mean(d$ess, na.rm = TRUE),
      mean_clipfrac = mean(d$clipfrac, na.rm = TRUE),
      row.names = NULL
    )
  }))
}


## Scenario-level summaries for manuscript tables.
aggregate_metrics_by_scenario <- function(df, alpha = 0.10) {
  group_cols <- intersect(c("setting", "mech", "n", "cens", "backend", "method"), names(df))
  if (!length(group_cols)) return(aggregate_metrics(df, alpha))
  key <- do.call(interaction, c(df[group_cols], drop = TRUE, sep = "|"))
  pieces <- split(df, key)
  ans <- do.call(rbind, lapply(pieces, function(d) {
    out <- aggregate_metrics(d, alpha)
    for (cc in rev(setdiff(group_cols, "method"))) out[[cc]] <- d[[cc]][1]
    out
  }))
  rownames(ans) <- NULL
  front <- intersect(c("setting", "mech", "n", "cens", "backend", "method"), names(ans))
  ans[, c(front, setdiff(names(ans), front)), drop = FALSE]
}
