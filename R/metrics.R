## =====================================================================
## metrics.R -- Performance metrics
## =====================================================================

## Per-replicate evaluation and Monte Carlo summaries used in 
## Simulations 1-3
##
## additions:
## * 5th and 10th percentiles of replicate-level coverage
## * exact-nominal PAC success retained for comparison
## * tolerance_based PAC success
## * realised clipping level
## * realised mean clipped weight nu_hat_c
## * runtime field retained if later populated

## per-replicate, per-method evaluation
eval_method <- function(res, test, nu, alpha = 0.10, event_lp = NULL) {
  lpb <- as.numeric(res$lpb)
  
## marginal empirical test-set coverage
  cov <- mean(test$T >= lpb, na.rm = TRUE)
  short <- max(0, (1 - alpha) - cov)
  abscal <- abs(cov - (1 - alpha))

## group-conditional / worst-slice coverage
  groups <- list(X3 = test$X$X3, X4 = test$X$X4)
  if (!is.null(event_lp)) {
    
## guard against duplicate empirical quantiles
    brks <- unique(quantile(event_lp, probs = seq(0, 1, 0.25),
        na.rm = TRUE,
        names = FALSE))
    
    if (length(brks) >= 2) {
    q <- cut(event_lp, breaks = brks, include.lowest = TRUE)
    
    groups$LPq <- q}
  }
  cov_g <- unlist(lapply(groups, function(g) {
    tapply(test$T >= lpb, g, mean, na.rm = TRUE)}))
  worst <- if (length(cov_g) && any (is.finite(cov_g))) {min(cov_g, 
       na.rm = TRUE)} else {NA_real_}
  
##  extract diagnostics safely
  get_diag <- function(name) {
    
    if (!is.null(res$diag) &&
      !is.null(res$diag[[name]])
    ) {
      res$diag[[name]]
    } else {
      NA_real_}}
  
## return one row per replicate-method combination
  data.frame(
    tau        = res$tau,
    coverage   = cov,
    shortfall  = short,
    abs_cal_err = abscal,
    mean_lpb   = mean(lpb, na.rm =TRUE),
    med_lpb    = median(lpb, na.rm = TRUE),
    mean_loglpb = mean(log(pmax(lpb, 1e-8)), na.rm = TRUE),
    worst_slice = worst,
    max_raw_wt  = get_diag("max_raw_wt"),
    max_norm_wt = get_diag("max_norm_wt"),
    ess         = get_diag("ess"),
    wt_cv       = get_diag("cv"),
    clipfrac    = get_diag("clipfrac"),
## diagnostics
    clip_level    = get_diag("clip_level"),
    nu_hat_c      = get_diag("nu_hat_c"),
     runtime       = get_diag("runtime")
   )
}

## Aggregate replicate_level results to method_level summaries
aggregate_metrics <- function(df, alpha = 0.10, pac_tol = 0.02) {
  safe_mean <- function(x) {
    if (!length(x) ||
        all(is.na(x))
      ) {
      NA_real_
    } else {
      mean(x, na.rm = TRUE)
    }
  }
  safe_median <- function(x) {
    if (
      !length(x) ||
      all(is.na(x))
    ) {
      NA_real_
    } else {
      median(x, na.rm = TRUE)
    }
  }
  safe_sd <- function(x) {
    x <- x[is.finite(x)]
    
    if (length(x) < 2) {
      NA_real_
    } else {sd(x)}
  }
  mcse <- function(x) {
    
    x <- x[is.finite(x)]
    
    if (length(x) < 2) {
      NA_real_
    } else {sd(x) /
        sqrt(length(x))}
  }
  
  safe_quantile <- function(
    x,
    p) {
    
    x <- x[is.finite(x)]
    
    if (!length(x)) {return(NA_real_)}
    
    as.numeric(
      quantile(x, probs = p, names = FALSE, na.rm = TRUE, type = 7))}

## split by method
  by_method <- split(df, df$method)
  
## calculate summaries
  ans <- do.call(rbind, lapply(names(by_method), function(m) {
        
        d <- by_method[[m]]
        
        data.frame(method = m,
          
          R = nrow(d),
          
          ## Mean empirical coverage
          mean_cov = safe_mean(d$coverage),
          
          ## Monte Carlo standard error of mean coverage
          mcse_cov = mcse(d$coverage),
          
          ## Median replicate coverage
          med_cov = safe_median(d$coverage),
    
  ##lower-tail reliability summaries
          cov_q05 = safe_quantile(d$coverage, 0.05),
          
          cov_q10 = safe_quantile(d$coverage, 0.10),
          
  ## Original PAC-success definition
            pac_success_nominal = safe_mean(d$coverage >= (1 - alpha)),
  
 ## tolerance-based PAC-style success
          pac_success_tol = safe_mean(d$coverage >= (1 - alpha - pac_tol)),
          
          pac_tol = pac_tol,
          
         
          ## Calibration error
         
          mean_short = safe_mean(d$shortfall),
          
          mean_abscal = safe_mean(d$abs_cal_err),
          
          ## Informativeness
          med_lpb = safe_median(d$med_lpb),
          
         ## Worst-slice reliability
          
          worst_slice = safe_mean(d$worst_slice),
          
          ## Across-replicate variability
          
          sd_tau = safe_sd(d$tau),
          
          sd_cov = safe_sd(d$coverage),
          
          ## Weight diagnostics
          mean_maxwt = safe_mean(d$max_raw_wt),
          
          mean_ess = safe_mean(d$ess),
          
          mean_clipfrac = safe_mean(d$clipfrac),
          
          ##clipping diagnostics
          mean_clip_level = safe_mean(d$clip_level),
          
          mean_nu_hat_c = safe_mean(d$nu_hat_c),
          ## computation-time diagnostic
          
          mean_runtime = safe_mean(d$runtime),
          
          row.names = NULL)}))
  
  rownames(ans) <- NULL
  
  ans
}

## Scenario-level summaries for manuscript tables

aggregate_metrics_by_scenario <- function(df, alpha = 0.10,
    pac_tol = 0.02) {
  
  group_cols <- intersect(
    c("setting", "mech", "n", "cens", "backend", "method"),
    names(df))
  
  ## If scenario identifiers are unavailable,
  ## fall back to overall method summaries.
  if (!length(group_cols)) {
    
    return(aggregate_metrics(
        df, alpha = alpha, pac_tol = pac_tol))}
  
  key <- do.call(interaction,
    c(df[group_cols], drop = TRUE, sep = "|"))
  
  pieces <- split(df, key)
  
  ans <- do.call(rbind, lapply(pieces, function(d) {
        
        out <- aggregate_metrics(d, alpha = alpha, pac_tol = pac_tol)
        
        ## Add scenario identifiers back to summary row
        for (cc in rev(setdiff(group_cols, "method"))
        ) {out[[cc]] <- d[[cc]][1]}
        
        out}))
  rownames(ans) <- NULL
  
  ## Put scenario identifiers at beginning of table
  front <- intersect(
    c("setting", "mech", "n", "cens", "backend", "method"),
    names(ans))
  
  ans[, c(front, setdiff(names(ans), front)),
    drop = FALSE]}
