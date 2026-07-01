## =====================================================================
## scripts/06_clip-sensitivity.R
##
## Sensitivity of the clipped methods to the clipping percentile q (the choice
## promised in Section 2.5 of the manuscript). For each q the SAME calibration
## code path is used (the library reads getOption("herg_clip_q")), with data,
## splits, seeds and nuisances held fixed across q so any difference is due to
## the clipping level alone. q = 1.00 recovers the unclipped Hajek baseline.
##
## Output:
##   results/clip_sensitivity/clip_sensitivity_raw.csv
##   results/clip_sensitivity/clip_sensitivity_summary.csv
##   figures/clip_sensitivity.png
##
## Run from the repository root:
##   source("scripts/06_clip-sensitivity.R")
## =====================================================================

source("R/source-code.R")
suppressMessages(library(ggplot2))

Q_GRID  <- c(1.00, 0.95, 0.90, 0.85, 0.80)   # 1.00 = no clipping (Hajek)
METHODS <- c("Clip-IPCW", "Clip-AIPCW")
R_SENS  <- 100                                # replicates per cell

## Cells: the Simulation 1 regime plus one hard stress cell, to show the knee
## holds in both the favourable and the difficult regime.
CELLS <- rbind(
  expand.grid(setting = "A", mech = "C1", n = c(600, 1200), cens = c(0.20, 0.40, 0.60),
              stringsAsFactors = FALSE),
  data.frame(setting = "C", mech = "C3", n = 600, cens = 0.80, stringsAsFactors = FALSE)
)

dir.create("results/clip_sensitivity", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

## One replicate: generate data once, fit nuisances once, sweep q on the same fit.
one_rep <- function(n, setting, mech, cens, seed) {
  lam <- get_lambda0(setting, mech, cens)
  set.seed(seed)
  obs   <- gen_data(n, setting, mech, lam)
  itr   <- sample.int(n, floor(0.6 * n))
  train <- subset_data(obs, itr)
  cal   <- subset_data(obs, setdiff(seq_len(n), itr))
  test  <- gen_data(NTEST, setting, mech, lam)
  nu    <- fit_nuisances(train, backend = "weibull_cox")
  ev_lp <- tryCatch(as.numeric(nu$qhat(test$X, 0.5)), error = function(e) NULL)

  do.call(rbind, lapply(Q_GRID, function(q) {
    old <- options(herg_clip_q = q); on.exit(options(old), add = TRUE)
    res <- run_internal_methods(nu, cal, TAUGRID, test$X, alpha = ALPHA,
                                methods = METHODS, aug_grid_size = 50)
    do.call(rbind, lapply(names(res), function(m)
      cbind(q = q, method = m,
            eval_method(res[[m]], test, nu, alpha = ALPHA, event_lp = ev_lp))))
  }))
}

## Driver
raw <- list()
for (i in seq_len(nrow(CELLS))) {
  cl <- CELLS[i, ]
  message(sprintf("[cell %d/%d] setting=%s mech=%s n=%d cens=%.2f",
                  i, nrow(CELLS), cl$setting, cl$mech, cl$n, cl$cens))
  for (r in seq_len(R_SENS)) {
    row <- tryCatch(
      cbind(rep = r, setting = cl$setting, mech = cl$mech, n = cl$n, cens = cl$cens,
            one_rep(cl$n, cl$setting, cl$mech, cl$cens, seed = 1000 * r)),
      error = function(e) { message("  rep ", r, " failed: ", conditionMessage(e)); NULL })
    if (!is.null(row)) raw[[length(raw) + 1]] <- row
  }
}
raw <- do.call(rbind, raw)

## Summary by q and method
agg <- function(d) data.frame(
  coverage    = mean(d$coverage),
  pac_success = mean(d$coverage >= 1 - ALPHA),
  abs_cal_err = mean(d$abs_cal_err),
  med_lpb     = mean(d$med_lpb),
  ess         = mean(d$ess),
  max_wt      = mean(d$max_raw_wt),
  clipfrac    = mean(d$clipfrac, na.rm = TRUE)
)
keys <- interaction(raw$q, raw$method, drop = TRUE, sep = "|")
summary_tbl <- do.call(rbind, lapply(split(raw, keys), function(d)
  cbind(q = d$q[1], method = d$method[1], agg(d))))
rownames(summary_tbl) <- NULL
summary_tbl <- summary_tbl[order(summary_tbl$method, -summary_tbl$q), ]

write.csv(raw, "results/clip_sensitivity/clip_sensitivity_raw.csv", row.names = FALSE)
write.csv(summary_tbl, "results/clip_sensitivity/clip_sensitivity_summary.csv", row.names = FALSE)

## Figure: the trade-off vs q (PAC-success and median LPB rescaled onto one panel)
long <- rbind(
  data.frame(q = summary_tbl$q, method = summary_tbl$method,
             metric = "PAC success", value = summary_tbl$pac_success),
  data.frame(q = summary_tbl$q, method = summary_tbl$method,
             metric = "Median LPB", value = summary_tbl$med_lpb),
  data.frame(q = summary_tbl$q, method = summary_tbl$method,
             metric = "log10 max weight", value = log10(pmax(summary_tbl$max_wt, 1))))
p <- ggplot(long, aes(x = q, y = value, colour = method)) +
  geom_line() + geom_point() +
  facet_wrap(~metric, scales = "free_y") +
  scale_x_reverse(breaks = Q_GRID) +
  labs(title = "Sensitivity to the clipping percentile q (1.00 = no clipping)",
       x = "Clipping percentile q", y = NULL, colour = NULL) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "bottom")
ggsave("figures/clip_sensitivity.png", p, width = 9, height = 3.8, dpi = 200)

message("Done. Wrote clip-sensitivity CSVs and figures/clip_sensitivity.png")
print(summary_tbl, row.names = FALSE, digits = 3)
