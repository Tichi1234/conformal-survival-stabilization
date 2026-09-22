## ============================================================
## Table 6: clipping-percentile sensitivity
## Uses existing replicate-level results; NO simulations rerun
## ============================================================

raw_file <- "results/clip_sensitivity/clip_sensitivity_raw.csv"
out_file <- "results/clip_sensitivity/Table6_clipping_sensitivity.csv"

dat <- read.csv(raw_file, stringsAsFactors = FALSE)

## Keep the three manuscript values:
## q = 1.00 : no clipping
## q = 0.90 : primary/default rule
## q = 0.80 : more aggressive clipping
dat <- dat[dat$q %in% c(1.00, 0.90, 0.80), ]

## Check required variables
required <- c(
  "q", "method", "coverage",
  "med_lpb", "ess",
  "max_raw_wt", "clipfrac"
)

missing_vars <- setdiff(required, names(dat))

if (length(missing_vars) > 0) {
  stop(
    "Missing variables in raw sensitivity file: ",
    paste(missing_vars, collapse = ", ")
  )
}

## Summarise replicate-level results
summarise_one <- function(d) {

  data.frame(
    Coverage = mean(d$coverage, na.rm = TRUE),

    Q05 = as.numeric(
      quantile(d$coverage, probs = 0.05, na.rm = TRUE)
    ),

    Q10 = as.numeric(
      quantile(d$coverage, probs = 0.10, na.rm = TRUE)
    ),

    PAC_090 = mean(
      d$coverage >= 0.90,
      na.rm = TRUE
    ),

    PAC_088 = mean(
      d$coverage >= 0.88,
      na.rm = TRUE
    ),

    Median_LPB = mean(
      d$med_lpb,
      na.rm = TRUE
    ),

    ESS = mean(
      d$ess,
      na.rm = TRUE
    ),

    Max_raw_wt = mean(
      d$max_raw_wt,
      na.rm = TRUE
    ),

    Clipped_frac = mean(
      d$clipfrac,
      na.rm = TRUE
    ),

    N = sum(
      is.finite(d$coverage)
    )
  )
}

groups <- split(
  dat,
  interaction(dat$q, dat$method, drop = TRUE)
)

tab6 <- do.call(
  rbind,
  lapply(
    groups,
    function(d) {
      cbind(
        q = d$q[1],
        Method = d$method[1],
        summarise_one(d)
      )
    }
  )
)

rownames(tab6) <- NULL

## Order methods and q values
method_order <- c(
  "Clip-AIPCW",
  "Clip-IPCW"
)

q_order <- c(
  1.00,
  0.90,
  0.80
)

tab6$Method <- factor(
  tab6$Method,
  levels = method_order
)

tab6$q <- factor(
  tab6$q,
  levels = q_order
)

tab6 <- tab6[
  order(tab6$Method, tab6$q),
]

## Convert q back to numeric-like display
tab6$q <- as.character(tab6$q)

## Round for manuscript display
tab6$Coverage     <- round(tab6$Coverage, 3)
tab6$Q05          <- round(tab6$Q05, 3)
tab6$Q10          <- round(tab6$Q10, 3)
tab6$PAC_090      <- round(tab6$PAC_090, 3)
tab6$PAC_088      <- round(tab6$PAC_088, 3)
tab6$Median_LPB   <- round(tab6$Median_LPB, 3)
tab6$ESS          <- round(tab6$ESS, 1)
tab6$Max_raw_wt   <- round(tab6$Max_raw_wt, 2)
tab6$Clipped_frac <- round(tab6$Clipped_frac, 3)

write.csv(
  tab6,
  out_file,
  row.names = FALSE
)

cat("\nRevised Table 6\n\n")
print(tab6, row.names = FALSE)

cat(
  "\nSaved to:",
  out_file,
  "\n"
)

