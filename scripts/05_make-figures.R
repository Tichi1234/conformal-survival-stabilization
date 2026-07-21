## =====================================================================
## scripts/05_make-figures.R
##
## Builds the manuscript figures from the CSVs written by scripts 01-04 and 03b.
## Reads deterministic paths under results/ and writes PNGs to figures/.
##
##   Simulation 1 (Setting A/C1 ablation)  -> figures/simulation1_coverage.png
##                                            figures/simulation1_stability.png
##   Simulation 2 (stress test)            -> figures/simulation2_coverage.png
##                                            figures/simulation2_maxweight.png
##   Simulation 3 (external benchmark)     -> figures/simulation3_coverage.png
##   GBSG real data                        -> figures/gbsg_panels.png
##
## Each figure is built only if its input CSV exists, so the script is safe to
## run after any subset of the experiments.
##
## Run from the repository root:
##   source("scripts/05_make-figures.R")
## =====================================================================

suppressMessages(library(ggplot2))

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

## Manuscript method orderings -----------------------------------------------------
LEV_SIM1 <- c("Naive-Y", "CC", "HT-IPCW", "H-IPCW", "Stab-IPCW", "Clip-IPCW", "H-AIPCW", "Clip-AIPCW")
LEV_SIM2 <- c("H-IPCW", "Stab-IPCW", "Clip-IPCW", "H-AIPCW", "Clip-AIPCW")
LEV_SIM3 <- c("H-IPCW", "Clip-IPCW", "H-AIPCW", "Clip-AIPCW",
              "DR-COSARC (fixed)", "DR-COSARC (adaptive)", "KM de-censoring", "Oracle")

TEAL <- "#1F6F6F"

read_raw <- function(path) {
  if (!file.exists(path)) { message("  (skip) not found: ", path); return(NULL) }
  df <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) { message("  (skip) empty: ", path); return(NULL) }
  df
}

## Stack chosen columns into long (method, panel, value) for facetted boxplots.
stack_panels <- function(df, cols, labels) {
  do.call(rbind, Map(function(cc, lab) {
    if (!cc %in% names(df)) return(NULL)
    data.frame(method = df$method, panel = lab, value = df[[cc]], stringsAsFactors = FALSE)
  }, cols, labels))
}

box_facets <- function(df, method_levels, cols, labels, title, hline = NULL) {
  long <- stack_panels(df, cols, labels)
  long <- long[long$method %in% method_levels, ]
  long$method <- factor(long$method, levels = method_levels)
  long$panel <- factor(long$panel, levels = labels[labels %in% long$panel])

  p <- ggplot(long, aes(x = method, y = value)) +
    geom_boxplot(
      outlier.size = 0.4,
      fill = TEAL,
      alpha = 0.55,
      colour = "grey25"
    ) +
    facet_wrap(~panel, scales = "free_y") +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "grey92")
    )

  if (!is.null(hline)) {
    p <- p +
      geom_hline(
        data = data.frame(
          panel = factor(labels[1], levels = base::levels(long$panel))
        ),
        aes(yintercept = hline),
        linetype = "dashed",
        colour = "grey30"
      )
  }

  p
}
## ---- Simulation 1 --------------------------------------------------------------
message("Simulation 1 figures ...")
s1 <- read_raw("results/simulation1/simulation1_raw.csv")
if (!is.null(s1)) {


#Plotting guard for Median LPB:
#remove non-finite or numerically exploded LPB values from the plotting variable only.
s1$med_lpb[!is.finite(s1$med_lpb) | s1$med_lpb < 0 | s1$med_lpb > 1] <- NA_real_
	
  ggsave("figures/simulation1_coverage.png",
         box_facets(s1, LEV_SIM1, c("coverage", "med_lpb"),
                    c("Coverage", "Median LPB"),
                    "Simulation 1: coverage and informativeness", hline = 0.90),
         width = 9, height = 4.2, dpi = 200)
  ggsave("figures/simulation1_stability.png",
         box_facets(transform(s1, log10_maxwt = log10(pmax(max_raw_wt, 1))),
                    LEV_SIM1, c("ess", "log10_maxwt"),
                    c("Effective sample size", "log10 maximum weight"),
                    "Simulation 1: stability of censoring-adjusted calibration"),
         width = 9, height = 4.2, dpi = 200)
}

## ---- Simulation 2 --------------------------------------------------------------
message("Simulation 2 figures ...")
s2 <- read_raw("results/simulation2/simulation2_weibull_cox_raw.csv")
if (!is.null(s2)) {

s2$med_lpb[!is.finite(s2$med_lpb) | s2$med_lpb < 0 | s2$med_lpb > 1] <- NA_real_

  ggsave("figures/simulation2_coverage.png",
         box_facets(s2, LEV_SIM2, c("coverage", "med_lpb"),
                    c("Coverage", "Median LPB"),
                    "Simulation 2: coverage and informativeness under stress test", hline = 0.90),
         width = 8, height = 4.2, dpi = 200)
  ## Figure 2: mean maximum raw inverse-censoring weight (log scale)
  agg <- aggregate(max_raw_wt ~ method, data = s2, FUN = mean)
  agg <- agg[agg$method %in% LEV_SIM2, ]
  agg$method <- factor(agg$method, levels = LEV_SIM2)
  p2 <- ggplot(agg, aes(x = method, y = max_raw_wt, fill = method)) +
    geom_col(width = 0.7, colour = "grey25", show.legend = FALSE) +
    scale_y_log10() +
    labs(title = "Simulation 2: clipping controls extreme inverse-censoring weights",
         x = NULL, y = "Mean maximum raw inverse-censoring weight (log scale)") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", size = 12))
  ggsave("figures/simulation2_maxweight.png", p2, width = 7, height = 4.5, dpi = 200)
}

## ---- Simulation 3 (external benchmark) -----------------------------------------
message("Simulation 3 figure ...")
s3 <- read_raw("results/simulation3/simulation3_external_raw_R100.csv")
if (is.null(s3)) s3 <- read_raw("results/simulation3/simulation3_internal_only_raw.csv")
if (!is.null(s3)) {
  d <- s3[s3$method %in% LEV_SIM3, ]
  d$method <- factor(d$method, levels = LEV_SIM3)
  p3 <- ggplot(d, aes(x = method, y = coverage)) +
    geom_boxplot(outlier.size = 0.4, fill = TEAL, alpha = 0.55, colour = "grey25") +
    geom_hline(yintercept = 0.90, linetype = "dashed", colour = "grey30") +
    labs(title = "Simulation 3: coverage across external-benchmark scenario-replicates",
         x = NULL, y = "Empirical coverage") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", size = 12))
  ggsave("figures/simulation3_coverage.png", p3, width = 8, height = 4.5, dpi = 200)
}

## ---- GBSG real data ------------------------------------------------------------
message("GBSG figure ...")
g <- read_raw("results/gbsg/gbsg_realdata_raw_R100.csv")
if (!is.null(g) && "method" %in% names(g)) {
  g$log10_maxwt <- log10(pmax(if ("max_raw_wt" %in% names(g)) g$max_raw_wt else NA, 1))
  cols   <- intersect(c("coverage", "med_lpb", "ess", "log10_maxwt"), names(g))
  labels <- c(coverage = "Coverage", med_lpb = "Median LPB",
              ess = "Effective sample size", log10_maxwt = "log10 maximum weight")[cols]
  gbsg_long <- do.call(rbind, lapply(cols, function(v) {
  data.frame(
    method = g$method,
    metric = unname(labels[v]),
    value = g[[v]]
  )
}))

gbsg_long <- gbsg_long[is.finite(gbsg_long$value), ]
gbsg_long$method <- factor(gbsg_long$method, levels = LEV_SIM1)
gbsg_long$metric <- factor(gbsg_long$metric, levels = unname(labels))

p_gbsg <- ggplot(gbsg_long, aes(x = method, y = value)) +
  geom_boxplot(outlier.size = 0.4, fill = TEAL, alpha = 0.55, colour = "grey25") +
  geom_hline(data = data.frame(metric = "Coverage", yintercept = 0.90),
             aes(yintercept = yintercept),
             linetype = "dashed", colour = "grey30") +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(title = "GBSG real-data application: coverage, informativeness and stability",
       x = "Method", y = NULL) +
  theme_bw(base_size = 11) +  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 12))

ggsave("figures/gbsg_panels.png", p_gbsg,
       width = 9, height = 6, dpi = 200)
}

message("Done. Figures written to figures/")


