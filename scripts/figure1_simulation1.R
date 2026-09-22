## ============================================================
## Figure 1: Simulation 1 internal ablation
## Panels:
## A = empirical coverage
## B = median lower predictive bound
## C = effective sample size
## D = maximum raw inverse-censoring weight (log scale)
## ============================================================

library(ggplot2)
library(dplyr)
library(patchwork)
source("R/figure_palette.R")

## ------------------------------------------------------------
## 1. Load corrected Simulation 1 data
## ------------------------------------------------------------


sim1 <- read.csv(
  "results/simulation1_pass3_R100/simulation1_raw_corrected.csv"
)

## Basic checks
dim(sim1)
table(sim1$method)

## ------------------------------------------------------------
## 2. Set method order
## ------------------------------------------------------------

method_order <- c(
  "Naive-Y",
  "CC",
  "HT-IPCW",
  "H-IPCW",
  "Stab-IPCW",
  "Clip-IPCW",
  "H-AIPCW",
  "Clip-AIPCW"
)

sim1$method <- factor(
  sim1$method,
  levels = method_order
)

## ------------------------------------------------------------
## 3. Common journal-style theme
## ------------------------------------------------------------

theme_sim <- theme_classic(base_size = 8) +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    plot.tag = element_text(
      size = 9,
      face = "bold"
    ),
    plot.margin = margin(4, 4, 4, 4)
  )

## ------------------------------------------------------------
## 4. Panel A: Empirical coverage
## ------------------------------------------------------------

pA <- ggplot(
  sim1,
  aes(
    x = method,
    y = coverage,
    fill = method
  )
) +
  geom_boxplot(
    width = 0.62,
    outlier.shape = NA,
    linewidth = 0.35,
    colour = "grey20",
    alpha = 0.80
  )+
  scale_fill_manual(
    values = METHOD_COLS,
    drop = FALSE,
    guide = "none"
  ) +
  geom_hline(
    yintercept = 0.90,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  labs(
    x = NULL,
    y = "Empirical coverage",
    tag = "A"
  ) +
  coord_cartesian(
    ylim = c(0.70, 1.00)
  ) +
  theme_sim +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

## ------------------------------------------------------------
## 5. Panel B: Median lower predictive bound
## ------------------------------------------------------------

pB <- ggplot(
  sim1,
  aes(
    x = method,
    y = med_lpb,
    fill = method
  )
) +
   geom_boxplot(
  width = 0.62,
  outlier.shape = NA,
  linewidth = 0.35,
  colour = "grey20",
  alpha = 0.80
) +
scale_fill_manual(
  values = METHOD_COLS,
  drop = FALSE,
  guide = "none"
) +
  labs(
    x = NULL,
    y = "Median lower predictive bound",
    tag = "B"
  ) +
  theme_sim +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

## ------------------------------------------------------------
## 6. Panel C: Effective sample size
## ------------------------------------------------------------

pC <- ggplot(
  sim1,
  aes(
    x = method,
    y = ess,
    fill = method
  )
) +
  geom_boxplot(
  width = 0.62,
  outlier.shape = NA,
  linewidth = 0.35,
  colour = "grey20",
  alpha = 0.80
) +
scale_fill_manual(
  values = METHOD_COLS,
  drop = FALSE,
  guide = "none"
) +
  labs(
    x = NULL,
    y = "Effective sample size",
    tag = "C"
  ) +
  theme_sim +
  theme(
   axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )

## ------------------------------------------------------------
## 7. Panel D: Maximum raw ICW
## ------------------------------------------------------------

pD <- ggplot(
  sim1,
  aes(
    x = method,
    y = max_raw_wt,
    fill = method
  )
) +
  geom_boxplot(
  width = 0.62,
  outlier.shape = NA,
  linewidth = 0.35,
  colour = "grey20",
  alpha = 0.80
) +
scale_fill_manual(
  values = METHOD_COLS,
  drop = FALSE,
  guide = "none"
) +
  scale_y_log10() +
  labs(
    x = NULL,
    y = "Maximum raw ICW",
    tag = "D"
  ) +
  theme_sim +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )

## ------------------------------------------------------------
## 8. Combine panels
## ------------------------------------------------------------

fig1 <- (pA | pB) /
        (pC | pD)

## Display
fig1

## ------------------------------------------------------------
## 9. Save figure
## ------------------------------------------------------------

dir.create(
  "figures",
  showWarnings = FALSE
)

## PNG for easy viewing
ggsave(
  filename = "figures/Figure1_Simulation1_ablation_corrected.png",
  plot = fig1,
  width = 183,
  height = 138,
  units = "mm",
  dpi = 300
)

## Publication-quality TIFF
ggsave(
  filename = "figures/Figure1_Simulation1_ablation_corrected.tiff",
  plot = fig1,
  width = 183,
  height = 138,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

## Vector PDF for manuscript submission
ggsave(
  filename = "figures/Figure1_Simulation1_ablation_corrected.pdf",
  plot = fig1,
  width = 183,
  height = 138,
  units = "mm"
)


