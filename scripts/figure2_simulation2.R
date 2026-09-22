## ============================================================
## Figure 2: Simulation 2 censoring stress
##
## Panels:
## A = Q05 coverage versus censoring
## B = PAC success heatmap
## C = mean effective sample size versus censoring
##
## Uses corrected Simulation 2 raw data.
## ============================================================

library(ggplot2)
library(dplyr)
library(patchwork)
source("R/figure_palette.R")

## ------------------------------------------------------------
## 1. Load corrected Simulation 2 data
## ------------------------------------------------------------


sim2 <- read.csv(
  "results/simulation2_pass3_R100_final/simulation2_raw_corrected.csv"
)

## Basic checks
dim(sim2)
table(sim2$method)

## ------------------------------------------------------------
## 2. Define method order
## ------------------------------------------------------------

method_order2 <- c(
  "H-IPCW",
  "H-AIPCW",
  "Stab-IPCW",
  "Clip-IPCW",
  "Clip-AIPCW"
)

sim2$method <- factor(
  sim2$method,
  levels = method_order2
)

## ------------------------------------------------------------
## 3. Recalculate censoring-stratified summary
## ------------------------------------------------------------

table2_sim2 <- sim2 %>%
  group_by(cens, method) %>%
  summarise(
    mean_cov = mean(
      coverage,
      na.rm = TRUE
    ),

    mcse_cov = sd(
      coverage,
      na.rm = TRUE
    ) / sqrt(
      sum(!is.na(coverage))
    ),

    q05_cov = quantile(
      coverage,
      probs = 0.05,
      na.rm = TRUE
    ),

    q10_cov = quantile(
      coverage,
      probs = 0.10,
      na.rm = TRUE
    ),

    pac_ge_090 = mean(
      coverage >= 0.90,
      na.rm = TRUE
    ),

    pac_ge_088 = mean(
      coverage >= 0.88,
      na.rm = TRUE
    ),

    med_lpb = median(
      med_lpb,
      na.rm = TRUE
    ),

    mean_ess = mean(
      ess,
      na.rm = TRUE
    ),

    mean_maxwt = mean(
      max_raw_wt,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

## ------------------------------------------------------------
## 4. Prepare plotting variables
## ------------------------------------------------------------

table2_sim2 <- table2_sim2 %>%
  mutate(
    method = factor(
      method,
      levels = method_order2
    ),

    cens_f = factor(
      cens,
      levels = c(0.4, 0.6, 0.8),
      labels = c("40%", "60%", "80%")
    )
  )

## ------------------------------------------------------------
## 5. Common theme
## ------------------------------------------------------------

theme_fig2 <- theme_classic(base_size = 8) +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    plot.tag = element_text(
      size = 9,
      face = "bold"
    ),
    plot.margin = margin(4, 4, 4, 4)
  )

## ------------------------------------------------------------
## 6. Panel A: Q05 coverage versus censoring
## ------------------------------------------------------------

p2A <- ggplot(
  table2_sim2,
  aes(
    x = cens * 100,
    y = q05_cov,
    group = method,
    colour = method,
    linetype = method,
    shape = method
  )
) +
  geom_line(
    linewidth = 0.55
  ) +
  geom_point(
    size = 2.2
  ) +
  scale_colour_manual(
  values = METHOD_COLS,
  drop = FALSE
) +
  geom_hline(
    yintercept = 0.90,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_x_continuous(
    breaks = c(40, 60, 80)
  ) +
  coord_cartesian(
    ylim = c(0.70, 0.92)
  ) +
  labs(
    x = "Censoring level (%)",
    y = "Q05 coverage",
    tag = "A"
  ) +
  theme_fig2

## ------------------------------------------------------------
## 7. Panel B: PAC success heatmap
## ------------------------------------------------------------

p2B <- ggplot(
  table2_sim2,
  aes(
    x = cens_f,
    y = method,
    fill = pac_ge_088
  )
) +
  geom_tile(
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        pac_ge_088
      )
    ),
    size = 2.5
  ) +
  scale_fill_gradient(
  low = "#F1F8F6",
  high = "#009E73",
  limits = c(0.65, 1.00),
  guide = "none"
  ) +
  labs(
    x = "Censoring level",
    y = NULL,
    tag = "B"
  ) +
  theme_fig2 +
  theme(
    axis.text.y = element_text(
      size = 7
    )
  )

## ------------------------------------------------------------
## 8. Panel C: ESS versus censoring
##
## Small horizontal offsets are used so methods with nearly
## identical ESS values remain visible rather than plotting
## directly on top of one another.
## ------------------------------------------------------------

p2C_df <- table2_sim2 %>%
  arrange(
    cens,
    method
  )

offset_map <- c(
  "H-IPCW"     = -0.8,
  "H-AIPCW"    = -0.4,
  "Stab-IPCW"  =  0.0,
  "Clip-IPCW"  =  0.4,
  "Clip-AIPCW" =  0.8
)

p2C_df$xpos <-
  p2C_df$cens * 100 +
  offset_map[
    as.character(
      p2C_df$method
    )
  ]

p2C <- ggplot(
  p2C_df,
  aes(
  x = xpos,
  y = mean_ess,
  group = method,
  colour = method,
  linetype = method,
  shape = method
)
) +
  geom_line(
    linewidth = 0.55
  ) +
  geom_point(
    size = 2.2
  ) +
 scale_colour_manual(
  values = METHOD_COLS,
  drop = FALSE
) +
 scale_x_continuous(
    breaks = c(40, 60, 80),
    labels = c(
      "40",
      "60",
      "80"
    )
  ) +
  labs(
    x = "Censoring level (%)",
    y = "Mean effective sample size",
    tag = "C"
  ) +
  theme_fig2

## ------------------------------------------------------------
## 9. Combine panels
## ------------------------------------------------------------

fig2 <- p2A | p2B | p2C

fig2 <- fig2 +
  plot_layout(
    widths = c(
      1.15,
      0.90,
      1.15
    ),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )

## Display
fig2

## ------------------------------------------------------------
## 10. Save figure
## ------------------------------------------------------------

dir.create(
  "figures",
  showWarnings = FALSE
)

## PNG for viewing
ggsave(
  filename = "figures/Figure2_Simulation2_censoring_stress.png",
  plot = fig2,
  width = 183,
  height = 80,
  units = "mm",
  dpi = 300
)

## Vector PDF for manuscript submission
ggsave(
  filename = "figures/Figure2_Simulation2_censoring_stress.pdf",
  plot = fig2,
  width = 183,
  height = 80,
  units = "mm"
)

## Publication-quality TIFF
ggsave(
  filename = "figures/Figure2_Simulation2_censoring_stress.tiff",
  plot = fig2,
  width = 183,
  height = 80,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

## ------------------------------------------------------------
## 11. Save the summary used for Figure 2 and Table 2
## ------------------------------------------------------------

write.csv(
  table2_sim2,
  "results/simulation2_pass3_R100_final/Table2_simulation2_by_censoring_corrected.csv",
  row.names = FALSE
)
