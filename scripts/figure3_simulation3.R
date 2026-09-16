## ============================================================
## Figure 3: Simulation 3 external benchmark
##
## A = Q05 coverage
## B = PAC success with 0.02 tolerance
## C = median lower predictive bound
## ============================================================

library(ggplot2)
library(dplyr)
library(patchwork)

setwd("~/clipipcw_project")

sim3 <- read.csv(
  "results/simulation3_pass3_R100/simulation3_raw.csv"
)

method_order3 <- c(
  "H-IPCW",
  "H-AIPCW",
  "Clip-IPCW",
  "Clip-AIPCW",
  "DR-COSARC"
)

sim3$method <- factor(
  sim3$method,
  levels = method_order3
)

## ------------------------------------------------------------
## Summary over successful evaluations
## ------------------------------------------------------------

fig3_df <- sim3 %>%
  group_by(method) %>%
  summarise(
    attempted = n(),

    successful = sum(
      fit_ok &
      is.finite(coverage)
    ),

    failure_rate = mean(
      !fit_ok
    ),

    q05_cov = quantile(
      coverage,
      0.05,
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

    .groups = "drop"
  )

print(fig3_df)


## ------------------------------------------------------------
## Common theme
## ------------------------------------------------------------

theme_fig3 <- theme_classic(base_size = 8) +
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
## Panel A: Q05 coverage
## ------------------------------------------------------------

p3A <- ggplot(
  fig3_df,
  aes(
    x = method,
    y = q05_cov
  )
) +
  geom_col(
    width = 0.65
  ) +
  geom_hline(
    yintercept = 0.90,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  coord_cartesian(
    ylim = c(0.65, 0.94)
  ) +
  labs(
    x = NULL,
    y = "Q05 coverage",
    tag = "A"
  ) +
  theme_fig3 +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )


## ------------------------------------------------------------
## Panel B: PAC >= 0.88
## ------------------------------------------------------------

p3B <- ggplot(
  fig3_df,
  aes(
    x = method,
    y = pac_ge_088
  )
) +
  geom_col(
    width = 0.65
  ) +
  coord_cartesian(
    ylim = c(0.60, 1.00)
  ) +
  labs(
    x = NULL,
    y = "PAC success (coverage \u2265 0.88)",
    tag = "B"
  ) +
  theme_fig3 +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )


## ------------------------------------------------------------
## Panel C: median LPB
## ------------------------------------------------------------

p3C <- ggplot(
  fig3_df,
  aes(
    x = method,
    y = med_lpb
  )
) +
  geom_col(
    width = 0.65
  ) +
  labs(
    x = NULL,
    y = "Median lower predictive bound",
    tag = "C"
  ) +
  theme_fig3 +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )


## ------------------------------------------------------------
## Combine
## ------------------------------------------------------------

fig3 <- p3A | p3B | p3C

fig3


## ------------------------------------------------------------
## Save
## ------------------------------------------------------------

dir.create(
  "figures",
  showWarnings = FALSE
)

ggsave(
  "figures/Figure3_Simulation3_external_benchmark.png",
  fig3,
  width = 183,
  height = 75,
  units = "mm",
  dpi = 300
)

ggsave(
  "figures/Figure3_Simulation3_external_benchmark.tiff",
  fig3,
  width = 183,
  height = 75,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)
