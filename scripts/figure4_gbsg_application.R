
## ============================================================
## Figure 4: GBSG real-data application
## Repeated train/calibration/test splits, R = 100
## ============================================================

library(ggplot2)
library(dplyr)
library(patchwork)

source("R/figure_palette.R")

## ------------------------------------------------------------
## 1. Locate most recent final GBSG R=100 run
## ------------------------------------------------------------

base_dir <- "results/real_data_gbsg_final_R100"

run_dirs <- list.dirs(
  base_dir,
  recursive = FALSE,
  full.names = TRUE
)

if (length(run_dirs) == 0) {
  stop("No GBSG final run directory found.")
}

run_dir <- run_dirs[which.max(file.info(run_dirs)$mtime)]

cat("Using:", run_dir, "\n")

raw_file <- file.path(
  run_dir,
  "gbsg_realdata_raw_R100.csv"
)

if (!file.exists(raw_file)) {
  stop("Cannot find: ", raw_file)
}

raw <- readr::read_csv(
  raw_file,
  show_col_types = FALSE
)

## ------------------------------------------------------------
## 2. Method order
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

raw <- raw %>%
  mutate(
    method = factor(
      method,
      levels = method_order
    )
  )

## ------------------------------------------------------------
## 3. Construct figure summary
##    Coverage and LPB summaries are conditional on genuine
##    feasible calibration solutions.
## ------------------------------------------------------------

fig4_df <- raw %>%
  group_by(method) %>%
  summarise(
    attempted = n(),

    feasible_n =
      sum(feasible_found %in% TRUE, na.rm = TRUE),

    feasible_rate =
      mean(feasible_found %in% TRUE, na.rm = TRUE),

    q05_cov_feasible =
      quantile(
        coverage_ipcw_hajek[
          feasible_found %in% TRUE
        ],
        0.05,
        na.rm = TRUE
      ),

    med_lpb_feasible =
      median(
        median_lpb[
          feasible_found %in% TRUE
        ],
        na.rm = TRUE
      ),

    mean_ess =
      mean(
        effective_sample_size,
        na.rm = TRUE
      ),

    mean_maxwt =
      mean(
        max_raw_weight,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

print(fig4_df)

## ------------------------------------------------------------
## 4. Publication theme
## ------------------------------------------------------------

theme_pub <- theme_classic(base_size = 8) +
  theme(
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8),
    plot.title = element_text(
      size = 9,
      face = "bold",
      hjust = 0
    ),
    axis.text.x = element_text(
      angle = 40,
      hjust = 1
    ),
    plot.margin = margin(4, 5, 4, 4)
  )

## ------------------------------------------------------------
## Panel A: Feasibility
## ------------------------------------------------------------

pA <- ggplot(
  fig4_df,
  aes(
  x = method,
  y = feasible_rate,
  fill = method
)
) +
  geom_col(
    width = 0.72,
    colour = "grey20",
    linewidth = 0.3,
    alpha = 0.85
  ) +
  scale_fill_manual(
    values = METHOD_COLS,
    drop = FALSE,
    guide = "none"
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.35
  ) +
  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, 0.2),
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  labs(
    title = "A",
    x = NULL,
    y = "Feasible calibration (%)"
  ) +
  theme_pub

## ------------------------------------------------------------
## Panel B: Lower-tail coverage among feasible splits
## ------------------------------------------------------------

pB <- ggplot(
  fig4_df,
  aes(
  x = method,
  y = q05_cov_feasible,
  colour = method
)
) +
  geom_point(
    size = 2.8
  ) +
  scale_colour_manual(
    values = METHOD_COLS,
    drop = FALSE,
    guide = "none"
  ) +
  geom_hline(
    yintercept = 0.90,
    linetype = "dashed",
    linewidth = 0.35
  ) +
  scale_y_continuous(
    breaks = seq(0.2, 1.0, 0.1)
  ) +
  coord_cartesian(
    ylim = c(0.20, 1.00)
  ) +
  labs(
    title = "B",
    x = NULL,
    y = expression(Q[0.05]~"held-out coverage")
  ) +
  theme_pub

##-------------------------------------------------------------
## Panel C: Informativeness
## ------------------------------------------------------------

pC <- ggplot(
  fig4_df,
  aes(
  x = method,
  y = med_lpb_feasible,
  fill = method
)
) +
  geom_col(
    width = 0.72,
    colour = "grey20",
    linewidth = 0.3,
    alpha = 0.85
  ) +
  scale_fill_manual(
    values = METHOD_COLS,
    drop = FALSE,
    guide = "none"
  ) +
  labs(
    title = "C",
    x = NULL,
    y = "Median LPB (years)"
  ) +
  theme_pub

## ------------------------------------------------------------
## Panel D: Effective sample size
## ------------------------------------------------------------

pD <- ggplot(
  fig4_df,
  aes(
  x = method,
  y = mean_ess,
  fill = method
)
) +
  geom_col(
    width = 0.72,
    colour = "grey20",
    linewidth = 0.3,
    alpha = 0.85
  ) +
  scale_fill_manual(
    values = METHOD_COLS,
    drop = FALSE,
    guide = "none"
  ) +
  labs(
    title = "D",
    x = NULL,
    y = "Mean effective sample size"
  ) +
  theme_pub

## ------------------------------------------------------------
## 5. Combine
## ------------------------------------------------------------

fig4 <- (pA | pB) /
        (pC | pD)

## ------------------------------------------------------------
## 6. Save
## ------------------------------------------------------------

dir.create(
  "figures",
  showWarnings = FALSE
)

ggsave(
  "figures/Figure4_GBSG_realdata_application.png",
  fig4,
  width = 183,
  height = 145,
  units = "mm",
  dpi = 600
)

ggsave(
  "figures/Figure4_GBSG_realdata_application.pdf",
  fig4,
  width = 183,
  height = 145,
  units = "mm"
)

ggsave(
  "figures/Figure4_GBSG_realdata_application.tiff",
  fig4,
  width = 183,
  height = 145,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

cat("\nFigure 4 saved successfully.\n")
