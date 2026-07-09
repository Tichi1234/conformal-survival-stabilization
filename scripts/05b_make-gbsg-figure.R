suppressPackageStartupMessages({
  library(ggplot2)
})

TEAL <- "#80B1B3"

g <- read.csv("results/gbsg/gbsg_realdata_raw_R100.csv")

method_order <- c("Naive-Y", "CC", "HT-IPCW", "H-IPCW",
                  "Stab-IPCW", "Clip-IPCW", "H-AIPCW", "Clip-AIPCW")
g$method <- factor(g$method, levels = method_order)

plot_df <- data.frame(
  method = g$method,
  Coverage = as.numeric(g$coverage_ipcw_hajek),
  `Median LPB` = as.numeric(g$median_lpb),
  `Effective sample size` = as.numeric(g$effective_sample_size),
  `log10 maximum weight` = log10(pmax(as.numeric(g$max_raw_weight), 1)),
  check.names = FALSE
)

gbsg_long <- reshape(
  plot_df,
  varying = c("Coverage", "Median LPB", "Effective sample size", "log10 maximum weight"),
  v.names = "value",
  timevar = "metric",
  times = c("Coverage", "Median LPB", "Effective sample size", "log10 maximum weight"),
  direction = "long"
)

gbsg_long <- gbsg_long[is.finite(gbsg_long$value) & !is.na(gbsg_long$method), ]
gbsg_long$metric <- factor(
  gbsg_long$metric,
  levels = c("Coverage", "Median LPB", "Effective sample size", "log10 maximum weight")
)

p <- ggplot(gbsg_long, aes(x = method, y = value)) +
  geom_boxplot(outlier.size = 0.4, fill = TEAL, alpha = 0.85, colour = "grey25") +
  geom_hline(
    data = data.frame(metric = factor("Coverage", levels = levels(gbsg_long$metric)),
                      yintercept = 0.90),
    aes(yintercept = yintercept),
    linetype = "dashed",
    colour = "grey30"
  ) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(
    title = "GBSG real-data application: coverage, informativeness and stability",
    x = "Method",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 12)
  )

ggsave("figures/gbsg_panels.png", p, width = 9, height = 6, dpi = 200)
ggsave("figures/GBSG_final_manuscript.png", p, width = 9, height = 6, dpi = 200)

cat("Wrote:\n")
cat("  figures/gbsg_panels.png\n")
cat("  figures/GBSG_final_manuscript.png\n")
