# ============================================================
# FIGURE 5
# Structural feasibility and emergent behavior of CF2006
# ============================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(viridis)
library(here)

# -------------------------------------------------------------------
# Load empirical global soil dataset
# -------------------------------------------------------------------

# EDIT THIS PATH if your processed Dataset 2 has a different filename
soil <- read.csv(
  here::here("data", "processed", "Dataset2.csv")
)

# Expected columns:
# MS = soil MurA concentration (mg g-1 soil)
# GS = soil GlcN concentration (mg g-1 soil)
# S  = soil organic C concentration (mg C g-1 soil)

# -------------------------------------------------------------------
# Calculate CF2006 necromass fraction
# Eq. 4:
#
# fnecC = 100 * (45*MS + 9*GS) / S
# -------------------------------------------------------------------

soil_fig5 <- soil %>%
  filter(
    is.finite(MS),
    is.finite(GS),
    is.finite(S),
    MS >= 0,
    GS > 0,
    S > 0
  ) %>%
  mutate(
    fnecC_CF2006 = 100 * (45 * MS + 9 * GS) / S,
    infeasible   = fnecC_CF2006 > 100
  )

# Quick diagnostic
soil_fig5 %>%
  summarise(
    n = n(),
    median_fnecC = median(fnecC_CF2006, na.rm = TRUE),
    pct_gt_100 = mean(infeasible, na.rm = TRUE) * 100,
    max_fnecC = max(fnecC_CF2006, na.rm = TRUE)
  )

# ============================================================
# PANEL A
# Analytical iso-fnecC contours in GS-S space
# ============================================================

# For the fungal CF term:
#
# fnecC = 100 * CFF * GS / S
#
# with CFF = 9 for CF2006.
#
# Rearranging:
#
# S = (100 * CFF * GS) / fnecC
#
# These lines therefore show the geometric scaling imposed by
# the fungal CF component of CF2006.
# -------------------------------------------------------------------

contour_levels <- c(10, 30, 50, 70, 100)

# Range based on observed data
GS_seq <- 10^seq(
  log10(min(soil_fig5$GS, na.rm = TRUE)),
  log10(max(soil_fig5$GS, na.rm = TRUE)),
  length.out = 500
)

iso_lines <- crossing(
  GS = GS_seq,
  fnecC = contour_levels
) %>%
  mutate(
    S = (100 * 9 * GS) / fnecC,
    contour = factor(
      fnecC,
      levels = contour_levels,
      labels = paste0(contour_levels, "%")
    )
  )

# Label positions toward upper portion of each contour
iso_labels <- iso_lines %>%
  group_by(contour) %>%
  slice_min(
    abs(
      log10(GS) -
        quantile(log10(GS), 0.78)
    ),
    n = 1
  ) %>%
  ungroup()


p5a <- ggplot() +
  
  # Observed soils
  geom_point(
    data = soil_fig5,
    aes(
      x = GS,
      y = S,
      color = S
    ),
    size = 1.5,
    alpha = 0.55
  ) +
  
  # 10-70% contours
  geom_line(
    data = iso_lines %>% filter(fnecC < 100),
    aes(
      x = GS,
      y = S,
      group = contour
    ),
    linewidth = 0.55,
    linetype = "dashed",
    color = "grey35"
  ) +
  
  # 100% physical boundary
  geom_line(
    data = iso_lines %>% filter(fnecC == 100),
    aes(
      x = GS,
      y = S
    ),
    linewidth = 1.0,
    color = "black"
  ) +
  
  # Contour labels
  geom_text(
    data = iso_labels,
    aes(
      x = GS,
      y = S,
      label = contour
    ),
    size = 3,
    hjust = -0.15,
    vjust = -0.2,
    color = "black"
  ) +
  
  scale_x_log10() +
  scale_y_log10() +
  
  scale_color_viridis_c(
    option = "mako",
    trans = "log10",
    name = expression(
      italic(S) ~ "(mg C " * g^{-1} * " soil)"
    )
  ) +
  
  labs(
    x = expression(
      italic(G)[S] ~ "(mg GlcN " * g^{-1} * " soil)"
    ),
    y = expression(
      italic(S) ~ "(mg C " * g^{-1} * " soil)"
    )
  ) +
  
  theme_classic() +
  theme_sbb_small() +
  
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )


# ============================================================
# PANEL B
# Observed CF2006 fnecC versus soil GlcN
# ============================================================

p5b <- ggplot(
  soil_fig5,
  aes(
    x = GS,
    y = fnecC_CF2006,
    color = S
  )
) +
  
  geom_point(
    size = 1.6,
    alpha = 0.6
  ) +
  
  # Physical upper bound
  geom_hline(
    yintercept = 100,
    linewidth = 0.7,
    linetype = "dashed",
    color = "black"
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = 100,
    label = "100% SOC",
    hjust = 1.05,
    vjust = -0.4,
    size = 3
  ) +
  
  scale_x_log10() +
  
  scale_color_viridis_c(
    option = "mako",
    trans = "log10",
    name = expression(
      italic(S) ~ "(mg C " * g^{-1} * " soil)"
    )
  ) +
  
  labs(
    x = expression(
      italic(G)[S] ~ "(mg GlcN " * g^{-1} * " soil)"
    ),
    y = expression(
      italic(f)[necC] ~ "(% of SOC)"
    )
  ) +
  
  theme_classic() +
  theme_sbb_small() +
  
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )


# ============================================================
# Combine panels
# ============================================================

fig5 <- p5a / p5b +
  plot_layout(
    heights = c(1, 1),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    legend.position = "bottom"
  )

fig5


# ============================================================
# Save
# ============================================================

dir.create(
  here::here("analysis", "figures"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  here::here("manuscript", "figures"),
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  filename = here::here(
    "analysis",
    "figures",
    "Fig5_panel.png"
  ),
  plot = fig5,
  width = 5,
  height = 8,
  dpi = 600
)

ggsave(
  filename = here::here(
    "manuscript",
    "figures",
    "Fig5.pdf"
  ),
  plot = fig5,
  width = 5,
  height = 8,
  device = grDevices::cairo_pdf
)

message(
  "✓ Saved Figure 5: analysis/figures/Fig5_panel.png & manuscript/figures/Fig5.pdf"
)