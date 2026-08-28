# ============================================================
# Figure 2 — Variation in microbial amino-sugar content
#
# Input:
#   data/processed/Dataset1.csv
#
# Outputs:
#   analysis/figures/Fig2.png
#   analysis/figures/Fig2.tiff
#   manuscript/figures/Fig2.pdf
#
# Figure shows:
#   - empirical distributions for major phyla
#   - individual observations
#   - P10-P90 interval
#   - median
#
# The canonical fungal dataset excludes Oomycota.
# ============================================================

library(tidyverse)
library(patchwork)
library(ggridges)
library(here)

options(scipen = 999)

# ------------------------------------------------------------
# Shared setup
# ------------------------------------------------------------

source(
  here::here(
    "analysis",
    "R",
    "01_setup.R"
  )
)

source(
  here::here(
    "analysis",
    "R",
    "02_trait_data.R"
  )
)

message(
  "Project root: ",
  here::here()
)

# ------------------------------------------------------------
# Load canonical datasets
# ------------------------------------------------------------

traits <- load_trait_data()

bac_traits <- traits$bacteria
fung_traits <- traits$fungi

# ------------------------------------------------------------
# Colors / groups
# ------------------------------------------------------------

pal_gram <- c(
  GN = "#43BBAD",
  GP = "#414388FF"
)

gram_map <- c(
  "Pseudomonadota" = "GN",
  "Firmicutes" = "GP",
  "Actinomycetota" = "GP"
)

phyl_A <- c(
  "Pseudomonadota",
  "Firmicutes",
  "Actinomycetota"
)

phyl_B <- c(
  "Ascomycota",
  "Basidiomycota"
)

# ============================================================
# PANEL A — Bacterial MurA
# ============================================================

bac_dat <- bac_traits %>%
  dplyr::filter(
    .data$Phylum %in% phyl_A
  ) %>%
  dplyr::transmute(
    Level = factor(
      .data$Phylum,
      levels = phyl_A
    ),
    MurA = .data$MurA,
    Gram = factor(
      dplyr::recode(
        .data$Phylum,
        !!!gram_map
      ),
      levels = c(
        "GN",
        "GP"
      )
    )
  )

nA <- bac_dat %>%
  dplyr::count(
    .data$Level,
    name = "n"
  )

ylabs_A <- stats::setNames(
  paste0(
    nA$Level,
    " (N=",
    nA$n,
    ")"
  ),
  nA$Level
)

bac_q <- bac_dat %>%
  dplyr::group_by(
    .data$Level
  ) %>%
  dplyr::summarise(
    p10 = stats::quantile(
      .data$MurA,
      0.10,
      na.rm = TRUE
    ),
    p50 = stats::median(
      .data$MurA,
      na.rm = TRUE
    ),
    p90 = stats::quantile(
      .data$MurA,
      0.90,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

pA <- ggplot2::ggplot(
  bac_dat,
  ggplot2::aes(
    x = .data$MurA,
    y = .data$Level,
    fill = .data$Gram
  )
) +
  ggridges::stat_density_ridges(
    jittered_points = TRUE,
    point_shape = "|",
    point_size = 2.5,
    point_alpha = 0.45,
    point_color = "grey30",
    rel_min_height = 0.001,
    scale = 0.95,
    alpha = 0.65
  ) +
  ggplot2::scale_fill_manual(
    values = pal_gram,
    name = " ",
    labels = c(
      "Gram-negative",
      "Gram-positive"
    )
  ) +
  ggplot2::scale_y_discrete(
    labels = ylabs_A,
    expand = ggplot2::expansion(
      mult = c(
        0.05,
        0.04
      )
    )
  ) +
  ggplot2::labs(
    tag = "A",
    y = "Bacterial phylum",
    x = expression(
      "mg MurA " ~
        g^{-1} ~
        " bacterial biomass"
    )
  ) +
  theme_sbb_small() +
  ggplot2::theme(
    legend.position = "top",
    panel.grid.major.x =
      ggplot2::element_blank(),
    panel.grid.minor.x =
      ggplot2::element_blank(),
    axis.title.x =
      ggplot2::element_text(
        margin = ggplot2::margin(
          t = 6
        )
      ),
    axis.title.y =
      ggplot2::element_text(
        margin = ggplot2::margin(
          r = 6
        )
      ),
    panel.border =
      ggplot2::element_blank()
  ) +
  ggplot2::geom_segment(
    data = bac_q,
    ggplot2::aes(
      x = .data$p10,
      xend = .data$p90,
      y = .data$Level,
      yend = .data$Level
    ),
    inherit.aes = FALSE,
    linewidth = 3.5,
    alpha = 0.35,
    color = "grey20"
  ) +
  ggplot2::geom_point(
    data = bac_q,
    ggplot2::aes(
      x = .data$p50,
      y = .data$Level
    ),
    inherit.aes = FALSE,
    shape = 124,
    size = 4,
    stroke = 0.8,
    color = "grey20"
  )

# ============================================================
# PANEL B — Fungal GlcN
# ============================================================

fung_dat <- fung_traits %>%
  dplyr::filter(
    .data$Phylum %in% phyl_B
  ) %>%
  dplyr::transmute(
    Level = factor(
      .data$Phylum,
      levels = phyl_B
    ),
    GlcN = .data$GlcN
  )

nB <- fung_dat %>%
  dplyr::count(
    .data$Level,
    name = "n"
  )

ylabs_B <- stats::setNames(
  paste0(
    nB$Level,
    " (N=",
    nB$n,
    ")"
  ),
  nB$Level
)

fung_q <- fung_dat %>%
  dplyr::group_by(
    .data$Level
  ) %>%
  dplyr::summarise(
    p10 = stats::quantile(
      .data$GlcN,
      0.10,
      na.rm = TRUE
    ),
    p50 = stats::median(
      .data$GlcN,
      na.rm = TRUE
    ),
    p90 = stats::quantile(
      .data$GlcN,
      0.90,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

pB <- ggplot2::ggplot(
  fung_dat,
  ggplot2::aes(
    x = .data$GlcN,
    y = forcats::fct_rev(
      .data$Level
    )
  )
) +
  ggridges::stat_density_ridges(
    jittered_points = TRUE,
    point_shape = "|",
    point_size = 2.5,
    point_alpha = 0.45,
    point_color = "grey30",
    rel_min_height = 0.001,
    scale = 0.95,
    fill = "#3482A4FF",
    alpha = 0.65
  ) +
  ggplot2::scale_y_discrete(
    labels = function(l) {
      ylabs_B[
        as.character(l)
      ]
    },
    expand = ggplot2::expansion(
      mult = c(
        0.05,
        0.04
      )
    )
  ) +
  ggplot2::labs(
    tag = "B",
    y = "Fungal phylum",
    x = expression(
      "mg GlcN " ~
        g^{-1} ~
        " fungal biomass"
    )
  ) +
  theme_sbb_small() +
  ggplot2::theme(
    legend.position = "none",
    panel.grid.major.x =
      ggplot2::element_blank(),
    panel.grid.minor.x =
      ggplot2::element_blank(),
    axis.title.x =
      ggplot2::element_text(
        margin = ggplot2::margin(
          t = 6
        )
      ),
    axis.title.y =
      ggplot2::element_text(
        margin = ggplot2::margin(
          r = 6
        )
      ),
    panel.border =
      ggplot2::element_blank()
  ) +
  ggplot2::geom_segment(
    data = fung_q,
    ggplot2::aes(
      x = .data$p10,
      xend = .data$p90,
      y = .data$Level,
      yend = .data$Level
    ),
    inherit.aes = FALSE,
    linewidth = 3.5,
    alpha = 0.35,
    color = "grey20"
  ) +
  ggplot2::geom_point(
    data = fung_q,
    ggplot2::aes(
      x = .data$p50,
      y = .data$Level
    ),
    inherit.aes = FALSE,
    shape = 124,
    size = 4,
    stroke = 0.8,
    color = "grey20"
  )

# ============================================================
# Combine and save
# ============================================================

p_fig <- pA +
  pB +
  patchwork::plot_layout(
    heights = c(
      3,
      2
    )
  )

print(
  p_fig
)

save_figure(
  plot = p_fig,
  stem = "Fig2",
  width = 8,
  height = 7,
  units = "in",
  dpi = 600,
  formats = c(
    "pdf",
    "png",
    "tiff"
  )
)

message(
  "\nSaved Figure 2."
)
