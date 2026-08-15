# ============================================================
# Phase 3 — Figure 5: Structural feasibility of CF2006
#
# Current manuscript figure:
#   Figure 5. Structural feasibility and emergent behavior
#   of the standard CF model (CF2006).
#
# Input:
#   data/derived/global_soils/global_soils_model_estimates_wide.csv
#
# IMPORTANT:
# - This input is already the canonical paired analytical sample.
# - MurA and GlcN are already in mg g^-1 soil.
# - NO additional outlier filtering is applied here.
# - CF2006 is NOT recalculated; the canonical value from Step 02 is used.
#
# Panel A:
#   Observed GlcN (GS) vs SOC (S), with reference iso-fnecC
#   contours for the dominant fungal term:
#
#       fnecC = 100 * (CFF * GS) / S
#
#   using CFF = 9.
#
#   These contours are a structural reference surface for the
#   fungal scaling term. They are not the exact full CF2006 surface,
#   because the full model also contains 45*MurA.
#
# Panel B:
#   Full canonical CF2006 fnecC vs observed GlcN, colored by SOC.
#
# Outputs:
#   analysis/figures/Fig5.png
#   analysis/figures/Fig5.tiff
#   manuscript/figures/Fig5.pdf
#   data/derived/global_soils/Fig5_panelA_points.csv
#   data/derived/global_soils/Fig5_panelA_contours.csv
#   data/derived/global_soils/Fig5_panelB_source.csv
# ============================================================

options(scipen = 999)

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(grid)
  library(viridis)
})

# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------

source(
  here::here(
    "analysis",
    "R",
    "01_setup.R"
  )
)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

infile <- here::here(
  "data",
  "derived",
  "global_soils",
  "global_soils_model_estimates_wide.csv"
)

data_outdir <- here::here(
  "data",
  "derived",
  "global_soils"
)

dir.create(
  data_outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(infile)) {
  stop(
    "Canonical model-estimate file not found: ",
    infile,
    "\nRun 02_Fig6_global_soils_model_comparisons.R first."
  )
}

# ------------------------------------------------------------
# Read canonical model estimates
# ------------------------------------------------------------

dat <- readr::read_csv(
  infile,
  show_col_types = FALSE
)

required_cols <- c(
  "obs_id",
  "MurA",
  "GlcN",
  "SOC",
  "CF2006"
)

missing_cols <- setdiff(
  required_cols,
  names(dat)
)

if (length(missing_cols) > 0) {
  stop(
    "Input file is missing required columns: ",
    paste(
      missing_cols,
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------
# Integrity checks
# ------------------------------------------------------------

if (anyDuplicated(dat$obs_id) > 0) {
  stop(
    "obs_id is not unique."
  )
}

numeric_check <- dat %>%
  select(
    MurA,
    GlcN,
    SOC,
    CF2006
  )

if (
  any(
    !is.finite(
      as.matrix(
        numeric_check
      )
    )
  )
) {
  stop(
    "Non-finite values found in the canonical Figure 5 input."
  )
}

if (any(dat$SOC <= 0)) {
  stop(
    "SOC must be > 0."
  )
}

if (any(dat$GlcN <= 0)) {
  stop(
    "Panel A uses a log10 y-axis, so GlcN must be > 0."
  )
}

message(
  "Figure 5 observations: ",
  nrow(dat)
)

message(
  "GlcN range (mg g^-1 soil): ",
  signif(
    min(dat$GlcN),
    5
  ),
  " - ",
  signif(
    max(dat$GlcN),
    5
  )
)

message(
  "SOC range (mg C g^-1 soil): ",
  signif(
    min(dat$SOC),
    5
  ),
  " - ",
  signif(
    max(dat$SOC),
    5
  )
)

message(
  "CF2006 range (% SOC): ",
  signif(
    min(dat$CF2006),
    5
  ),
  " - ",
  signif(
    max(dat$CF2006),
    5
  )
)

# ============================================================
# Source data
# ============================================================

src_panelA_points <- dat %>%
  transmute(
    obs_id,
    SOC,
    MurA,
    GS = GlcN,
    CF2006,
    over100 = (
      CF2006 > 100
    )
  )

src_panelB <- dat %>%
  transmute(
    obs_id,
    SOC,
    MurA,
    soilGlcN = GlcN,
    CF2006
  )

# ============================================================
# Panel A: structural reference surface
#
# Dominant fungal term:
#
#   fnecC = 100 * CFF * GS / S
#
# Solve for GS:
#
#   GS = (fnecC / 100) * S / CFF
#
# For CF2006, CFF = 9.
# ============================================================

levels_f <- c(
  10,
  30,
  50,
  70,
  100
)

CFF_std <- 9

soc_grid <- tibble(
  SOC = 10^seq(
    log10(
      min(
        dat$SOC
      )
    ),
    log10(
      max(
        dat$SOC
      )
    ),
    length.out = 300
  )
)

contours <- purrr::map_dfr(
  levels_f,
  function(f) {
    
    soc_grid %>%
      mutate(
        fnecC = f,
        GS = (
          (f / 100) *
            SOC /
            CFF_std
        )
      )
  }
) %>%
  mutate(
    fnecC = factor(
      fnecC,
      levels = levels_f
    )
  )

contours_100 <- contours %>%
  filter(
    fnecC == "100"
  )

contours_other <- contours %>%
  filter(
    fnecC != "100"
  )

lab_100 <- contours_100 %>%
  slice(
    round(
      0.80 *
        n()
    )
  ) %>%
  transmute(
    SOC,
    GS,
    label =
      "italic(f)[necC] == 100*'% (reference limit)'"
  )

# ------------------------------------------------------------
# Save source data
# ------------------------------------------------------------

readr::write_csv(
  src_panelA_points,
  file.path(
    data_outdir,
    "Fig5_panelA_points.csv"
  )
)

readr::write_csv(
  contours,
  file.path(
    data_outdir,
    "Fig5_panelA_contours.csv"
  )
)

readr::write_csv(
  src_panelB,
  file.path(
    data_outdir,
    "Fig5_panelB_source.csv"
  )
)

# ============================================================
# Shared SOC fill scale
# ============================================================

soc_min <- min(
  dat$SOC
)

soc_max <- max(
  dat$SOC
)

soc_break_candidates <- c(
  2,
  5,
  10,
  25,
  50,
  100,
  250,
  500
)

soc_breaks <- soc_break_candidates[
  soc_break_candidates >= soc_min &
    soc_break_candidates <= soc_max
]

fill_name <- expression(
  italic(S) ~
    "(mg C " *
    g^{-1} *
    " soil)"
)

make_soc_fill_scale <- function() {
  scale_fill_viridis_c(
    option = "mako",
    direction = 1,
    begin = 0.2,
    end = 1,
    trans = "log10",
    breaks = soc_breaks,
    limits = c(
      soc_min,
      soc_max
    ),
    oob = scales::squish,
    name = fill_name
  )
}

# ============================================================
# Panel A
# ============================================================

pA <- ggplot(
  src_panelA_points,
  aes(
    x = SOC,
    y = GS
  )
) +
  
  geom_point(
    aes(
      fill = SOC
    ),
    shape = 21,
    color = "black",
    stroke = 0.18,
    alpha = 0.72,
    size = 1.6
  ) +
  
  geom_line(
    data = contours_other,
    aes(
      x = SOC,
      y = GS,
      linetype = fnecC
    ),
    linewidth = 0.55,
    color = "grey50",
    alpha = 0.9,
    inherit.aes = FALSE
  ) +
  
  geom_line(
    data = contours_100,
    aes(
      x = SOC,
      y = GS
    ),
    linewidth = 0.95,
    color = "grey12",
    inherit.aes = FALSE
  ) +
  
  geom_text(
    data = lab_100,
    aes(
      x = SOC,
      y = GS,
      label = label
    ),
    parse = TRUE,
    hjust = 1.02,
    vjust = -0.45,
    size = 2.7,
    color = "grey12",
    inherit.aes = FALSE
  ) +
  
  scale_x_log10(
    labels = scales::label_number(
      accuracy = 1
    )
  ) +
  
  scale_y_log10(
    labels = scales::label_number(
      accuracy = 0.01
    )
  ) +
  
  make_soc_fill_scale() +
  
  scale_linetype_manual(
    values = c(
      "10" = "dotted",
      "30" = "longdash",
      "50" = "dashed",
      "70" = "dotdash"
    ),
    name = expression(
      "Fungal-term iso-" *
        italic(f)[necC] *
        " (%)"
    )
  ) +
  
  labs(
    x = expression(
      italic(S) ~
        "(mg C " *
        g^{-1} *
        " soil)"
    ),
    y = expression(
      italic(G)[S] ~
        "(mg GlcN " *
        g^{-1} *
        " soil)"
    )
  ) +
  
  theme_sbb_small() +
  
  guides(
    fill = guide_colorbar(
      title.position = "top",
      barwidth = unit(
        36,
        "mm"
      ),
      barheight = unit(
        3,
        "mm"
      )
    ),
    
    linetype = guide_legend(
      title.position = "top",
      nrow = 1,
      byrow = TRUE
    )
  ) +
  
  theme(
    panel.grid.major =
      element_blank(),
    panel.grid.minor =
      element_blank()
  )

# ============================================================
# Panel B: full CF2006 outcomes vs GlcN
# ============================================================

x_break_candidates <- c(
  0,
  0.5,
  1,
  2,
  4,
  8,
  16
)

x_breaks <- x_break_candidates[
  x_break_candidates <=
    max(
      src_panelB$soilGlcN
    )
]

pB <- ggplot(
  src_panelB,
  aes(
    x = soilGlcN,
    y = CF2006,
    fill = SOC
  )
) +
  
  geom_point(
    shape = 21,
    color = "black",
    stroke = 0.22,
    size = 1.6,
    alpha = 0.72
  ) +
  
  geom_hline(
    yintercept = 100,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey35"
  ) +
  
  make_soc_fill_scale() +
  
  scale_x_continuous(
    trans = "log1p",
    breaks = x_breaks,
    labels = scales::label_number(
      accuracy = 0.1
    )
  ) +
  
  scale_y_continuous(
    trans = "sqrt",
    breaks = c(
      0,
      25,
      50,
      100,
      200,
      500,
      1000,
      1500
    ),
    labels = function(x) {
      paste0(scales::comma(x), "%")
    }
  ) +
  
  labs(
    y = expression(
      italic(f)[necC] ~ "(% SOC)"
    ),
    x = expression(
      italic(G)[S] ~
        "(mg GlcN " *
        g^{-1} *
        " soil)"
    )
  ) +
  
  theme_sbb_small()
# ============================================================
# Combine panels
# ============================================================

final_fig <- (
  pA /
    pB
) +
  plot_layout(
    guides = "collect",
    heights = c(
      1,
      1
    )
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.margin = margin(
      t = 6
    )
  )

print(
  final_fig
)

# ============================================================
# Save
# ============================================================

save_figure(
  plot = final_fig,
  stem = "Fig5",
  width = 7.2,
  height = 6.4,
  units = "in",
  dpi = 600,
  formats = c(
    "pdf",
    "png",
    "tiff"
  )
)

message(
  "\nSaved Figure 5:\n",
  "  analysis/figures/Fig5.png\n",
  "  analysis/figures/Fig5.tiff\n",
  "  manuscript/figures/Fig5.pdf\n"
)

message(
  "\nSaved Figure 5 source data:\n",
  "  ",
  file.path(
    data_outdir,
    "Fig5_panelA_points.csv"
  ),
  "\n  ",
  file.path(
    data_outdir,
    "Fig5_panelA_contours.csv"
  ),
  "\n  ",
  file.path(
    data_outdir,
    "Fig5_panelB_source.csv"
  )
)