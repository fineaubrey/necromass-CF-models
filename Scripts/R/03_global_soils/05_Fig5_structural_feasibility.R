# ============================================================
# Phase 3 — Figure 5: Structural feasibility of CF2006
#
# Figure 5:
# Structural feasibility and emergent behavior of the
# standard CF model (CF2006).
#
# Input:
#   data/derived/global_soils/global_soils_model_estimates_wide.csv
#
# IMPORTANT:
# - Input is the canonical paired analytical sample.
# - MurA and GlcN are already in mg g^-1 soil.
# - NO additional outlier filtering is applied.
# - The canonical CF2006 value generated in Step 02 is retained.
#
# Panel A:
#   Exact full-CF2006 structural relationship:
#
#       N_CF2006 = 45*MS + 9*GS
#
#       fnecC = 100 * N_CF2006 / S
#
#   Therefore:
#
#       N_CF2006 = (fnecC / 100) * S
#
#   The 100% contour is the exact physical boundary for CF2006.
#
# Panel B:
#   Full canonical CF2006 fnecC vs observed soil GlcN,
#   colored by SOC.
#
#   This illustrates that increasing GS does not map uniquely
#   onto increasing fnecC because SOC normalization and MurA
#   also affect the final estimate.
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
    "\nRun 02_Fig4_global_soils_model_comparisons.R first."
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
    "GlcN must be > 0."
  )
}


# ============================================================
# Construct exact CF2006 numerator
#
# CF2006:
#
#   fnecC = 100 * (45*MS + 9*GS) / S
#
# Therefore:
#
#   N_CF2006 = 45*MS + 9*GS
#
# ============================================================

dat <- dat %>%
  mutate(
    CF2006_numerator =
      45 * MurA +
      9 * GlcN,
    
    CF2006_check =
      100 *
      CF2006_numerator /
      SOC,
    
    over100 =
      CF2006 > 100
  )


# ------------------------------------------------------------
# Verify that the formula reproduces canonical CF2006
# ------------------------------------------------------------

max_cf_difference <- max(
  abs(
    dat$CF2006 -
      dat$CF2006_check
  )
)


message(
  "Maximum difference between canonical CF2006 and ",
  "formula-derived CF2006: ",
  signif(
    max_cf_difference,
    6
  )
)


if (max_cf_difference > 1e-6) {
  stop(
    "Canonical CF2006 values do not match ",
    "100*(45*MurA + 9*GlcN)/SOC."
  )
}


# ------------------------------------------------------------
# Diagnostics
# ------------------------------------------------------------

message(
  "Figure 5 observations: ",
  nrow(dat)
)


message(
  "CF2006 observations >100% SOC: ",
  sum(dat$over100),
  " (",
  round(
    100 * mean(dat$over100),
    2
  ),
  "%)"
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
    CF2006_numerator,
    CF2006,
    over100
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
# Panel A — exact full-CF2006 feasibility surface
#
# Full CF2006:
#
#   fnecC =
#       100 * (45*MS + 9*GS) / S
#
# Define:
#
#   N_CF2006 = 45*MS + 9*GS
#
# For a specified fnecC:
#
#   N_CF2006 =
#       (fnecC / 100) * S
#
# Thus the 100% boundary is exactly:
#
#   N_CF2006 = S
#
# ============================================================

levels_f <- c(
  10,
  30,
  50,
  70,
  100
)


soc_grid <- tibble(
  SOC = 10^seq(
    log10(
      min(dat$SOC)
    ),
    log10(
      max(dat$SOC)
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
        
        CF2006_numerator =
          (f / 100) *
          SOC
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


# ------------------------------------------------------------
# Position label on 100% boundary
# ------------------------------------------------------------

lab_100 <- contours_100 %>%
  slice(
    round(
      0.68 *
        n()
    )
  ) %>%
  transmute(
    SOC,
    CF2006_numerator,
    label =
      "italic(f)[necC] == 100*'% limit'"
  )


# ------------------------------------------------------------
# Exact classification check
#
# Every observation above the 100% line should correspond
# exactly to canonical CF2006 > 100%.
# ------------------------------------------------------------

boundary_check <- src_panelA_points %>%
  mutate(
    above_100_line =
      CF2006_numerator > SOC
  )


if (
  !all(
    boundary_check$above_100_line ==
    boundary_check$over100
  )
) {
  stop(
    "Panel A 100% boundary does not exactly reproduce ",
    "the canonical CF2006 >100% classification."
  )
}


message(
  "Panel A exact-boundary check passed: ",
  sum(boundary_check$above_100_line),
  " observations above 100% line."
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
# SOC fill scale for Panel B
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
    name = fill_name,
    
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(
        42,
        "mm"
      ),
      barheight = unit(
        3,
        "mm"
      )
    )
  )
}


# ============================================================
# Panel A
# ============================================================

pA <- ggplot(
  src_panelA_points,
  aes(
    x = SOC,
    y = CF2006_numerator
  )
) +
  
  geom_point(
    shape = 21,
    fill = "grey75",
    color = "black",
    stroke = 0.18,
    alpha = 0.52,
    size = 1.6
  ) +
  
  # 10–70% contours
  geom_line(
    data = contours_other,
    aes(
      x = SOC,
      y = CF2006_numerator,
      linetype = fnecC
    ),
    linewidth = 0.55,
    color = "grey50",
    alpha = 0.9,
    inherit.aes = FALSE
  ) +
  
  # Exact 100% physical boundary
  geom_line(
    data = contours_100,
    aes(
      x = SOC,
      y = CF2006_numerator
    ),
    linewidth = 0.95,
    color = "grey12",
    inherit.aes = FALSE
  ) +
  
  geom_text(
    data = lab_100,
    aes(
      x = SOC,
      y = CF2006_numerator,
      label = label
    ),
    parse = TRUE,
    hjust = 1.02,
    vjust = -0.5,
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
  
  scale_linetype_manual(
    values = c(
      "10" = "dotted",
      "30" = "longdash",
      "50" = "dashed",
      "70" = "dotdash"
    ),
    
    breaks = c(
      "10",
      "30",
      "50",
      "70"
    ),
    
    name = expression(
      "Full CF2006 iso-" *
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
      "Necromass C" ~
        "(mg C " * g^{-1} * " soil)"
    )
  ) +
  
  theme_sbb_small() +
  
  guides(
    linetype = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      keywidth = unit(
        12,
        "mm"
      )
    )
  ) +
  
  theme(
    panel.grid.major =
      element_blank(),
    
    panel.grid.minor =
      element_blank()
  )


# ============================================================
# Panel B
#
# Full CF2006 outcomes vs observed soil GlcN
#
# This panel demonstrates that higher GlcN does not map
# uniquely onto higher fnecC because SOC normalization and
# MurA also influence the final estimate.
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
      paste0(
        scales::comma(x),
        "%"
      )
    }
  ) +
  
  labs(
    y = expression(
      italic(f)[necC] ~
        "(% SOC)"
    ),
    
    x = expression(
      italic(G)[S] ~
        "(mg GlcN " *
        g^{-1} *
        " soil)"
    )
  ) +
  
  theme_sbb_small() +
  
  theme(
    panel.grid.major =
      element_blank(),
    
    panel.grid.minor =
      element_blank()
  )


# ============================================================
# Combine panels
#
# Panel A:
#   exact full-CF2006 feasibility geometry
#
# Panel B:
#   empirical GS–fnecC relationship
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
    legend.box = "vertical",
    legend.box.just = "center",
    
    legend.box.margin = margin(
      t = 6,
      r = 4,
      b = 2,
      l = 4
    ),
    
    legend.margin = margin(
      t = 1,
      r = 2,
      b = 1,
      l = 2
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