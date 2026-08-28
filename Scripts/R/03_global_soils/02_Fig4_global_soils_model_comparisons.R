# ============================================================
# Phase 3 — Empirical model comparison and Figure 4
#
# Input:
#   data/derived/global_soils/Dataset2_analytical_sample.csv
#
# IMPORTANT:
#   MurA and GlcN in this canonical analytical file are already
#   expressed in mg g^-1 soil. DO NOT divide by 1000 again.
#
# Methods:
#   CF2006
#   CF2024
#   Cmass_MurAGlcN
#   Cmass_GlcN
#
# Model choices:
# - CF2006 uses CFB = 45 and CFF = 9; no bacterial GlcN correction.
# - CF2024 reproduces the published Hu et al. (2024) formulation,
#   including the fixed 1.16 bacterial GlcN correction.
# - The 1.16 correction is specific to CF2024 and is NOT used in the
#   five-input composite GSA formulation.
# - Carbon-mass calculations use molecular C fractions:
#       MurA = 0.430
#       GlcN = 0.402
# - Empirical outputs expressed as % SOC are not capped at 100%;
#   exceedances are retained and reported as diagnostics.
#
# Outputs:
#   data/derived/global_soils/global_soils_model_estimates_wide.csv
#   data/derived/global_soils/global_soils_model_estimates_long.csv
#   results/tables/global_soils_model_summary.csv
#   analysis/figures/Fig4_panel.png
#   manuscript/figures/Fig4.pdf
# ============================================================

options(scipen = 999)

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(scales)
  library(ggplot2)
})

# ------------------------------------------------------------
# Theme
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
  "Dataset2_analytical_sample.csv"
)

data_outdir <- here::here(
  "data",
  "derived",
  "global_soils"
)

table_outdir <- here::here(
  "results",
  "tables"
)

fig_dir_analysis <- here::here(
  "analysis",
  "figures"
)

fig_dir_manuscript <- here::here(
  "manuscript",
  "figures"
)

dir.create(
  data_outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  fig_dir_analysis,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  fig_dir_manuscript,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(infile)) {
  stop(
    "Canonical analytical dataset not found: ",
    infile,
    "\nRun 01_prepare_global_soils.R first."
  )
}

# ------------------------------------------------------------
# Read canonical analytical dataset
# ------------------------------------------------------------

dat <- readr::read_csv(
  infile,
  show_col_types = FALSE
)

required_cols <- c(
  "MurA",
  "GlcN",
  "SOC"
)

missing_cols <- setdiff(
  required_cols,
  names(dat)
)

if (length(missing_cols) > 0) {
  stop(
    "Analytical dataset is missing required columns: ",
    paste(
      missing_cols,
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------
# Sanity checks
# ------------------------------------------------------------

if (any(!is.finite(dat$MurA))) {
  stop("Non-finite MurA values found in canonical analytical dataset.")
}

if (any(!is.finite(dat$GlcN))) {
  stop("Non-finite GlcN values found in canonical analytical dataset.")
}

if (any(!is.finite(dat$SOC))) {
  stop("Non-finite SOC values found in canonical analytical dataset.")
}

if (any(dat$SOC <= 0)) {
  stop("SOC <= 0 found in canonical analytical dataset.")
}

# Preserve canonical row identity if available.
if ("Dataset2_row" %in% names(dat)) {
  dat <- dat %>%
    mutate(
      obs_id = Dataset2_row
    )
} else {
  dat <- dat %>%
    mutate(
      obs_id = row_number()
    )
}

message(
  "Canonical soil observations: ",
  nrow(dat)
)

message(
  "MurA range (mg g^-1 soil): ",
  signif(min(dat$MurA), 5),
  " - ",
  signif(max(dat$MurA), 5)
)

message(
  "GlcN range (mg g^-1 soil): ",
  signif(min(dat$GlcN), 5),
  " - ",
  signif(max(dat$GlcN), 5)
)

message(
  "SOC range (mg C g^-1 soil): ",
  signif(min(dat$SOC), 5),
  " - ",
  signif(max(dat$SOC), 5)
)

# ============================================================
# Calculate four empirical model implementations
# ============================================================

dat_wide <- dat %>%
  transmute(
    obs_id,
    MurA,
    GlcN,
    SOC,
    
    # ----------------------------------------------------------
    # Standard conversion-factor model:
    # Appuhn & Joergensen (2006)
    #
    # No bacterial GlcN correction.
    # ----------------------------------------------------------
    
    CF2006 =
      100 *
      (
        45 * MurA +
          9 * GlcN
      ) /
      SOC,
    
    # ----------------------------------------------------------
    # Hu et al. (2024) published alternative formulation
    #
    # This fixed 1.16 correction is retained ONLY because it is
    # part of the published CF2024 implementation being compared
    # empirically. It is not the five-input composite GSA model.
    # ----------------------------------------------------------
    
    CF2024 =
      100 *
      (
        31.3 * MurA +
          10.8 *
          pmax(
            GlcN -
              1.16 * MurA,
            0
          )
      ) /
      SOC,
    
    # ----------------------------------------------------------
    # Direct molecular carbon: MurA + GlcN
    # ----------------------------------------------------------
    
    Cmass_MurAGlcN =
      100 *
      (
        0.430 * MurA +
          0.402 * GlcN
      ) /
      SOC,
    
    # ----------------------------------------------------------
    # Direct molecular carbon: GlcN only
    # ----------------------------------------------------------
    
    Cmass_GlcN =
      100 *
      (
        0.402 * GlcN
      ) /
      SOC
  )

# Check outputs before reshaping.
model_cols <- c(
  "CF2006",
  "CF2024",
  "Cmass_MurAGlcN",
  "Cmass_GlcN"
)

if (
  any(
    !is.finite(
      as.matrix(
        dat_wide[
          model_cols
        ]
      )
    )
  )
) {
  stop(
    "Non-finite empirical model outputs were generated. ",
    "Inspect the canonical analytical dataset before proceeding."
  )
}

# ------------------------------------------------------------
# Long-form data
#
# The shared output column is named generically because the
# CF-based approaches estimate necromass carbon, whereas the
# carbon-mass approaches quantify measured amino-sugar carbon.
# All outputs are expressed as percentages of SOC.
# ------------------------------------------------------------

dat_long <- dat_wide %>%
  pivot_longer(
    cols = all_of(model_cols),
    names_to = "method",
    values_to = "output_pct_SOC"
  ) %>%
  mutate(
    method = factor(
      method,
      levels = model_cols
    )
  )

# ============================================================
# Summary statistics
# ============================================================

sum_F4 <- dat_long %>%
  group_by(method) %>%
  summarise(
    N = n(),
    
    median = median(
      output_pct_SOC
    ),
    
    q25 = quantile(
      output_pct_SOC,
      0.25,
      names = FALSE
    ),
    
    q75 = quantile(
      output_pct_SOC,
      0.75,
      names = FALSE
    ),
    
    p95 = quantile(
      output_pct_SOC,
      0.95,
      names = FALSE
    ),
    
    p99 = quantile(
      output_pct_SOC,
      0.99,
      names = FALSE
    ),
    
    p995 = quantile(
      output_pct_SOC,
      0.995,
      names = FALSE
    ),
    
    over100_n = sum(
      output_pct_SOC > 100
    ),
    
    over100_pct =
      100 *
      mean(
        output_pct_SOC > 100
      ),
    
    maximum = max(
      output_pct_SOC
    ),
    
    .groups = "drop"
  )

# ------------------------------------------------------------
# Bootstrap 95% percentile CI for each method median
# ------------------------------------------------------------

boot_median <- function(
    x,
    indices
) {
  median(
    x[indices]
  )
}

set.seed(42)

median_ci <- dat_long %>%
  group_split(method) %>%
  purrr::map_dfr(
    function(df) {
      
      x <- df$output_pct_SOC
      
      b <- boot::boot(
        data = x,
        statistic = boot_median,
        R = 10000
      )
      
      ci <- boot::boot.ci(
        b,
        type = "perc"
      )
      
      tibble(
        method = unique(df$method),
        
        median_lo =
          ci$percent[4],
        
        median_hi =
          ci$percent[5]
      )
    }
  )

sum_F4 <- sum_F4 %>%
  left_join(
    median_ci,
    by = "method"
  )

# ============================================================
# Save canonical model outputs before plotting
# ============================================================

readr::write_csv(
  dat_wide,
  file.path(
    data_outdir,
    "global_soils_model_estimates_wide.csv"
  )
)

readr::write_csv(
  dat_long,
  file.path(
    data_outdir,
    "global_soils_model_estimates_long.csv"
  )
)

readr::write_csv(
  sum_F4,
  file.path(
    table_outdir,
    "global_soils_model_summary.csv"
  )
)

# ============================================================
# Figure labels
# ============================================================

method_labels_expr <- c(
  expression(CF[2006]),
  expression(CF[2024]),
  expression(plain("C,mass")[plain("MurA + GlcN")]),
  expression(plain("C,mass")[plain("GlcN")])
)

# ============================================================
# Dynamic plotting range
#
# The upper display limit is 110% of the pooled 99.5th
# percentile, with a minimum upper limit of 110% SOC.
# All observations remain in the saved data and numerical
# summary statistics.
# ============================================================

plot_max <- quantile(
  dat_long$output_pct_SOC,
  0.995,
  na.rm = TRUE
)

plot_max <- max(
  plot_max * 1.10,
  110
)

median_label_y <- plot_max * 0.78
ci_label_y <- plot_max * 0.71

# ============================================================
# Plot
# ============================================================

mako_cols <- c(
  "#312142ff",
  "#414388ff",
  "#3482a4ff",
  "#43BBAD"
)

p4 <- ggplot(
  dat_long,
  aes(
    x = method,
    y = output_pct_SOC,
    fill = method
  )
) +
  
  geom_violin(
    width = 0.90,
    alpha = 0.75,
    color = NA,
    trim = FALSE
  ) +
  
  geom_boxplot(
    width = 0.18,
    outlier.shape = 16,
    outlier.size = 0.8,
    fill = "white",
    color = "grey25",
    fatten = 1,
    outlier.alpha = 0.7,
    outlier.colour = "black"
  ) +
  
  geom_hline(
    yintercept = 100,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey40"
  ) +
  
  scale_fill_manual(
    values = mako_cols
  ) +
  
  scale_x_discrete(
    labels = method_labels_expr
  ) +
  
  scale_y_continuous(
    trans = "sqrt",
    
    breaks = c(
      0,
      1,
      10,
      25,
      50,
      100,
      200,
      500,
      1000,
      2000
    ),
    
    limits = c(
      0,
      plot_max
    ),
    
    labels = function(x) {
      paste0(
        scales::comma(x),
        "%"
      )
    },
    
    expand = expansion(
      mult = c(
        0.12,
        0.03
      )
    )
  ) +
  
  theme_sbb_small() +
  
  labs(
    x = NULL,
    y = "Model output (% SOC)"
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  theme(
    legend.position = "none",
    
    panel.grid.major =
      element_line(
        color = "grey96",
        linewidth = 0.25
      ),
    
    panel.grid.minor =
      element_blank(),
    
    panel.border =
      element_blank(),
    
    axis.text.x =
      element_text(
        size = 9,
        colour = "black",
        margin = margin(
          t = 5
        )
      ),
    
    axis.title.y =
      element_text(
        margin = margin(
          r = 6
        )
      ),
    
    plot.margin =
      margin(
        5,
        10,
        10,
        8
      )
  )

# ------------------------------------------------------------
# Percent exceeding the physical SOC boundary
# ------------------------------------------------------------

p4 <- p4 +
  geom_text(
    data = sum_F4,
    
    aes(
      x = method,
      y = 0,
      label =
        paste0(
          round(
            over100_pct,
            1
          ),
          "% > 100%"
        )
    ),
    
    vjust = 1.2,
    size = 2.7,
    fontface = "italic",
    inherit.aes = FALSE
  )

# ------------------------------------------------------------
# Median and bootstrap 95% CI labels
# ------------------------------------------------------------

p4 <- p4 +
  geom_text(
    data = sum_F4,
    
    aes(
      x = method,
      y = median_label_y,
      
      label =
        paste0(
          round(
            median,
            2
          ),
          "%"
        )
    ),
    
    size = 3.0,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  
  geom_text(
    data =
      sum_F4 %>%
      filter(
        !is.na(median_lo)
      ),
    
    aes(
      x = method,
      y = ci_label_y,
      
      label =
        paste0(
          "95% CI: ",
          round(
            median_lo,
            1
          ),
          "\u2013",
          round(
            median_hi,
            1
          ),
          "%"
        )
    ),
    
    size = 2.1,
    inherit.aes = FALSE
  )

# ------------------------------------------------------------
# Display and save
# ------------------------------------------------------------

print(p4)

ggsave(
  filename = here::here(
    "analysis",
    "figures",
    "Fig4_panel.png"
  ),
  plot = p4,
  width = 7,
  height = 5,
  units = "in",
  dpi = 600
)

ggsave(
  filename = here::here(
    "manuscript",
    "figures",
    "Fig4.pdf"
  ),
  plot = p4,
  width = 8,
  height = 4.5,
  units = "in",
  device = grDevices::cairo_pdf
)

# ============================================================
# Print final summary
# ============================================================

message(
  "\n===== Figure 4 / empirical model summary ====="
)

print(
  sum_F4 %>%
    select(
      method,
      N,
      median,
      q25,
      q75,
      median_lo,
      median_hi,
      p95,
      p99,
      p995,
      over100_n,
      over100_pct,
      maximum
    )
)

message(
  "\nSaved canonical empirical outputs to:\n",
  "  ",
  file.path(
    data_outdir,
    "global_soils_model_estimates_wide.csv"
  ),
  "\n  ",
  file.path(
    data_outdir,
    "global_soils_model_estimates_long.csv"
  ),
  "\n  ",
  file.path(
    table_outdir,
    "global_soils_model_summary.csv"
  )
)