# ============================================================
# Phase 3 — Agreement analysis: Lin's CCC + Bland–Altman
# CORRECTED CCC EXTRACTION VERSION
#
# Input:
#   data/derived/global_soils/global_soils_model_estimates_wide.csv
#
# Methods:
#   CF2006
#   CF2024
#   Cmass_MurAGlcN
#   Cmass_GlcN
#
# Analyses:
#   1. Lin's concordance correlation coefficient (CCC)
#      for all six unique method pairs
#   2. Bland–Altman analysis for all six unique method pairs
#   3. Bootstrap 95% CIs (10,000 resamples) for CCC, bias,
#      and limits of agreement
#
# Bland–Altman orientation
# ------------------------
# For every comparison:
#
#   difference = method1 - method2
#   pair_mean  = (method1 + method2) / 2
#
# Therefore:
#   positive bias -> method1 tends to give larger estimates
#   negative bias -> method2 tends to give larger estimates
#
# Outputs:
#   results/tables/global_soils_CCC_pairwise.csv
#   results/tables/global_soils_bland_altman_summary.csv
#   data/derived/global_soils/global_soils_bland_altman_long.csv
#   analysis/figures/BlandAltman_grid.png
#   manuscript/figures/BlandAltman_grid.pdf
#
# IMPORTANT
# ---------
# This script does NOT recalculate the four models.
# It uses the canonical estimates created by:
#   02_Fig6_global_soils_model_comparisons.R
# ============================================================

options(scipen = 999)

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(ggplot2)
  library(patchwork)
  library(DescTools)
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
  "global_soils_model_estimates_wide.csv"
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
    "Canonical model-estimate file not found: ",
    infile,
    "\nRun 02_Fig6_global_soils_model_comparisons.R first."
  )
}

# ------------------------------------------------------------
# Read canonical estimates
# ------------------------------------------------------------

wide <- readr::read_csv(
  infile,
  show_col_types = FALSE
)

methods <- c(
  "CF2006",
  "CF2024",
  "Cmass_MurAGlcN",
  "Cmass_GlcN"
)

required_cols <- c(
  "obs_id",
  methods
)

missing_cols <- setdiff(
  required_cols,
  names(wide)
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

if (anyDuplicated(wide$obs_id) > 0) {
  stop(
    "obs_id is not unique. Agreement analyses require one row per soil."
  )
}

finite_check <- vapply(
  wide[methods],
  function(x) all(is.finite(x)),
  logical(1)
)

if (!all(finite_check)) {
  stop(
    "Non-finite model outputs found in: ",
    paste(
      names(finite_check)[!finite_check],
      collapse = ", "
    )
  )
}

message(
  "Paired soil observations: ",
  nrow(wide)
)

# ------------------------------------------------------------
# Human-readable labels
# ------------------------------------------------------------

method_labels <- c(
  CF2006 = "CF2006",
  CF2024 = "CF2024",
  Cmass_MurAGlcN = "Cmass (MurA + GlcN)",
  Cmass_GlcN = "Cmass (GlcN)"
)

# ------------------------------------------------------------
# All six unique pairs
# ------------------------------------------------------------

pair_list <- combn(
  methods,
  2,
  simplify = FALSE
)

# ============================================================
# Lin's concordance correlation coefficient
# ============================================================

# DescTools::CCC() returns rho.c in matrix-like form.
# Extract only the scalar point estimate for each bootstrap replicate.
extract_ccc <- function(x, y) {
  fit <- suppressWarnings(
    DescTools::CCC(x, y)
  )
  
  value <- as.numeric(
    fit$rho.c[1, 1]
  )
  
  if (
    length(value) != 1 ||
    !is.finite(value)
  ) {
    return(NA_real_)
  }
  
  value
}


ccc_pair <- function(
    df,
    method1,
    method2,
    R = 10000,
    seed = 42
) {
  
  x <- as.numeric(df[[method1]])
  y <- as.numeric(df[[method2]])
  
  ok <- is.finite(x) & is.finite(y)
  
  x <- x[ok]
  y <- y[ok]
  
  n <- length(x)
  
  if (n < 10) {
    return(
      tibble(
        method1 = method1,
        method2 = method2,
        N = n,
        CCC = NA_real_,
        CCC_low = NA_real_,
        CCC_high = NA_real_
      )
    )
  }
  
  ccc0 <- extract_ccc(x, y)
  
  set.seed(seed)
  
  boot_vals <- replicate(
    R,
    {
      ii <- sample.int(
        n,
        n,
        replace = TRUE
      )
      
      extract_ccc(
        x[ii],
        y[ii]
      )
    }
  )
  
  boot_vals <- as.numeric(boot_vals)
  boot_vals <- boot_vals[is.finite(boot_vals)]
  
  if (length(boot_vals) >= 100) {
    ci <- quantile(
      boot_vals,
      probs = c(0.025, 0.975),
      names = FALSE,
      na.rm = TRUE
    )
  } else {
    ci <- c(NA_real_, NA_real_)
  }
  
  tibble(
    method1 = method1,
    method2 = method2,
    N = n,
    CCC = ccc0,
    CCC_low = as.numeric(ci[1]),
    CCC_high = as.numeric(ci[2])
  )
}

ccc_tbl <- purrr::map2_dfr(
  pair_list,
  seq_along(pair_list),
  function(
    pair,
    pair_index
  ) {
    
    ccc_pair(
      df = wide,
      method1 = pair[1],
      method2 = pair[2],
      R = 10000,
      seed = 42 + pair_index - 1
    )
  }
) %>%
  mutate(
    comparison = paste(
      recode(
        method1,
        !!!method_labels
      ),
      "vs",
      recode(
        method2,
        !!!method_labels
      )
    ),
    .before = 1
  ) %>%
  arrange(
    desc(CCC)
  )

# ============================================================
# Bland–Altman helpers
# ============================================================

make_ba_data <- function(
    df,
    method1,
    method2
) {
  
  x <- df[[method1]]
  y <- df[[method2]]
  
  ok <- is.finite(x) &
    is.finite(y)
  
  tibble(
    obs_id = df$obs_id[ok],
    method1 = method1,
    method2 = method2,
    comparison = paste(
      recode(
        method1,
        !!!method_labels
      ),
      "vs",
      recode(
        method2,
        !!!method_labels
      )
    ),
    estimate1 = x[ok],
    estimate2 = y[ok],
    pair_mean = (
      x[ok] + y[ok]
    ) / 2,
    difference = (
      x[ok] - y[ok]
    )
  )
}

ba_long <- purrr::map_dfr(
  pair_list,
  function(pair) {
    make_ba_data(
      df = wide,
      method1 = pair[1],
      method2 = pair[2]
    )
  }
)

# ------------------------------------------------------------
# Bland–Altman summary + bootstrap CIs
# ------------------------------------------------------------

ba_summary_pair <- function(
    df_pair,
    R = 10000,
    seed = 42
) {
  
  d <- df_pair$difference
  n <- length(d)
  
  bias <- mean(d)
  sd_diff <- sd(d)
  
  lower_loa <- (
    bias -
      1.96 * sd_diff
  )
  
  upper_loa <- (
    bias +
      1.96 * sd_diff
  )
  
  set.seed(seed)
  
  boot_stats <- replicate(
    R,
    {
      ii <- sample.int(
        n,
        n,
        replace = TRUE
      )
      
      db <- d[ii]
      
      bias_b <- mean(db)
      sd_b <- sd(db)
      
      c(
        bias = bias_b,
        lower_loa = bias_b - 1.96 * sd_b,
        upper_loa = bias_b + 1.96 * sd_b
      )
    }
  )
  
  boot_stats <- t(
    boot_stats
  )
  
  ci_bias <- quantile(
    boot_stats[, "bias"],
    probs = c(
      0.025,
      0.975
    ),
    names = FALSE,
    na.rm = TRUE
  )
  
  ci_lower <- quantile(
    boot_stats[, "lower_loa"],
    probs = c(
      0.025,
      0.975
    ),
    names = FALSE,
    na.rm = TRUE
  )
  
  ci_upper <- quantile(
    boot_stats[, "upper_loa"],
    probs = c(
      0.025,
      0.975
    ),
    names = FALSE,
    na.rm = TRUE
  )
  
  tibble(
    method1 = unique(
      df_pair$method1
    ),
    method2 = unique(
      df_pair$method2
    ),
    comparison = unique(
      df_pair$comparison
    ),
    N = n,
    
    bias = bias,
    bias_ci_low = ci_bias[1],
    bias_ci_high = ci_bias[2],
    
    SD_difference = sd_diff,
    
    lower_LoA = lower_loa,
    lower_LoA_ci_low = ci_lower[1],
    lower_LoA_ci_high = ci_lower[2],
    
    upper_LoA = upper_loa,
    upper_LoA_ci_low = ci_upper[1],
    upper_LoA_ci_high = ci_upper[2],
    
    median_difference = median(
      d
    ),
    
    min_difference = min(
      d
    ),
    
    max_difference = max(
      d
    )
  )
}

ba_summary <- purrr::map2_dfr(
  pair_list,
  seq_along(pair_list),
  function(
    pair,
    pair_index
  ) {
    
    df_pair <- ba_long %>%
      filter(
        method1 == pair[1],
        method2 == pair[2]
      )
    
    ba_summary_pair(
      df_pair = df_pair,
      R = 10000,
      seed = 100 + pair_index - 1
    )
  }
)

# ============================================================
# Save agreement tables
# ============================================================

ccc_path <- file.path(
  table_outdir,
  "global_soils_CCC_pairwise.csv"
)

ba_summary_path <- file.path(
  table_outdir,
  "global_soils_bland_altman_summary.csv"
)

ba_long_path <- file.path(
  data_outdir,
  "global_soils_bland_altman_long.csv"
)

readr::write_csv(
  ccc_tbl,
  ccc_path
)

readr::write_csv(
  ba_summary,
  ba_summary_path
)

readr::write_csv(
  ba_long,
  ba_long_path
)

# ============================================================
# Bland–Altman plots
# ============================================================

plot_bland_altman <- function(
    df_pair,
    summary_row
) {
  
  bias <- summary_row$bias
  lower <- summary_row$lower_LoA
  upper <- summary_row$upper_LoA
  
  ggplot(
    df_pair,
    aes(
      x = pair_mean,
      y = difference
    )
  ) +
    
    geom_point(
      shape = 21,
      size = 1.15,
      alpha = 0.40,
      fill = "grey45",
      color = "grey15",
      stroke = 0.2
    ) +
    
    geom_hline(
      yintercept = 0,
      color = "grey75",
      linewidth = 0.35
    ) +
    
    geom_hline(
      yintercept = bias,
      linewidth = 0.6
    ) +
    
    geom_hline(
      yintercept = upper,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    
    geom_hline(
      yintercept = lower,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    
    labs(
      title = unique(
        df_pair$comparison
      ),
      x = "Mean of paired estimates (% SOC)",
      y = expression(
        Delta * italic(f)[necC] ~ "(% SOC)"
      )
    ) +
    
    theme_sbb_small() +
    
    theme(
      plot.title =
        element_text(
          size = 9,
          hjust = 0
        ),
      
      axis.title =
        element_text(
          size = 8.5
        ),
      
      axis.text =
        element_text(
          size = 8
        ),
      
      panel.grid.minor =
        element_blank(),
      
      plot.margin =
        margin(
          4,
          5,
          4,
          5
        )
    )
}

ba_plots <- purrr::map(
  seq_along(
    pair_list
  ),
  function(i) {
    
    pair <- pair_list[[i]]
    
    df_pair <- ba_long %>%
      filter(
        method1 == pair[1],
        method2 == pair[2]
      )
    
    summary_row <- ba_summary %>%
      filter(
        method1 == pair[1],
        method2 == pair[2]
      )
    
    plot_bland_altman(
      df_pair = df_pair,
      summary_row = summary_row
    )
  }
)

grid_plot <- (
  ba_plots[[1]] |
    ba_plots[[2]] |
    ba_plots[[3]]
) /
  (
    ba_plots[[4]] |
      ba_plots[[5]] |
      ba_plots[[6]]
  ) +
  patchwork::plot_annotation(
    tag_levels = "A"
  )

print(
  grid_plot
)

# ------------------------------------------------------------
# Save Bland–Altman figure
# ------------------------------------------------------------

out_png <- here::here(
  "analysis",
  "figures",
  "BlandAltman_grid.png"
)

out_pdf <- here::here(
  "manuscript",
  "figures",
  "BlandAltman_grid.pdf"
)

ggsave(
  filename = out_png,
  plot = grid_plot,
  width = 10,
  height = 6.5,
  units = "in",
  dpi = 600
)

ggsave(
  filename = out_pdf,
  plot = grid_plot,
  width = 10,
  height = 6.5,
  units = "in",
  device = grDevices::cairo_pdf
)

# ============================================================
# Console output
# ============================================================

message(
  "\n===== Lin's concordance correlation coefficient ====="
)

print(
  ccc_tbl,
  n = Inf
)

message(
  "\n===== Bland–Altman summary ====="
)

print(
  ba_summary %>%
    select(
      comparison,
      N,
      bias,
      bias_ci_low,
      bias_ci_high,
      lower_LoA,
      upper_LoA,
      median_difference
    ),
  n = Inf
)

message(
  "\nSaved:\n  ",
  ccc_path,
  "\n  ",
  ba_summary_path,
  "\n  ",
  ba_long_path,
  "\n  ",
  out_png,
  "\n  ",
  out_pdf
)