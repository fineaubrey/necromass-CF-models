# ============================================================
# Phase 3 — Paired numerical comparison of empirical outputs
#
# Input:
#   data/derived/global_soils/global_soils_model_estimates_wide.csv
#
# Models:
#   CF2006
#   CF2024
#   Cmass_MurAGlcN
#   Cmass_GlcN
#
# Analysis:
#   Median paired numerical differences with bootstrap 95% CIs
#
# Output:
#   results/tables/global_soils_paired_numerical_differences.csv
#
# Notes:
# - Each soil observation is evaluated with all four methods.
# - Pairwise difference is defined as method1 - method2.
# - Bootstrap CIs use 10,000 resamples of the paired
#   difference vector.
# - Cross-family differences are numerical contrasts, not
#   estimates of measurement bias.
# - No model outputs are capped or filtered here.
# - No null-hypothesis tests are performed.
# ============================================================
options(scipen = 999)

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

infile <- here::here(
  "data",
  "derived",
  "global_soils",
  "global_soils_model_estimates_wide.csv"
)

table_outdir <- here::here(
  "results",
  "tables"
)

dir.create(
  table_outdir,
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

model_cols <- c(
  "CF2006",
  "CF2024",
  "Cmass_MurAGlcN",
  "Cmass_GlcN"
)

missing_cols <- setdiff(
  c("obs_id", model_cols),
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
# Data integrity checks
# ------------------------------------------------------------

if (anyDuplicated(dat$obs_id) > 0) {
  stop(
    "obs_id is not unique. Paired model comparisons require ",
    "one row per soil observation."
  )
}

finite_check <- vapply(
  dat[model_cols],
  function(x) all(is.finite(x)),
  logical(1)
)

if (!all(finite_check)) {
  stop(
    "Non-finite values found in model outputs: ",
    paste(
      names(finite_check)[!finite_check],
      collapse = ", "
    )
  )
}

n_obs <- nrow(dat)
k_approaches <- length(model_cols)

message(
  "Paired soil observations: ",
  n_obs
)

message(
  "Approaches: ",
  paste(
    model_cols,
    collapse = ", "
  )
)

# ============================================================
# Paired numerical differences among model outputs
#
# These contrasts describe numerical differences in outputs
# expressed as percentage points of SOC. Cross-family
# differences are not interpreted as measurement bias because
# CF-based and molecule-level approaches quantify different
# quantities.
#
# Difference orientation:
#
#   method1 - method2
# ============================================================

pair_matrix <- combn(
  model_cols,
  2,
  simplify = FALSE
)

bootstrap_median_difference <- function(
    difference,
    R = 10000,
    seed = 42
) {
  
  difference <- as.numeric(
    difference
  )
  
  difference <- difference[
    is.finite(
      difference
    )
  ]
  
  n <- length(
    difference
  )
  
  if (n < 2) {
    return(
      tibble(
        median_difference = NA_real_,
        median_difference_ci_low = NA_real_,
        median_difference_ci_high = NA_real_
      )
    )
  }
  
  set.seed(
    seed
  )
  
  boot_medians <- replicate(
    R,
    median(
      sample(
        difference,
        size = n,
        replace = TRUE
      )
    )
  )
  
  tibble(
    median_difference = median(
      difference
    ),
    median_difference_ci_low = unname(
      quantile(
        boot_medians,
        0.025,
        names = FALSE
      )
    ),
    median_difference_ci_high = unname(
      quantile(
        boot_medians,
        0.975,
        names = FALSE
      )
    )
  )
}

cf_methods <- c(
  "CF2006",
  "CF2024"
)

mass_methods <- c(
  "Cmass_MurAGlcN",
  "Cmass_GlcN"
)

paired_difference_table <- purrr::map2_dfr(
  pair_matrix,
  seq_along(
    pair_matrix
  ),
  function(
    pair,
    pair_index
  ) {
    
    method1 <- pair[1]
    method2 <- pair[2]
    
    difference <- (
      dat[[method1]] -
        dat[[method2]]
    )
    
    bootstrap_median_difference(
      difference = difference,
      R = 10000,
      seed = 42 + pair_index - 1
    ) %>%
      mutate(
        method1 = method1,
        method2 = method2,
        comparison = paste(
          method1,
          "minus",
          method2
        ),
        N_pairs = length(
          difference
        ),
        comparison_type = case_when(
          method1 %in% cf_methods &
            method2 %in% cf_methods ~
            "within CF family",
          
          method1 %in% mass_methods &
            method2 %in% mass_methods ~
            "within molecule-level family",
          
          TRUE ~
            "cross-family numerical contrast"
        ),
        .before = 1
      )
  }
) %>%
  select(
    comparison,
    comparison_type,
    method1,
    method2,
    N_pairs,
    median_difference,
    median_difference_ci_low,
    median_difference_ci_high
  )

# ============================================================
# Save output
# ============================================================

paired_difference_path <- file.path(
  table_outdir,
  "global_soils_paired_numerical_differences.csv"
)

readr::write_csv(
  paired_difference_table,
  paired_difference_path
)

# ============================================================
# Console output
# ============================================================

message(
  "\n===== Paired numerical differences ====="
)

message(
  "Difference orientation: method1 - method2"
)

print(
  paired_difference_table,
  n = Inf
)

message(
  "\nSaved:\n  ",
  paired_difference_path
)
