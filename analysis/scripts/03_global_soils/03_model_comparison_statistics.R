# ============================================================
# Phase 3 — Statistical comparison among empirical calculation
# approaches
#
# Input:
#   data/derived/global_soils/global_soils_model_estimates_wide.csv
#
# Approaches:
#   CF2006
#   CF2024
#   Cmass_MurAGlcN
#   Cmass_GlcN
#
# Analyses:
#   1. Friedman test of calculated outputs expressed as % SOC,
#      treating soil observations as blocks
#   2. Pairwise Wilcoxon signed-rank tests
#   3. Holm correction across all six pairwise comparisons
#   4. Median paired numerical differences with bootstrap 95% CIs
#
# Outputs:
#   results/tables/global_soils_friedman_test.csv
#   results/tables/global_soils_pairwise_comparisons.csv
#
# Notes:
# - All comparisons are paired because each soil observation is
#   evaluated using all four calculation approaches.
# - The CF-based and carbon-mass approaches quantify different
#   quantities. These tests compare their numerical outputs when
#   expressed on the common scale of % SOC; they do not evaluate
#   accuracy or agreement with a common underlying measurand.
# - Pairwise difference is always defined as:
#
#       method1 - method2
#
# - Bootstrap CIs are based on 10,000 resamples of the paired
#   difference vector.
# - No calculated outputs are capped or filtered here.
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
# 1. Friedman test
# ============================================================

# Matrix rows = soil observations; columns = methods.
friedman_matrix <- as.matrix(
  dat[
    model_cols
  ]
)

friedman_fit <- stats::friedman.test(
  friedman_matrix
)

friedman_chisq <- unname(
  friedman_fit$statistic
)

friedman_df <- unname(
  friedman_fit$parameter
)

friedman_p <- friedman_fit$p.value

# Kendall's W for Friedman repeated-measures design:
# W = chi-square / [N * (k - 1)]
kendalls_W <- (
  friedman_chisq /
    (
      n_obs *
        (k_approaches - 1)
    )
)

friedman_table <- tibble(
  N = n_obs,
  Approaches = k_approaches,
  Chi_square = friedman_chisq,
  df = friedman_df,
  p_value = friedman_p,
  Kendalls_W = kendalls_W
)

# ============================================================
# 2. Pairwise Wilcoxon signed-rank tests
# ============================================================

pair_matrix <- combn(
  model_cols,
  2,
  simplify = FALSE
)

pairwise_results <- purrr::map_dfr(
  pair_matrix,
  function(pair) {
    
    method1 <- pair[1]
    method2 <- pair[2]
    
    x <- dat[[method1]]
    y <- dat[[method2]]
    
    diff <- x - y
    
    wt <- suppressWarnings(
      stats::wilcox.test(
        x,
        y,
        paired = TRUE,
        exact = FALSE,
        correct = TRUE,
        alternative = "two.sided"
      )
    )
    
    tibble(
      method1 = method1,
      method2 = method2,
      comparison = paste(
        method1,
        "vs",
        method2
      ),
      N_pairs = length(diff),
      N_nonzero_differences = sum(
        diff != 0
      ),
      Wilcoxon_V = unname(
        wt$statistic
      ),
      p_value_raw = wt$p.value
    )
  }
) %>%
  mutate(
    p_value_Holm = p.adjust(
      p_value_raw,
      method = "holm"
    ),
    significant_0.05 = (
      p_value_Holm < 0.05
    ),
    significant_0.001 = (
      p_value_Holm < 0.001
    )
  )

# ============================================================
# 3. Median paired differences + bootstrap 95% CI
# ============================================================

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
  
  # Use a pair-specific deterministic seed supplied by the caller.
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

difference_results <- purrr::map2_dfr(
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
    
    diff <- (
      dat[[method1]] -
        dat[[method2]]
    )
    
    bootstrap_median_difference(
      difference = diff,
      R = 10000,
      seed = 42 + pair_index - 1
    ) %>%
      mutate(
        method1 = method1,
        method2 = method2,
        comparison = paste(
          method1,
          "vs",
          method2
        ),
        .before = 1
      )
  }
)

# ------------------------------------------------------------
# Combine inferential and paired-difference results
# ------------------------------------------------------------

pairwise_table <- pairwise_results %>%
  left_join(
    difference_results,
    by = c(
      "method1",
      "method2",
      "comparison"
    )
  ) %>%
  select(
    comparison,
    method1,
    method2,
    N_pairs,
    N_nonzero_differences,
    median_difference,
    median_difference_ci_low,
    median_difference_ci_high,
    Wilcoxon_V,
    p_value_raw,
    p_value_Holm,
    significant_0.05,
    significant_0.001
  )

# ============================================================
# Save outputs
# ============================================================

friedman_path <- file.path(
  table_outdir,
  "global_soils_friedman_test.csv"
)

pairwise_path <- file.path(
  table_outdir,
  "global_soils_pairwise_comparisons.csv"
)

readr::write_csv(
  friedman_table,
  friedman_path
)

readr::write_csv(
  pairwise_table,
  pairwise_path
)

# ============================================================
# Console output
# ============================================================

message(
  "\n===== Friedman test ====="
)

print(
  friedman_table
)

message(
  "\n===== Pairwise Wilcoxon signed-rank tests ====="
)

print(
  pairwise_table %>%
    select(
      comparison,
      N_pairs,
      median_difference,
      median_difference_ci_low,
      median_difference_ci_high,
      Wilcoxon_V,
      p_value_raw,
      p_value_Holm,
      significant_0.001
    ),
  n = Inf
)

# ------------------------------------------------------------
# Explicit significance check
# ------------------------------------------------------------

n_sig_005 <- sum(
  pairwise_table$significant_0.05
)

n_sig_001 <- sum(
  pairwise_table$significant_0.001
)

n_comparisons <- nrow(
  pairwise_table
)

message(
  "\nPairwise significance summary:"
)

message(
  "  Holm-adjusted P < 0.05: ",
  n_sig_005,
  " / ",
  n_comparisons
)

message(
  "  Holm-adjusted P < 0.001: ",
  n_sig_001,
  " / ",
  n_comparisons
)

if (
  n_sig_001 ==
  n_comparisons
) {
  message(
    "  All pairwise comparisons are significant at P < 0.001."
  )
} else if (
  n_sig_005 ==
  n_comparisons
) {
  message(
    "  All pairwise comparisons are significant at P < 0.05, ",
    "but not all are significant at P < 0.001."
  )
} else {
  message(
    "  Not all pairwise comparisons are significant at P < 0.05."
  )
}

message(
  "\nSaved:\n  ",
  friedman_path,
  "\n  ",
  pairwise_path
)