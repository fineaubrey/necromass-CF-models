# ============================================================
# Trait variability — descriptive summaries
#
# Input:
#   data/processed/Dataset1.csv
#
# Outputs:
#   results/tables/trait_overall_summary.csv
#   results/tables/trait_parameter_ranges.csv
#   results/tables/trait_group_summary.csv
#
# Purpose:
#   1. Load the canonical bacterial and fungal trait datasets.
#   2. Calculate descriptive statistics and bootstrap 95% CIs.
#   3. Save P10-P90 ranges used in downstream uncertainty and
#      global sensitivity analyses.
#
# Fungal analyses exclude Oomycota, yielding N = 781.
# ============================================================

library(tidyverse)
library(here)

options(scipen = 999)

# ------------------------------------------------------------
# Shared trait-data definition
# ------------------------------------------------------------

source(
  here::here(
    "scripts",
    "R",
    "01_traits",
    "02_trait_data.R"
  )
)

# ------------------------------------------------------------
# Reproducibility
# ------------------------------------------------------------

SEED <- 42
BOOT_R <- 10000

# ------------------------------------------------------------
# Output directory
# ------------------------------------------------------------

out_dir <- here::here(
  "results",
  "tables"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Load canonical datasets
# ------------------------------------------------------------

traits <- load_trait_data()

bac_traits <- traits$bacteria
fung_traits <- traits$fungi

message(
  "Bacterial MurA N = ",
  nrow(bac_traits)
)

message(
  "Fungal GlcN N = ",
  nrow(fung_traits),
  " (Oomycota excluded)"
)

# ============================================================
# Helper functions
# ============================================================

bootstrap_ci <- function(
  x,
  stat_fun,
  R = BOOT_R,
  seed = SEED
) {

  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0) {
    return(
      c(
        lower = NA_real_,
        upper = NA_real_
      )
    )
  }

  set.seed(seed)

  boot_stats <- replicate(
    R,
    stat_fun(
      sample(
        x,
        size = length(x),
        replace = TRUE
      )
    )
  )

  stats::quantile(
    boot_stats,
    probs = c(
      0.025,
      0.975
    ),
    na.rm = TRUE,
    names = FALSE
  ) %>%
    stats::setNames(
      c(
        "lower",
        "upper"
      )
    )
}

summarize_trait <- function(
  x,
  kingdom,
  trait,
  units = "mg g^-1 biomass"
) {

  x <- x[
    is.finite(x)
  ]

  mean_ci <- bootstrap_ci(
    x,
    mean
  )

  median_ci <- bootstrap_ci(
    x,
    median
  )

  tibble::tibble(
    kingdom = kingdom,
    trait = trait,
    units = units,
    n = length(x),
    min = min(x),
    p10 = stats::quantile(
      x,
      0.10,
      names = FALSE
    ),
    mean = mean(x),
    mean_ci_lower = unname(
      mean_ci["lower"]
    ),
    mean_ci_upper = unname(
      mean_ci["upper"]
    ),
    median = stats::median(x),
    median_ci_lower = unname(
      median_ci["lower"]
    ),
    median_ci_upper = unname(
      median_ci["upper"]
    ),
    p90 = stats::quantile(
      x,
      0.90,
      names = FALSE
    ),
    max = max(x),
    sd = stats::sd(x),
    q25 = stats::quantile(
      x,
      0.25,
      names = FALSE
    ),
    q75 = stats::quantile(
      x,
      0.75,
      names = FALSE
    ),
    iqr = stats::IQR(x)
  )
}

# ============================================================
# Overall summaries
# ============================================================

trait_overall_summary <- dplyr::bind_rows(

  summarize_trait(
    bac_traits$MurA,
    kingdom = "Bacteria",
    trait = "MurA"
  ),

  summarize_trait(
    fung_traits$GlcN,
    kingdom = "Fungi",
    trait = "GlcN"
  )

)

# ============================================================
# P10-P90 parameter ranges used downstream
# ============================================================

bacterial_ranges <- bac_traits %>%
  dplyr::filter(
    .data$Gram %in% c(
      "GN",
      "GP"
    )
  ) %>%
  dplyr::group_by(
    .data$Gram
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    p10 = stats::quantile(
      .data$MurA,
      0.10,
      na.rm = TRUE
    ),
    median = stats::median(
      .data$MurA,
      na.rm = TRUE
    ),
    p90 = stats::quantile(
      .data$MurA,
      0.90,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    parameter = dplyr::recode(
      .data$Gram,
      GN = "MGN",
      GP = "MGP"
    ),
    description = dplyr::recode(
      .data$Gram,
      GN = "Gram-negative bacterial MurA content",
      GP = "Gram-positive bacterial MurA content"
    ),
    group = dplyr::recode(
      .data$Gram,
      GN = "Gram-negative bacteria",
      GP = "Gram-positive bacteria"
    ),
    trait = "MurA"
  ) %>%
  dplyr::select(
    .data$parameter,
    .data$description,
    .data$group,
    .data$trait,
    .data$n,
    .data$p10,
    .data$median,
    .data$p90
  )

fungal_range <- fung_traits %>%
  dplyr::summarise(
    parameter = "GF",
    description = "Fungal GlcN content",
    group = "Fungi",
    trait = "GlcN",
    n = dplyr::n(),
    p10 = stats::quantile(
      .data$GlcN,
      0.10,
      na.rm = TRUE
    ),
    median = stats::median(
      .data$GlcN,
      na.rm = TRUE
    ),
    p90 = stats::quantile(
      .data$GlcN,
      0.90,
      na.rm = TRUE
    )
  )

trait_parameter_ranges <- dplyr::bind_rows(
  bacterial_ranges,
  fungal_range
)

# ============================================================
# Descriptive group summaries
# ============================================================

bacterial_group_summary <- bac_traits %>%
  dplyr::filter(
    .data$Gram %in% c(
      "GN",
      "GP"
    )
  ) %>%
  dplyr::group_by(
    .data$Gram
  ) %>%
  dplyr::summarise(
    kingdom = "Bacteria",
    group_type = "Gram",
    group = dplyr::first(
      .data$Gram
    ),
    trait = "MurA",
    n = dplyr::n(),
    mean = mean(
      .data$MurA,
      na.rm = TRUE
    ),
    median = stats::median(
      .data$MurA,
      na.rm = TRUE
    ),
    p10 = stats::quantile(
      .data$MurA,
      0.10,
      na.rm = TRUE
    ),
    p90 = stats::quantile(
      .data$MurA,
      0.90,
      na.rm = TRUE
    ),
    min = min(
      .data$MurA,
      na.rm = TRUE
    ),
    max = max(
      .data$MurA,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

fungal_group_summary <- fung_traits %>%
  dplyr::group_by(
    .data$Phylum
  ) %>%
  dplyr::summarise(
    kingdom = "Fungi",
    group_type = "Phylum",
    group = dplyr::first(
      .data$Phylum
    ),
    trait = "GlcN",
    n = dplyr::n(),
    mean = mean(
      .data$GlcN,
      na.rm = TRUE
    ),
    median = stats::median(
      .data$GlcN,
      na.rm = TRUE
    ),
    p10 = stats::quantile(
      .data$GlcN,
      0.10,
      na.rm = TRUE
    ),
    p90 = stats::quantile(
      .data$GlcN,
      0.90,
      na.rm = TRUE
    ),
    min = min(
      .data$GlcN,
      na.rm = TRUE
    ),
    max = max(
      .data$GlcN,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

trait_group_summary <- dplyr::bind_rows(
  bacterial_group_summary,
  fungal_group_summary
)

# ============================================================
# Save outputs
# ============================================================

readr::write_csv(
  trait_overall_summary,
  file.path(
    out_dir,
    "trait_overall_summary.csv"
  )
)

readr::write_csv(
  trait_parameter_ranges,
  file.path(
    out_dir,
    "trait_parameter_ranges.csv"
  )
)

readr::write_csv(
  trait_group_summary,
  file.path(
    out_dir,
    "trait_group_summary.csv"
  )
)

message(
  "\nSaved descriptive trait outputs to results/tables/."
)
