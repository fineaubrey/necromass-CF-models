# ============================================================
# Trait variability — inferential statistics
#
# Input:
#   data/processed/Dataset1.csv
#
# Outputs:
#   results/tables/trait_gram_test.csv
#   results/tables/trait_kruskal_tests.csv
#   results/tables/trait_dunn_tests.csv
#   results/tables/trait_pairwise_summary.csv
#   results/tables/trait_taxa_included.csv
#   results/tables/trait_gp_kruskal_tests.csv
#   results/tables/trait_gp_dunn_tests.csv
#   results/tables/trait_outlier_flags.csv
#
# Analyses:
#   - Mann-Whitney U test: bacterial MurA, GN vs GP
#   - Kruskal-Wallis tests across taxonomic ranks
#   - Dunn post-hoc tests with Holm adjustment following
#     significant Kruskal-Wallis tests
#   - Rank-biserial correlations for the Mann-Whitney and
#     pairwise Dunn contrasts
#   - Epsilon-squared effect sizes for Kruskal-Wallis tests
#   - GP-only taxonomic comparisons
#   - 3 x median absolute deviation (MAD) outlier flags
#
# Taxonomic minimum sample-size rules follow the manuscript:
#   Phylum  : n > 50
#   Class   : n > 20
#   Order   : n > 20
#   Genus   : n > 10
#   Species : n > 10
#
# All potential MAD outliers are RETAINED in statistical tests.
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
# Paths
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

# ============================================================
# Settings
# ============================================================

alpha <- 0.05

# Bootstrap settings for the primary GN-vs-GP rank-biserial
# correlation. Dunn rank-biserial correlations are reported as
# point estimates to avoid running a separate 10,000-resample
# bootstrap for every taxonomic contrast.
effect_boot_R <- 10000
effect_boot_seed <- 42

rank_thresholds <- c(
  Phylum = 50,
  Class = 20,
  Order = 20,
  Genus = 10,
  Species = 10
)

# ============================================================
# Helper: raw median absolute deviation
#
# The manuscript describes a threshold of three times the
# median absolute deviation. Here MAD is calculated literally as
# median(|x - median(x)|), without the normal-consistency
# scaling factor used by stats::mad() by default.
# ============================================================

flag_mad_outliers <- function(
    data,
    value_col,
    kingdom,
    trait
) {

  x <- data[[value_col]]

  med <- stats::median(
    x,
    na.rm = TRUE
  )

  mad_raw <- stats::median(
    abs(
      x - med
    ),
    na.rm = TRUE
  )

  lower <- med - 3 * mad_raw
  upper <- med + 3 * mad_raw

  data %>%
    dplyr::mutate(
      kingdom = kingdom,
      trait = trait,
      trait_value = .data[[value_col]],
      trait_median = med,
      mad = mad_raw,
      lower_threshold = lower,
      upper_threshold = upper,
      potential_outlier =
        abs(
          .data[[value_col]] - med
        ) > 3 * mad_raw
    ) %>%
    dplyr::filter(
      .data$potential_outlier
    ) %>%
    dplyr::select(
      dplyr::any_of(
        c(
          "Number",
          "AuthorID",
          "Reference",
          "Kingdom",
          "Phylum",
          "Class",
          "Order",
          "Genus",
          "Species",
          "Gram"
        )
      ),
      .data$kingdom,
      .data$trait,
      .data$trait_value,
      .data$trait_median,
      .data$mad,
      .data$lower_threshold,
      .data$upper_threshold,
      .data$potential_outlier
    )
}

# ============================================================
# Helpers: rank-biserial correlation
#
# For two independent groups, the rank-biserial correlation is
# calculated from the Mann-Whitney U statistic:
#
#   r_rb = 2U / (n1 * n2) - 1
#
# Average ranks are used for ties. Positive values indicate that
# values in group 1 tend to exceed values in group 2; negative
# values indicate that group 2 tends to exceed group 1.
# The statistic ranges from -1 to 1.
# ============================================================

rank_biserial_independent <- function(
    x,
    y
) {

  x <- as.numeric(x)
  y <- as.numeric(y)

  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  n1 <- length(x)
  n2 <- length(y)

  if (
    n1 == 0 ||
    n2 == 0
  ) {
    return(NA_real_)
  }

  pooled_ranks <- rank(
    c(
      x,
      y
    ),
    ties.method = "average"
  )

  U1 <- sum(
    pooled_ranks[
      seq_len(n1)
    ]
  ) -
    n1 *
    (
      n1 + 1
    ) /
    2

  r_rb <-
    2 * U1 /
    (
      n1 * n2
    ) -
    1

  # Protect against negligible floating-point excursions beyond
  # the theoretical range.
  max(
    -1,
    min(
      1,
      r_rb
    )
  )
}

bootstrap_rank_biserial <- function(
    x,
    y,
    R = 10000,
    seed = 42
) {

  x <- as.numeric(x)
  y <- as.numeric(y)

  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  n1 <- length(x)
  n2 <- length(y)

  if (
    n1 < 2 ||
    n2 < 2
  ) {
    return(
      tibble::tibble(
        rank_biserial = NA_real_,
        rank_biserial_ci_low = NA_real_,
        rank_biserial_ci_high = NA_real_
      )
    )
  }

  estimate <- rank_biserial_independent(
    x,
    y
  )

  set.seed(seed)

  boot_values <- replicate(
    R,
    rank_biserial_independent(
      sample(
        x,
        size = n1,
        replace = TRUE
      ),
      sample(
        y,
        size = n2,
        replace = TRUE
      )
    )
  )

  ci <- stats::quantile(
    boot_values,
    probs = c(
      0.025,
      0.975
    ),
    names = FALSE,
    na.rm = TRUE
  )

  tibble::tibble(
    rank_biserial = estimate,
    rank_biserial_ci_low = as.numeric(
      ci[1]
    ),
    rank_biserial_ci_high = as.numeric(
      ci[2]
    )
  )
}

# ============================================================
# Helper: Dunn test with Holm correction
#
# Implemented directly to minimize package dependencies.
# Uses pooled ranks, standard tie correction, two-sided
# normal approximation, and Holm-adjusted P values.
# ============================================================

dunn_holm <- function(
    data,
    value_col,
    group_col
) {

  d <- data %>%
    dplyr::select(
      value = dplyr::all_of(
        value_col
      ),
      group = dplyr::all_of(
        group_col
      )
    ) %>%
    dplyr::filter(
      !is.na(.data$value),
      is.finite(.data$value),
      !is.na(.data$group),
      .data$group != ""
    ) %>%
    dplyr::mutate(
      group = as.character(
        .data$group
      ),
      rank = rank(
        .data$value,
        ties.method = "average"
      )
    )

  group_stats <- d %>%
    dplyr::group_by(
      .data$group
    ) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_rank = mean(
        .data$rank
      ),
      .groups = "drop"
    )

  if (nrow(group_stats) < 2) {
    return(
      tibble::tibble()
    )
  }

  n_total <- nrow(d)

  tie_counts <- table(
    d$value
  )

  tie_term <- if (n_total > 1) {
    sum(
      tie_counts^3 -
        tie_counts
    ) /
      (
        12 *
          (n_total - 1)
      )
  } else {
    0
  }

  rank_variance <- (
    n_total *
      (n_total + 1) /
      12
  ) - tie_term

  pairs <- utils::combn(
    group_stats$group,
    2,
    simplify = FALSE
  )

  out <- purrr::map_dfr(
    pairs,
    function(pair) {

      g1 <- group_stats %>%
        dplyr::filter(
          .data$group == pair[1]
        )

      g2 <- group_stats %>%
        dplyr::filter(
          .data$group == pair[2]
        )

      values1 <- d %>%
        dplyr::filter(
          .data$group == pair[1]
        ) %>%
        dplyr::pull(
          .data$value
        )

      values2 <- d %>%
        dplyr::filter(
          .data$group == pair[2]
        ) %>%
        dplyr::pull(
          .data$value
        )

      se <- sqrt(
        rank_variance *
          (
            1 / g1$n +
              1 / g2$n
          )
      )

      z <- (
        g1$mean_rank -
          g2$mean_rank
      ) / se

      p <- 2 * stats::pnorm(
        abs(z),
        lower.tail = FALSE
      )

      r_rb <- rank_biserial_independent(
        values1,
        values2
      )

      tibble::tibble(
        group1 = pair[1],
        group2 = pair[2],
        n1 = g1$n,
        n2 = g2$n,
        mean_rank1 = g1$mean_rank,
        mean_rank2 = g2$mean_rank,
        z = z,
        p = p,
        rank_biserial = r_rb,
        rank_biserial_abs = abs(
          r_rb
        ),
        effect_direction =
          "positive = group1 > group2; negative = group2 > group1"
      )
    }
  )

  out %>%
    dplyr::mutate(
      p_adj = stats::p.adjust(
        .data$p,
        method = "holm"
      ),
      significant = .data$p_adj < alpha
    )
}

# ============================================================
# Helper: taxonomic analysis at one rank
# ============================================================

analyse_rank <- function(
    data,
    value_col,
    rank_col,
    threshold,
    kingdom,
    trait,
    subset_label = "all"
) {

  taxon_counts <- data %>%
    dplyr::filter(
      !is.na(
        .data[[rank_col]]
      ),
      .data[[rank_col]] != "",
      !is.na(
        .data[[value_col]]
      ),
      is.finite(
        .data[[value_col]]
      )
    ) %>%
    dplyr::count(
      taxon = .data[[rank_col]],
      name = "n"
    ) %>%
    dplyr::mutate(
      kingdom = kingdom,
      trait = trait,
      subset = subset_label,
      rank = rank_col,
      threshold_rule = paste0(
        "n > ",
        threshold
      ),
      included = .data$n > threshold
    ) %>%
    dplyr::select(
      .data$kingdom,
      .data$trait,
      .data$subset,
      .data$rank,
      .data$taxon,
      .data$n,
      .data$threshold_rule,
      .data$included
    )

  eligible_taxa <- taxon_counts %>%
    dplyr::filter(
      .data$included
    ) %>%
    dplyr::pull(
      .data$taxon
    )

  d <- data %>%
    dplyr::filter(
      .data[[rank_col]] %in%
        eligible_taxa,
      !is.na(
        .data[[value_col]]
      ),
      is.finite(
        .data[[value_col]]
      )
    )

  if (
    length(
      unique(
        d[[rank_col]]
      )
    ) < 2
  ) {

    return(
      list(
        counts = taxon_counts,
        kruskal = tibble::tibble(),
        dunn = tibble::tibble()
      )
    )
  }

  formula_obj <- stats::reformulate(
    rank_col,
    response = value_col
  )

  kw <- stats::kruskal.test(
    formula_obj,
    data = d
  )

  k <- dplyr::n_distinct(
    d[[rank_col]]
  )

  n_total <- nrow(d)

  h <- unname(
    kw$statistic
  )

  epsilon_squared <- if (
    n_total > k
  ) {
    max(
      0,
      (
        h - k + 1
      ) /
        (
          n_total - k
        )
    )
  } else {
    NA_real_
  }

  kw_tbl <- tibble::tibble(
    kingdom = kingdom,
    trait = trait,
    subset = subset_label,
    rank = rank_col,
    threshold_rule = paste0(
      "n > ",
      threshold
    ),
    n_taxa = k,
    n_observations = n_total,
    statistic_H = h,
    df = unname(
      kw$parameter
    ),
    p = kw$p.value,
    epsilon_squared = epsilon_squared,
    significant = kw$p.value < alpha
  )

  # Dunn post-hoc comparisons are performed only when the
  # omnibus Kruskal-Wallis test is significant.
  if (kw$p.value < alpha) {

    dunn_tbl <- dunn_holm(
      d,
      value_col = value_col,
      group_col = rank_col
    ) %>%
      dplyr::mutate(
        kingdom = kingdom,
        trait = trait,
        subset = subset_label,
        rank = rank_col,
        threshold_rule = paste0(
          "n > ",
          threshold
        ),
        .before = 1
      )

  } else {

    dunn_tbl <- tibble::tibble()

  }

  list(
    counts = taxon_counts,
    kruskal = kw_tbl,
    dunn = dunn_tbl
  )
}

# ============================================================
# Potential outliers
# ============================================================

trait_outlier_flags <- dplyr::bind_rows(

  flag_mad_outliers(
    bac_traits,
    value_col = "MurA",
    kingdom = "Bacteria",
    trait = "MurA"
  ),

  flag_mad_outliers(
    fung_traits,
    value_col = "GlcN",
    kingdom = "Fungi",
    trait = "GlcN"
  )

)

# ============================================================
# Gram-status comparison
# ============================================================

bac_gram <- bac_traits %>%
  dplyr::filter(
    .data$Gram %in% c(
      "GN",
      "GP"
    )
  ) %>%
  dplyr::mutate(
    Gram = factor(
      .data$Gram,
      levels = c(
        "GN",
        "GP"
      )
    )
  )

gram_test <- stats::wilcox.test(
  MurA ~ Gram,
  data = bac_gram,
  exact = FALSE
)

# The effect is oriented GN relative to GP to match the order in
# the comparison label and the Mann-Whitney U statistic.
gram_effect <- bootstrap_rank_biserial(
  x = bac_gram %>%
    dplyr::filter(
      .data$Gram == "GN"
    ) %>%
    dplyr::pull(
      .data$MurA
    ),
  y = bac_gram %>%
    dplyr::filter(
      .data$Gram == "GP"
    ) %>%
    dplyr::pull(
      .data$MurA
    ),
  R = effect_boot_R,
  seed = effect_boot_seed
)

gram_summary <- bac_gram %>%
  dplyr::group_by(
    .data$Gram
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    median = stats::median(
      .data$MurA
    ),
    mean = mean(
      .data$MurA
    ),
    .groups = "drop"
  )

trait_gram_test <- tibble::tibble(
  comparison = "GN vs GP",
  trait = "MurA",
  test = "Mann-Whitney U",
  n_GN = gram_summary$n[
    gram_summary$Gram == "GN"
  ],
  n_GP = gram_summary$n[
    gram_summary$Gram == "GP"
  ],
  median_GN = gram_summary$median[
    gram_summary$Gram == "GN"
  ],
  median_GP = gram_summary$median[
    gram_summary$Gram == "GP"
  ],
  statistic_U = unname(
    gram_test$statistic
  ),
  p = gram_test$p.value,
  rank_biserial =
    gram_effect$rank_biserial,
  rank_biserial_ci_low =
    gram_effect$rank_biserial_ci_low,
  rank_biserial_ci_high =
    gram_effect$rank_biserial_ci_high,
  rank_biserial_abs = abs(
    gram_effect$rank_biserial
  ),
  effect_direction =
    "positive = GN > GP; negative = GP > GN",
  significant = gram_test$p.value < alpha
)

# ============================================================
# Taxonomic-rank analyses
# ============================================================

analyse_dataset <- function(
    data,
    value_col,
    kingdom,
    trait,
    subset_label = "all",
    ranks = names(
      rank_thresholds
    )
) {

  purrr::map(
    ranks,
    function(rank_col) {
      analyse_rank(
        data = data,
        value_col = value_col,
        rank_col = rank_col,
        threshold =
          rank_thresholds[
            rank_col
          ],
        kingdom = kingdom,
        trait = trait,
        subset_label = subset_label
      )
    }
  ) %>%
    stats::setNames(
      ranks
    )
}

bac_results <- analyse_dataset(
  data = bac_traits,
  value_col = "MurA",
  kingdom = "Bacteria",
  trait = "MurA"
)

fung_results <- analyse_dataset(
  data = fung_traits,
  value_col = "GlcN",
  kingdom = "Fungi",
  trait = "GlcN"
)

trait_taxa_included <- dplyr::bind_rows(
  purrr::map_dfr(
    bac_results,
    "counts"
  ),
  purrr::map_dfr(
    fung_results,
    "counts"
  )
)

trait_kruskal_tests <- dplyr::bind_rows(
  purrr::map_dfr(
    bac_results,
    "kruskal"
  ),
  purrr::map_dfr(
    fung_results,
    "kruskal"
  )
)

trait_dunn_tests <- dplyr::bind_rows(
  purrr::map_dfr(
    bac_results,
    "dunn"
  ),
  purrr::map_dfr(
    fung_results,
    "dunn"
  )
)

trait_pairwise_summary <- trait_dunn_tests %>%
  dplyr::group_by(
    .data$kingdom,
    .data$trait,
    .data$subset,
    .data$rank,
    .data$threshold_rule
  ) %>%
  dplyr::summarise(
    n_pairwise = dplyr::n(),
    n_significant = sum(
      .data$significant,
      na.rm = TRUE
    ),
    percent_significant =
      100 *
      .data$n_significant /
      .data$n_pairwise,
    median_abs_rank_biserial = stats::median(
      abs(
        .data$rank_biserial
      ),
      na.rm = TRUE
    ),
    median_abs_rank_biserial_significant = if (
      any(
        .data$significant,
        na.rm = TRUE
      )
    ) {
      stats::median(
        abs(
          .data$rank_biserial[
            .data$significant
          ]
        ),
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    maximum_abs_rank_biserial = max(
      abs(
        .data$rank_biserial
      ),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# ============================================================
# GP-only taxonomic analyses
#
# These provide a direct check of taxonomic variation within
# Gram-positive bacteria, including the phylum-level
# Actinomycetota vs Firmicutes and class-level
# Actinomycetia vs Bacilli comparisons discussed in the text.
# ============================================================

gp_traits <- bac_traits %>%
  dplyr::filter(
    .data$Gram == "GP"
  )

gp_results <- analyse_dataset(
  data = gp_traits,
  value_col = "MurA",
  kingdom = "Bacteria",
  trait = "MurA",
  subset_label = "GP only"
)

trait_gp_kruskal_tests <- purrr::map_dfr(
  gp_results,
  "kruskal"
)

trait_gp_dunn_tests <- purrr::map_dfr(
  gp_results,
  "dunn"
)

# ============================================================
# Save outputs
# ============================================================

readr::write_csv(
  trait_gram_test,
  file.path(
    out_dir,
    "trait_gram_test.csv"
  )
)

readr::write_csv(
  trait_kruskal_tests,
  file.path(
    out_dir,
    "trait_kruskal_tests.csv"
  )
)

readr::write_csv(
  trait_dunn_tests,
  file.path(
    out_dir,
    "trait_dunn_tests.csv"
  )
)

readr::write_csv(
  trait_pairwise_summary,
  file.path(
    out_dir,
    "trait_pairwise_summary.csv"
  )
)

readr::write_csv(
  trait_taxa_included,
  file.path(
    out_dir,
    "trait_taxa_included.csv"
  )
)

readr::write_csv(
  trait_gp_kruskal_tests,
  file.path(
    out_dir,
    "trait_gp_kruskal_tests.csv"
  )
)

readr::write_csv(
  trait_gp_dunn_tests,
  file.path(
    out_dir,
    "trait_gp_dunn_tests.csv"
  )
)

readr::write_csv(
  trait_outlier_flags,
  file.path(
    out_dir,
    "trait_outlier_flags.csv"
  )
)

# ============================================================
# Console report
# ============================================================

message(
  "\nGram-status test:"
)

print(
  trait_gram_test
)

message(
  "\nKruskal-Wallis tests:"
)

print(
  trait_kruskal_tests
)

message(
  "\nPairwise Dunn summary:"
)

print(
  trait_pairwise_summary
)

message(
  "\nGP-only Kruskal-Wallis tests:"
)

print(
  trait_gp_kruskal_tests
)

message(
  "\nPotential MAD outliers flagged (retained): ",
  nrow(
    trait_outlier_flags
  )
)

message(
  "\nSaved inferential trait outputs to results/tables/."
)