# ============================================================
# Phase 3 — Prepare and summarize Dataset 2 for empirical analysis
#
# FINAL VERSION
#
# Source units:
#   MurA = ug MurA g^-1 soil
#   GlcN = ug GlcN g^-1 soil
#   SOC  = mg C g^-1 soil
#
# Canonical derived units:
#   MurA = mg MurA g^-1 soil
#   GlcN = mg GlcN g^-1 soil
#   SOC  = mg C g^-1 soil
#
# Final choices:
# - MurA, GlcN, and SOC must be finite.
# - SOC must be > 0.
# - NO GS >= MS constraint.
# - NO statistical outlier filtering.
# - NO fnecC <= 100% filtering here.
# - NO automatic duplicate removal.
# - Downstream scripts must NOT divide MurA or GlcN by 1000 again.
# ============================================================

library(tidyverse)
library(here)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

infile <- here::here(
  "data",
  "processed",
  "Dataset2.csv"
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

# ------------------------------------------------------------
# Read Dataset 2
# ------------------------------------------------------------

dat <- readr::read_csv(
  infile,
  show_col_types = FALSE
)

required_cols <- c("MurA", "GlcN", "SOC")

missing_cols <- setdiff(
  required_cols,
  names(dat)
)

if (length(missing_cols) > 0) {
  stop(
    "Dataset2.csv is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ------------------------------------------------------------
# Convert units ONCE
#
# Preserve original source values in ug g^-1.
# Replace MurA and GlcN with canonical mg g^-1 values.
# ------------------------------------------------------------

dat <- dat %>%
  mutate(
    Dataset2_row = row_number(),
    
    MurA_ug_g = suppressWarnings(as.numeric(MurA)),
    GlcN_ug_g = suppressWarnings(as.numeric(GlcN)),
    SOC = suppressWarnings(as.numeric(SOC)),
    
    MurA = MurA_ug_g / 1000,
    GlcN = GlcN_ug_g / 1000
  )

cat(
  "\nRaw Dataset 2 rows:",
  nrow(dat),
  "\n"
)

# ------------------------------------------------------------
# Screening audit
# ------------------------------------------------------------

screening_audit <- tibble(
  Metric = c(
    "Raw rows",
    "Non-finite MurA",
    "Non-finite GlcN",
    "Non-finite SOC",
    "SOC <= 0 among finite SOC",
    "MurA < 0 among finite MurA",
    "GlcN < 0 among finite GlcN",
    "Exact duplicate rows in processed input"
  ),
  N = c(
    nrow(dat),
    sum(!is.finite(dat$MurA)),
    sum(!is.finite(dat$GlcN)),
    sum(!is.finite(dat$SOC)),
    sum(is.finite(dat$SOC) & dat$SOC <= 0),
    sum(is.finite(dat$MurA) & dat$MurA < 0),
    sum(is.finite(dat$GlcN) & dat$GlcN < 0),
    sum(
      duplicated(
        dat %>%
          select(
            -Dataset2_row,
            -MurA_ug_g,
            -GlcN_ug_g
          )
      )
    )
  )
)

cat("\n===== Screening audit =====\n")
print(screening_audit, n = Inf)

# ------------------------------------------------------------
# Canonical analytical sample
# ------------------------------------------------------------

dat_analysis <- dat %>%
  filter(
    is.finite(MurA),
    is.finite(GlcN),
    is.finite(SOC),
    SOC > 0
  )

cat(
  "\nAnalytical sample size:",
  nrow(dat_analysis),
  "\n"
)

# ------------------------------------------------------------
# Relationship diagnostics only
#
# NO filtering is based on any of these relationships.
# ------------------------------------------------------------

relationship_diagnostics <- tibble(
  Metric = c(
    "GS < MS",
    "GS >= MS"
  ),
  N = c(
    sum(dat_analysis$GlcN < dat_analysis$MurA),
    sum(dat_analysis$GlcN >= dat_analysis$MurA)
  ),
  Percent = 100 * N / nrow(dat_analysis)
)

cat("\n===== Relationship diagnostics (no filtering) =====\n")
print(relationship_diagnostics, n = Inf)

# ------------------------------------------------------------
# Molecular-carbon feasibility diagnostic
#
# Carbon fractions used in this project:
#   MurA = 0.430
#   GlcN = 0.402
#
# Diagnostic only; does NOT filter rows.
# ------------------------------------------------------------

dat_analysis <- dat_analysis %>%
  mutate(
    AminoSugar_C_mg_g =
      0.430 * MurA +
      0.402 * GlcN,
    
    AminoSugar_C_pct_SOC =
      100 * AminoSugar_C_mg_g / SOC
  )

molecular_C_diagnostics <- tibble(
  Metric = c(
    "Amino-sugar molecular C > SOC",
    "Amino-sugar molecular C <= SOC"
  ),
  N = c(
    sum(dat_analysis$AminoSugar_C_mg_g > dat_analysis$SOC),
    sum(dat_analysis$AminoSugar_C_mg_g <= dat_analysis$SOC)
  ),
  Percent = 100 * N / nrow(dat_analysis)
)

cat("\n===== Molecular-carbon diagnostics (no filtering) =====\n")
print(molecular_C_diagnostics, n = Inf)

# ------------------------------------------------------------
# Negative-value warnings
# ------------------------------------------------------------

if (any(dat_analysis$MurA < 0)) {
  warning(
    "Negative MurA concentrations are present in the analytical sample. ",
    "They have NOT been removed."
  )
}

if (any(dat_analysis$GlcN < 0)) {
  warning(
    "Negative GlcN concentrations are present in the analytical sample. ",
    "They have NOT been removed."
  )
}

# ------------------------------------------------------------
# Descriptive statistics
# ------------------------------------------------------------

summary_one <- function(x, variable, units) {
  tibble(
    Variable = variable,
    Units = units,
    N = sum(is.finite(x)),
    Min = min(x, na.rm = TRUE),
    Q1 = quantile(x, 0.25, na.rm = TRUE, names = FALSE),
    Median = median(x, na.rm = TRUE),
    Q3 = quantile(x, 0.75, na.rm = TRUE, names = FALSE),
    Max = max(x, na.rm = TRUE),
    Mean = mean(x, na.rm = TRUE),
    SD = sd(x, na.rm = TRUE)
  )
}

summary_table <- bind_rows(
  summary_one(
    dat_analysis$MurA,
    "MurA",
    "mg MurA g^-1 soil"
  ),
  summary_one(
    dat_analysis$GlcN,
    "GlcN",
    "mg GlcN g^-1 soil"
  ),
  summary_one(
    dat_analysis$SOC,
    "SOC",
    "mg C g^-1 soil"
  )
)

range_table <- summary_table %>%
  select(
    Variable,
    Units,
    N,
    Minimum = Min,
    Maximum = Max
  )

cat("\n===== Analytical ranges =====\n")
print(range_table, n = Inf)

cat("\n===== Analytical descriptive statistics =====\n")
print(summary_table, n = Inf)

# ------------------------------------------------------------
# Ecosystem and climate counts
# ------------------------------------------------------------

if ("Ecosystem.type" %in% names(dat_analysis)) {
  ecosystem_counts <- dat_analysis %>%
    count(
      Ecosystem.type,
      name = "N",
      sort = TRUE
    ) %>%
    mutate(
      Percent = 100 * N / nrow(dat_analysis)
    )
} else {
  ecosystem_counts <- tibble()
}

if ("Climate" %in% names(dat_analysis)) {
  climate_counts <- dat_analysis %>%
    count(
      Climate,
      name = "N",
      sort = TRUE
    ) %>%
    mutate(
      Percent = 100 * N / nrow(dat_analysis)
    )
} else {
  climate_counts <- tibble()
}

cat("\n===== Ecosystem type =====\n")

if (nrow(ecosystem_counts) > 0) {
  print(ecosystem_counts, n = Inf)
} else {
  cat("Ecosystem.type column not present.\n")
}

cat("\n===== Climate =====\n")

if (nrow(climate_counts) > 0) {
  print(climate_counts, n = Inf)
} else {
  cat("Climate column not present.\n")
}

# ------------------------------------------------------------
# Save canonical outputs
# ------------------------------------------------------------

readr::write_csv(
  dat_analysis,
  file.path(
    data_outdir,
    "Dataset2_analytical_sample.csv"
  )
)

readr::write_csv(
  summary_table,
  file.path(
    table_outdir,
    "Dataset2_analytical_summary.csv"
  )
)

readr::write_csv(
  range_table,
  file.path(
    table_outdir,
    "Dataset2_analytical_ranges.csv"
  )
)

readr::write_csv(
  screening_audit,
  file.path(
    table_outdir,
    "Dataset2_screening_audit.csv"
  )
)

readr::write_csv(
  relationship_diagnostics,
  file.path(
    table_outdir,
    "Dataset2_relationship_diagnostics.csv"
  )
)

readr::write_csv(
  molecular_C_diagnostics,
  file.path(
    table_outdir,
    "Dataset2_molecular_C_diagnostics.csv"
  )
)

if (nrow(ecosystem_counts) > 0) {
  readr::write_csv(
    ecosystem_counts,
    file.path(
      table_outdir,
      "Dataset2_ecosystem_counts.csv"
    )
  )
}

if (nrow(climate_counts) > 0) {
  readr::write_csv(
    climate_counts,
    file.path(
      table_outdir,
      "Dataset2_climate_counts.csv"
    )
  )
}

cat(
  "\nSaved canonical analytical sample to:\n  ",
  file.path(
    data_outdir,
    "Dataset2_analytical_sample.csv"
  ),
  "\n",
  sep = ""
)

cat(
  "\nSaved summary/audit tables to:\n  ",
  table_outdir,
  "\n",
  sep = ""
)