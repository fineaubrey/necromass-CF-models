# ============================================================
# CHECK UPDATED LHS + GSA RESULTS
# After restricting both MS and GS to empirical P10-P90
# ============================================================

library(tidyverse)
library(here)

options(scipen = 999)

# ------------------------------------------------------------
# Helper: stop if an expected file is missing
# ------------------------------------------------------------

read_checked <- function(path) {
  
  if (!file.exists(path)) {
    stop(
      "File not found:\n",
      path
    )
  }
  
  readr::read_csv(
    path,
    show_col_types = FALSE
  )
}

# ============================================================
# 1. Confirm soil sampling ranges from analytical Dataset 2
# ============================================================

soil <- read_checked(
  here::here(
    "data",
    "derived",
    "global_soils",
    "Dataset2_analytical_sample.csv"
  )
)

soil_ranges <- soil %>%
  summarise(
    
    MS_P10 = quantile(
      MurA,
      0.10,
      na.rm = TRUE
    ),
    
    MS_P90 = quantile(
      MurA,
      0.90,
      na.rm = TRUE
    ),
    
    GS_P10 = quantile(
      GlcN,
      0.10,
      na.rm = TRUE
    ),
    
    GS_P90 = quantile(
      GlcN,
      0.90,
      na.rm = TRUE
    ),
    
    S_min = min(
      SOC,
      na.rm = TRUE
    ),
    
    S_max = max(
      SOC,
      na.rm = TRUE
    )
  )

cat(
  "\n============================================================\n",
  "SOIL SAMPLING DOMAIN\n",
  "============================================================\n"
)

print(soil_ranges)


# ============================================================
# 2. Read LHS summary statistics
# ============================================================

lhs_summary <- read_checked(
  here::here(
    "data",
    "derived",
    "lhs",
    "LHS_summary_statistics.csv"
  )
)

cat(
  "\n============================================================\n",
  "LHS SUMMARY: NF + fnecC MODELS\n",
  "============================================================\n"
)

lhs_summary %>%
  filter(
    Model %in% c(
      "NF",
      "fnecC composite",
      "fnecC trait-resolved"
    )
  ) %>%
  print(
    n = Inf,
    width = Inf
  )


# ============================================================
# 3. Read raw composite and trait-resolved LHS outputs
# ============================================================

f5 <- read_checked(
  here::here(
    "data",
    "derived",
    "lhs",
    "fnecC5_lhs.csv"
  )
)

f10 <- read_checked(
  here::here(
    "data",
    "derived",
    "lhs",
    "fnecC10_lhs.csv"
  )
)

# Verify expected output columns
if (!"fnecC5_pct_uncapped" %in% names(f5)) {
  stop(
    "Column 'fnecC5_pct_uncapped' not found in fnecC5_lhs.csv"
  )
}

if (!"fnecC10_pct_uncapped" %in% names(f10)) {
  stop(
    "Column 'fnecC10_pct_uncapped' not found in fnecC10_lhs.csv"
  )
}


# ============================================================
# 4. Direct fnecC distribution statistics
# ============================================================

fnec_summary <- bind_rows(
  
  f5 %>%
    summarise(
      
      Formulation = "Composite",
      
      N = n(),
      
      Min = min(
        fnecC5_pct_uncapped,
        na.rm = TRUE
      ),
      
      Q1 = quantile(
        fnecC5_pct_uncapped,
        0.25,
        na.rm = TRUE
      ),
      
      Median = median(
        fnecC5_pct_uncapped,
        na.rm = TRUE
      ),
      
      Q3 = quantile(
        fnecC5_pct_uncapped,
        0.75,
        na.rm = TRUE
      ),
      
      Max = max(
        fnecC5_pct_uncapped,
        na.rm = TRUE
      ),
      
      Mean = mean(
        fnecC5_pct_uncapped,
        na.rm = TRUE
      ),
      
      SD = sd(
        fnecC5_pct_uncapped,
        na.rm = TRUE
      )
    ),
  
  
  f10 %>%
    summarise(
      
      Formulation = "Trait-resolved",
      
      N = n(),
      
      Min = min(
        fnecC10_pct_uncapped,
        na.rm = TRUE
      ),
      
      Q1 = quantile(
        fnecC10_pct_uncapped,
        0.25,
        na.rm = TRUE
      ),
      
      Median = median(
        fnecC10_pct_uncapped,
        na.rm = TRUE
      ),
      
      Q3 = quantile(
        fnecC10_pct_uncapped,
        0.75,
        na.rm = TRUE
      ),
      
      Max = max(
        fnecC10_pct_uncapped,
        na.rm = TRUE
      ),
      
      Mean = mean(
        fnecC10_pct_uncapped,
        na.rm = TRUE
      ),
      
      SD = sd(
        fnecC10_pct_uncapped,
        na.rm = TRUE
      )
    )
)

cat(
  "\n============================================================\n",
  "DIRECT fnecC DISTRIBUTION CHECK\n",
  "============================================================\n"
)

print(
  fnec_summary,
  n = Inf,
  width = Inf
)


# ============================================================
# 5. Exceedance rates
# ============================================================

exceedance <- bind_rows(
  
  tibble(
    
    Formulation = "Composite",
    
    `>100%` =
      100 *
      mean(
        f5$fnecC5_pct_uncapped > 100,
        na.rm = TRUE
      ),
    
    `>200%` =
      100 *
      mean(
        f5$fnecC5_pct_uncapped > 200,
        na.rm = TRUE
      ),
    
    `>500%` =
      100 *
      mean(
        f5$fnecC5_pct_uncapped > 500,
        na.rm = TRUE
      )
  ),
  
  
  tibble(
    
    Formulation = "Trait-resolved",
    
    `>100%` =
      100 *
      mean(
        f10$fnecC10_pct_uncapped > 100,
        na.rm = TRUE
      ),
    
    `>200%` =
      100 *
      mean(
        f10$fnecC10_pct_uncapped > 200,
        na.rm = TRUE
      ),
    
    `>500%` =
      100 *
      mean(
        f10$fnecC10_pct_uncapped > 500,
        na.rm = TRUE
      )
  )
)

cat(
  "\n============================================================\n",
  "fnecC EXCEEDANCE RATES\n",
  "============================================================\n"
)

print(
  exceedance,
  n = Inf,
  width = Inf
)


# ============================================================
# 6. Read affected Sobol analyses
# ============================================================

sobol_nf <- read_checked(
  here::here(
    "data",
    "derived",
    "gsa",
    "sobol_NF_saltelli_summary.csv"
  )
)

sobol_f5 <- read_checked(
  here::here(
    "data",
    "derived",
    "gsa",
    "sobol_fnecC5_saltelli_summary.csv"
  )
)

sobol_f10 <- read_checked(
  here::here(
    "data",
    "derived",
    "gsa",
    "sobol_fnecC10_saltelli_summary.csv"
  )
)


# ============================================================
# 7. Rank Sobol results by total-order sensitivity
# ============================================================

sobol_nf_ranked <- sobol_nf %>%
  arrange(
    desc(ST)
  )

sobol_f5_ranked <- sobol_f5 %>%
  arrange(
    desc(ST)
  )

sobol_f10_ranked <- sobol_f10 %>%
  arrange(
    desc(ST)
  )


cat(
  "\n============================================================\n",
  "SOBOL: FUNGAL NECROMASS (NF)\n",
  "============================================================\n"
)

print(
  sobol_nf_ranked,
  n = Inf,
  width = Inf
)


cat(
  "\n============================================================\n",
  "SOBOL: COMPOSITE fnecC\n",
  "============================================================\n"
)

print(
  sobol_f5_ranked,
  n = Inf,
  width = Inf
)


cat(
  "\n============================================================\n",
  "SOBOL: TRAIT-RESOLVED fnecC\n",
  "============================================================\n"
)

print(
  sobol_f10_ranked,
  n = Inf,
  width = Inf
)


# ============================================================
# 8. Compact manuscript-oriented summary
# ============================================================

cat(
  "\n============================================================\n",
  "COMPACT RESULTS FOR MANUSCRIPT CHECK\n",
  "============================================================\n\n"
)

cat("fnecC distributions:\n")

fnec_summary %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 2)
    )
  ) %>%
  print(
    n = Inf,
    width = Inf
  )


cat("\nExceedance rates (% of simulations):\n")

exceedance %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 1)
    )
  ) %>%
  print(
    n = Inf,
    width = Inf
  )


cat("\nNF Sobol total-order indices:\n")

sobol_nf_ranked %>%
  select(
    param,
    S1,
    S1_conf,
    ST,
    ST_conf
  ) %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 3)
    )
  ) %>%
  print(
    n = Inf,
    width = Inf
  )


cat("\nComposite fnecC Sobol total-order indices:\n")

sobol_f5_ranked %>%
  select(
    param,
    S1,
    S1_conf,
    ST,
    ST_conf
  ) %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 3)
    )
  ) %>%
  print(
    n = Inf,
    width = Inf
  )


cat("\nTrait-resolved fnecC Sobol total-order indices:\n")

sobol_f10_ranked %>%
  select(
    param,
    S1,
    S1_conf,
    ST,
    ST_conf
  ) %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 3)
    )
  ) %>%
  print(
    n = Inf,
    width = Inf
  )


cat(
  "\n============================================================\n",
  "CHECK COMPLETE\n",
  "============================================================\n"
)