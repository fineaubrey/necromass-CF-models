library(tidyverse)
library(here)

lhs_sum <- readr::read_csv(
  here::here(
    "data",
    "derived",
    "lhs",
    "LHS_summary_statistics.csv"
  ),
  show_col_types = FALSE
)

names(lhs_sum)
print(lhs_sum)

lhs_fnecC_table <- lhs_sum %>%
  filter(
    str_detect(
      tolower(Model),
      "fnecc|composite|trait"
    )
  ) %>%
  select(
    Model,
    N,
    Min,
    Q1,
    Median,
    Q3,
    Max,
    Mean,
    SD
  ) %>%
  mutate(
    across(
      c(Min, Q1, Median, Q3, Max, Mean, SD),
      ~ round(.x, 2)
    )
  )

lhs_fnecC_table
readr::write_csv(
  lhs_fnecC_table,
  here::here(
    "results",
    "tables",
    "LHS_fnecC_summary.csv"
  )
)


library(tidyverse)
library(here)

# ------------------------------------------------------------
# Read regenerated LHS outputs
# ------------------------------------------------------------

lhs_comp <- readr::read_csv(
  here::here(
    "data",
    "derived",
    "lhs",
    "fnecC5_lhs.csv"
  ),
  show_col_types = FALSE
)

lhs_trait <- readr::read_csv(
  here::here(
    "data",
    "derived",
    "lhs",
    "fnecC10_lhs.csv"
  ),
  show_col_types = FALSE
)

# Check names
names(lhs_comp)
names(lhs_trait)
exceedance_table <- tibble(
  Formulation = c(
    "Composite",
    "Trait-resolved"
  ),
  
  `>100%` = c(
    mean(lhs_comp$fnecC5_pct_uncapped > 100) * 100,
    mean(lhs_trait$fnecC10_pct_uncapped > 100) * 100
  ),
  
  `>200%` = c(
    mean(lhs_comp$fnecC5_pct_uncapped > 200) * 100,
    mean(lhs_trait$fnecC10_pct_uncapped > 200) * 100
  ),
  
  `>500%` = c(
    mean(lhs_comp$fnecC5_pct_uncapped > 500) * 100,
    mean(lhs_trait$fnecC10_pct_uncapped > 500) * 100
  )
) %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 1)
    )
  )

exceedance_table

readr::write_csv(
  exceedance_table,
  here::here(
    "results",
    "tables",
    "LHS_fnecC_exceedance_summary.csv"
  )
)