# ============================================================
# Shared helpers — microbial trait dataset
#
# Purpose:
#   Define a single canonical import and filtering step for
#   Dataset 1 so all trait analyses use identical observations.
#
# Canonical datasets:
#   - bacterial MurA: N = 444
#   - fungal GlcN:    N = 781, excluding Oomycota
# ============================================================

load_trait_data <- function(
  path = here::here(
    "data",
    "processed",
    "Dataset1.csv"
  ),
  check_n = TRUE
) {

  if (!file.exists(path)) {
    stop(
      "Trait input file not found: ",
      path
    )
  }

  dat <- readr::read_csv(
    path,
    locale = readr::locale(
      encoding = "ISO-8859-1"
    ),
    show_col_types = FALSE
  ) %>%
    dplyr::mutate(
      # Remove a stray encoding artifact present in one
      # Cyanobacteriota label in the source CSV.
      Phylum = stringr::str_remove_all(
        .data$Phylum,
        "\u00C2"
      )
    )

  required_cols <- c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Genus",
    "Species",
    "Gram",
    "MurA",
    "GlcN"
  )

  missing_cols <- setdiff(
    required_cols,
    names(dat)
  )

  if (length(missing_cols) > 0) {
    stop(
      "Dataset1.csv is missing required columns: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }

  bacteria <- dat %>%
    dplyr::filter(
      .data$Kingdom == "Bacteria"
    ) %>%
    dplyr::mutate(
      MurA = as.numeric(
        .data$MurA
      )
    ) %>%
    dplyr::filter(
      is.finite(
        .data$MurA
      )
    )

  fungi <- dat %>%
    dplyr::filter(
      .data$Kingdom == "Fungi",
      .data$Phylum != "Oomycota"
    ) %>%
    dplyr::mutate(
      GlcN = as.numeric(
        .data$GlcN
      )
    ) %>%
    dplyr::filter(
      is.finite(
        .data$GlcN
      )
    )

  if (isTRUE(check_n)) {

    if (nrow(bacteria) != 444) {
      warning(
        "Expected 444 bacterial MurA observations; found ",
        nrow(bacteria),
        "."
      )
    }

    if (nrow(fungi) != 781) {
      warning(
        paste0(
          "Expected 781 fungal GlcN observations after ",
          "excluding Oomycota; found "
        ),
        nrow(fungi),
        "."
      )
    }
  }

  list(
    all = dat,
    bacteria = bacteria,
    fungi = fungi
  )
}
