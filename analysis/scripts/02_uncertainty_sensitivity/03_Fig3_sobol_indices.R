# ============================================================
# Figure 3 — Global sensitivity analysis
#
# Inputs:
#   data/derived/gsa/sobol_CFB_saltelli_summary.csv
#   data/derived/gsa/sobol_NF_saltelli_summary.csv
#   data/derived/gsa/sobol_fnecC5_saltelli_summary.csv
#   data/derived/gsa/sobol_fnecC10_saltelli_summary.csv
#
# Supplemental:
#   data/derived/gsa/sobol_CFF_saltelli_summary.csv
#
# Outputs:
#   analysis/figures/Fig3_panel.png
#   manuscript/figures/Fig3.pdf
#   analysis/figures/FigS_CFF.png
#   manuscript/figures/FigS_CFF.pdf
#   results/tables/TableS3.csv
# ============================================================

options(scipen = 999)

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(here)
})

source(
  here::here(
    "analysis",
    "R",
    "01_setup.R"
  )
)

gsa_dir <- here::here(
  "data",
  "derived",
  "gsa"
)

table_dir <- here::here(
  "results",
  "tables"
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

read_sobol <- function(filename) {

  path <- file.path(
    gsa_dir,
    filename
  )

  if (!file.exists(path)) {
    stop(
      "Missing Sobol result: ",
      path
    )
  }

  readr::read_csv(
    path,
    show_col_types = FALSE
  )
}

prepare_sobol <- function(df) {

  required <- c(
    "param",
    "S1",
    "S1_conf",
    "ST",
    "ST_conf"
  )

  missing <- setdiff(
    required,
    names(df)
  )

  if (length(missing) > 0) {
    stop(
      "Missing columns: ",
      paste(
        missing,
        collapse = ", "
      )
    )
  }

  df %>%
    transmute(
      parameter = as.character(.data$param),
      S1 = as.numeric(.data$S1),
      S1_conf = as.numeric(.data$S1_conf),
      S1_lo95 = .data$S1 - .data$S1_conf,
      S1_hi95 = .data$S1 + .data$S1_conf,
      ST = as.numeric(.data$ST),
      ST_conf = as.numeric(.data$ST_conf),
      ST_lo95 = .data$ST - .data$ST_conf,
      ST_hi95 = .data$ST + .data$ST_conf
    )
}

sobol_cfb <- prepare_sobol(
  read_sobol(
    "sobol_CFB_saltelli_summary.csv"
  )
)

sobol_cff <- prepare_sobol(
  read_sobol(
    "sobol_CFF_saltelli_summary.csv"
  )
)

sobol_nf <- prepare_sobol(
  read_sobol(
    "sobol_NF_saltelli_summary.csv"
  )
)

sobol_f5 <- prepare_sobol(
  read_sobol(
    "sobol_fnecC5_saltelli_summary.csv"
  )
)

sobol_f10 <- prepare_sobol(
  read_sobol(
    "sobol_fnecC10_saltelli_summary.csv"
  )
)

labels_cfb <- c(
  cB = expression(italic(c)[B]),
  MGP = expression(italic(M)[GP]),
  MGN = expression(italic(M)[GN]),
  fGP = expression(italic(f)[GP])
)

labels_cff <- c(
  cF = expression(italic(c)[F]),
  GF = expression(italic(G)[F])
)

labels_nf <- c(
  MS = expression(italic(M)[S]),
  GS = expression(italic(G)[S]),
  cF = expression(italic(c)[F]),
  GF = expression(italic(G)[F]),
  rB = expression(italic(r)[B])
)

labels_f5 <- c(
  MS = expression(italic(M)[S]),
  GS = expression(italic(G)[S]),
  S = expression(italic(S)),
  CFB = expression(italic(CF)[B]),
  CFF = expression(italic(CF)[F])
)

labels_f10 <- c(
  MS = expression(italic(M)[S]),
  GS = expression(italic(G)[S]),
  S = expression(italic(S)),
  cB = expression(italic(c)[B]),
  MGP = expression(italic(M)[GP]),
  MGN = expression(italic(M)[GN]),
  fGP = expression(italic(f)[GP]),
  cF = expression(italic(c)[F]),
  GF = expression(italic(G)[F]),
  rB = expression(italic(r)[B])
)

build_panel <- function(
    df,
    labels,
    title
) {

  df2 <- df %>%
    arrange(
      desc(.data$ST)
    ) %>%
    mutate(
      parameter = factor(
        .data$parameter,
        levels = .data$parameter
      )
    )

  ggplot(
    df2,
    aes(
      x = .data$parameter
    )
  ) +
    geom_col(
      aes(
        y = .data$ST
      ),
      width = 0.62,
      fill = "#61C0B8",
      color = "grey35",
      linewidth = 0.35
    ) +
    geom_errorbar(
      aes(
        ymin = .data$ST_lo95,
        ymax = .data$ST_hi95
      ),
      width = 0.20,
      linewidth = 0.45,
      color = "grey35"
    ) +
    geom_errorbar(
      aes(
        ymin = .data$S1_lo95,
        ymax = .data$S1_hi95
      ),
      width = 0.17,
      linewidth = 0.55,
      color = "black"
    ) +
    geom_point(
      aes(
        y = .data$S1
      ),
      shape = 21,
      size = 1.75,
      stroke = 0.35,
      fill = "black",
      color = "black"
    ) +
    coord_flip() +
    scale_x_discrete(
      labels = labels
    ) +
    scale_y_continuous(
      limits = c(
        0,
        1.05
      ),
      breaks = seq(
        0,
        1,
        by = 0.25
      ),
      expand = expansion(
        mult = c(
          0,
          0.025
        )
      )
    ) +
    labs(
      title = title,
      x = NULL,
      y = "Sobol index"
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      plot.title = element_text(
        size = 10,
        face = "plain",
        hjust = 0,
        margin = margin(
          b = 5
        )
      ),
      axis.text.y = element_text(
        size = 9,
        color = "grey25"
      ),
      axis.text.x = element_text(
        size = 8,
        color = "grey25"
      ),
      axis.title.x = element_text(
        size = 9,
        margin = margin(
          t = 4
        )
      ),
      axis.ticks.y = element_blank(),
      axis.line = element_line(
        linewidth = 0.4,
        color = "grey35"
      )
    )
}

pA <- build_panel(
  sobol_cfb,
  labels_cfb,
  expression(
    "Bacterial conversion factor" ~ (CF[B])
  )
)

pB <- build_panel(
  sobol_nf,
  labels_nf,
  expression(
    "Fungal necromass carbon" ~ (N[F])
  )
)

pC <- build_panel(
  sobol_f5,
  labels_f5,
  expression(
    "Composite" ~ italic(f)[necC]
  )
)

pD <- build_panel(
  sobol_f10,
  labels_f10,
  expression(
    "Trait-resolved" ~ italic(f)[necC]
  )
)

plots <- (
  (pA | pB) /
    (pC | pD)
) +
  plot_layout(
    widths = c(
      1,
      1.15
    ),
    heights = c(
      0.95,
      1.05
    )
  ) +
  plot_annotation(
    tag_levels = "A"
  )

print(
  plots
)

save_figure(
  plot = plots,
  stem = "Fig3",
  width = 9,
  height = 7.5,
  units = "in",
  dpi = 600,
  formats = c(
    "pdf",
    "png"
  )
)

# Supplemental CFF panel
p_cff <- build_panel(
  sobol_cff,
  labels_cff,
  expression(
    "Fungal conversion factor" ~ (CF[F])
  )
)

save_figure(
  plot = p_cff,
  stem = "FigS_CFF",
  width = 4.5,
  height = 3.5,
  units = "in",
  dpi = 600,
  formats = c(
    "pdf",
    "png"
  )
)

TableS3 <- bind_rows(
  sobol_cfb %>%
    mutate(
      model = "CFB"
    ),
  sobol_nf %>%
    mutate(
      model = "NF"
    ),
  sobol_f5 %>%
    mutate(
      model = "fnecC composite"
    ),
  sobol_f10 %>%
    mutate(
      model = "fnecC trait-resolved"
    )
) %>%
  select(
    .data$model,
    .data$parameter,
    .data$S1,
    .data$S1_conf,
    .data$S1_lo95,
    .data$S1_hi95,
    .data$ST,
    .data$ST_conf,
    .data$ST_lo95,
    .data$ST_hi95
  )

readr::write_csv(
  TableS3,
  file.path(
    table_dir,
    "TableS3.csv"
  )
)

message(
  "Saved Figure 3, supplemental CFF figure, and Table S3."
)
