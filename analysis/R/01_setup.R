# -------------------------
# Figure saving helper
# -------------------------

# -------------------------
# Theme (polished, small)
# -------------------------
theme_sbb_small <- function() {
  theme_minimal(base_size = 7, base_family = "sans") %+replace%
    theme(
      plot.title   = element_text(size = 9, face = "plain", hjust = 0),
      axis.title   = element_text(size = 8),
      axis.text    = element_text(size = 7, colour = "black"),
      legend.title = element_text(size = 7),
      legend.text  = element_text(size = 7),
      strip.text   = element_text(size = 7),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.35),
      axis.ticks   = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks.length = unit(1.8, "pt"),
      legend.key.width  = unit(14, "mm"),
      legend.key.height = unit(3, "mm"),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin  = margin(4, 6, 4, 6),
      panel.spacing = unit(6, "pt"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "#EEEEEE")
    )
}

FIG_DIR_ANALYSIS   <- here::here("analysis", "figures")
FIG_DIR_MANUSCRIPT <- here::here("manuscript", "figures")

dir.create(FIG_DIR_ANALYSIS,   recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR_MANUSCRIPT, recursive = TRUE, showWarnings = FALSE)

save_figure <- function(plot,
                        stem,                    # e.g., "Fig2" (no extension)
                        width  = 6.5,
                        height = 4.5,
                        units  = "in",
                        dpi    = 600,
                        formats = c("pdf","png","tiff"),
                        path_analysis   = FIG_DIR_ANALYSIS,
                        path_manuscript = FIG_DIR_MANUSCRIPT,
                        pdf_to_manuscript = TRUE,
                        raster_to_analysis = TRUE) {

  stopifnot(is.character(stem), length(stem) == 1)

  formats <- tolower(formats)
  formats <- intersect(formats, c("pdf","png","tif","tiff","jpeg","jpg"))
  if (length(formats) == 0) stop("No valid formats requested.")

  # normalize tif -> tiff
  formats[formats == "tif"] <- "tiff"

  # ensure dirs
  dir.create(path_analysis,   recursive = TRUE, showWarnings = FALSE)
  dir.create(path_manuscript, recursive = TRUE, showWarnings = FALSE)

  out_files <- character(0)

  for (ext in formats) {
    filename <- paste0(stem, ".", ext)

    # Decide where to write each file
    out_path <- if (ext == "pdf" && isTRUE(pdf_to_manuscript)) {
      file.path(path_manuscript, filename)
    } else if (ext != "pdf" && isTRUE(raster_to_analysis)) {
      file.path(path_analysis, filename)
    } else {
      # fallback: analysis
      file.path(path_analysis, filename)
    }

    if (ext == "pdf") {
      ggplot2::ggsave(
        filename = out_path,
        plot = plot,
        width = width, height = height, units = units,
        device = grDevices::cairo_pdf
      )
    } else if (ext == "png") {
      ggplot2::ggsave(
        filename = out_path,
        plot = plot,
        width = width, height = height, units = units,
        dpi = dpi,
        device = "png"
      )
    } else if (ext == "tiff") {
      ggplot2::ggsave(
        filename = out_path,
        plot = plot,
        width = width, height = height, units = units,
        dpi = dpi,
        device = "tiff",
        compression = "lzw"
      )
    } else if (ext %in% c("jpg","jpeg")) {
      ggplot2::ggsave(
        filename = out_path,
        plot = plot,
        width = width, height = height, units = units,
        dpi = dpi,
        device = "jpeg"
      )
    }

    out_files <- c(out_files, out_path)
  }

  invisible(out_files)
}
