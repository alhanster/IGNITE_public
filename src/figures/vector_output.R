# Shared PNG and PDF export for the stage-2 renderers, sourced alongside palette.R.
#
# Writes a 300 dpi ragg PNG for review and a vector PDF for typesetting through one entry point.
# Paths are passed in full rather than as a stem, since supplementary output goes to
# final_plots/supplementary/, a location a stem cannot express.

.repo_rel <- function(path) sub(".*/(final_plots/)", "\\1", path)

# The committed figure_data/ tables were produced on Apple Silicon macOS, the only platform where they
# reproduce exactly. Checks that tie a table to that run fail there and warn elsewhere.
# See REPRODUCIBILITY.md.
on_reference_platform <- function() {
  si <- Sys.info()
  identical(si[["sysname"]], "Darwin") && identical(si[["machine"]], "arm64")
}
check_or_warn <- function(msg) {
  if (on_reference_platform()) stop(msg, call. = FALSE)
  message("WARNING (not the reference platform): ", msg)
}

.quartz_pdf <- function(filename, width, height, ...) {
  grDevices::quartz(type = "pdf", file = filename, width = width, height = height, ...)
}

save_figure <- function(plot, out_png, out_pdf, width_cm, height_cm) {
  dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "cairo"
  ggplot2::ggsave(out_png, plot, width = width_cm, height = height_cm, units = "cm",
                  dpi = 300, bg = "white", device = dev)
  cat("wrote", .repo_rel(out_png), "\n")

  # PDFs need macOS's quartz device; everywhere else the PNG is the whole deliverable.
  # IGNITE_SKIP_PDF=1 forces this path on macOS so it can be exercised there.
  if (!isTRUE(capabilities("aqua")) || nzchar(Sys.getenv("IGNITE_SKIP_PDF"))) {
    # A stale PDF from an earlier render must not ship beside a fresh PNG.
    if (file.exists(out_pdf)) unlink(out_pdf)
    message("skipped ", .repo_rel(out_pdf), ": the PDF deliverables need macOS's ",
            "quartz(type = \"pdf\") device",
            if (nzchar(Sys.getenv("IGNITE_SKIP_PDF"))) " (IGNITE_SKIP_PDF is set)" else
              " and this R has no aqua capability",
            ". The PNG above, figure_data/ and every number are unaffected.")
    return(invisible(NULL))
  }
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(out_pdf, plot, width = width_cm, height = height_cm, units = "cm",
                  bg = "white", device = .quartz_pdf)
  cat("wrote", .repo_rel(out_pdf), "\n")
}
