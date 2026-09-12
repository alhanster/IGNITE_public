#!/usr/bin/env Rscript
# See REPRODUCIBILITY.md, R (stage 2, `make figures`).

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()
repo_root <- normalizePath(file.path(script_dir, ".."))
req_file <- file.path(repo_root, "R-requirements.txt")

if (!file.exists(req_file)) {
  cat("  R-requirements.txt not found at", req_file, "\n")
  quit(status = 1)
}

# Only pkg==X.Y.Z lines are treated as version pins; comments and bare package names are advisory and skipped, matching _version_guard.py's ^pkg==version$ parsing.
raw <- readLines(req_file, warn = FALSE)
lines <- trimws(sub("#.*$", "", raw))
pins <- lines[grepl("^[A-Za-z0-9._-]+==[0-9]", lines)]

if (!length(pins)) {
  cat("  no exact == pins in R-requirements.txt; nothing enforced\n")
  quit(status = 0)
}

problems <- character(0)

for (p in pins) {
  parts <- strsplit(p, "==", fixed = TRUE)[[1]]
  pkg <- trimws(parts[1]); want <- trimws(parts[2])
  if (!requireNamespace(pkg, quietly = TRUE)) {
    problems <- c(problems, sprintf("  %-10s pinned %-8s NOT INSTALLED", pkg, want))
    next
  }
  have <- as.character(utils::packageVersion(pkg))
  if (have != want)
    problems <- c(problems, sprintf("  %-10s pinned %-8s installed %s", pkg, want, have))
}

# Reads `raw` here, not `lines`, since `lines` is comment-stripped and the pattern would not match.
want_r <- sub(".*R ([0-9.]+).*", "\\1", grep("^#[[:space:]]*R [0-9]", raw, value = TRUE)[1])
have_r <- as.character(getRversion())
if (!is.na(want_r) && nzchar(want_r) && want_r != have_r)
  problems <- c(problems, sprintf("  %-10s pinned %-8s installed %s", "R", want_r, have_r))

# Every renderer requests the font by name. Substitution is silent, shifting text layout without raising an error.
font_note <- NULL
if (requireNamespace("systemfonts", quietly = TRUE)) {
  # match_font() is deprecated since systemfonts 1.1.0 in favor of match_fonts(); call the new name with a fallback.
  resolved <- if ("match_fonts" %in% getNamespaceExports("systemfonts")) {
    systemfonts::match_fonts("Arial Unicode MS")$path[1]
  } else {
    systemfonts::match_font("Arial Unicode MS")$path
  }
  if (!grepl("Arial Unicode", basename(resolved), ignore.case = TRUE)) {
    msg <- sprintf("  %-10s 'Arial Unicode MS' resolves to %s -- glyph metrics will differ",
                   "font", basename(resolved))
    # The font ships with macOS, where its absence signals a broken install. Elsewhere it is
    # expected to be missing: PNG pixels then differ from the macOS baseline (a soft class in
    # verify-figures) while every number still reproduces, so it is reported, not fatal.
    if (identical(Sys.info()[["sysname"]], "Darwin")) {
      problems <- c(problems, msg)
    } else {
      font_note <- paste0(msg, "\n  (not fatal off macOS: PNGs will differ from the macOS baseline; ",
                          "numbers are unaffected)")
    }
  }
} else {
  font_note <- "  (systemfonts absent: could not check the 'Arial Unicode MS' requirement)"
}

if (length(problems)) {
  cat(strrep("=", 72), "\n")
  cat("R ENVIRONMENT MISMATCH -- rendered figures will NOT byte-match the baseline\n")
  cat(strrep("=", 72), "\n")
  cat(problems, sep = "\n"); cat("\n\n")
  cat("The figures are compared byte-for-byte against final_plots.sha256. Rendering\n")
  cat("libraries change pixel output between releases, so a mismatch here shows up as a\n")
  cat("figure regression with no other symptom. See R-requirements.txt for exact installs.\n")
  cat(strrep("=", 72), "\n")
  quit(status = 1)
}

cat("  R version guard: OK -- R", have_r, "and all pinned packages match\n")
if (!is.null(font_note)) cat(font_note, "\n")
quit(status = 0)
