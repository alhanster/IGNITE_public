#!/usr/bin/env Rscript
# fig_leakage_controlled.R
#
# Produces the leakage-controlled comparison against GPS and missense-z reported as Supplementary S3.

for (loc in c("en_US.UTF-8", "C.UTF-8", "UTF-8")) {
  if (suppressWarnings(Sys.setlocale("LC_ALL", loc)) != "") break
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(jsonlite)
})

args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir  <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()
repo_root   <- normalizePath(file.path(script_dir, "..", ".."))
stopifnot(file.exists(file.path(repo_root, "Makefile")))
source(file.path(script_dir, "vector_output.R"))
need <- function(f) {
  p <- file.path(repo_root, "figure_data", f)
  if (!file.exists(p))
    stop("Missing ", p, "\nRun `make tables` first.")
  p
}

COL_OURS  <- "#4477AA"
COL_OTHER <- "#C3B8AE"

font_family <- "Helvetica"
update_geom_defaults("text", list(family = font_family))

AXIS_TITLE <- 8
AXIS_TEXT  <- 7
TITLE_SIZE <- 8
DENSE      <- 6.5
GEOM_PT <- function(pt) pt / .pt
LW       <- function(pt) pt / .pt
HAIRLINE <- 0.8
HATCH_LW <- 0.6

base_theme <- theme_classic(base_size = AXIS_TITLE) +
  theme(text       = element_text(family = font_family),
        axis.title = element_text(size = AXIS_TITLE),
        axis.text  = element_text(size = AXIS_TEXT),
        axis.line  = element_line(linewidth = LW(HAIRLINE)),
        axis.ticks = element_line(linewidth = LW(HAIRLINE)))

# --- Data -------------------------------------------------------------------
cmp <- read.csv(need("leakage_controlled_comparison.csv"), stringsAsFactors = FALSE)
dl  <- jsonlite::fromJSON(need("leakage_controlled_delong.json"))

ORDER <- c("PU full model", "GPS overall", "missense-z")
XLAB  <- c("PU full model" = "IGNITE", "GPS overall" = "GPS", "missense-z" = "Missense Z-score")

if ("L2G" %in% cmp$method)
  stop("L2G left this analysis on 2026-08-31; a sheet still carrying it is a stale figure_data/")

std <- cmp[cmp$evaluation == "standard_taskA", ]
ctl <- cmp[cmp$evaluation == "leakage_controlled", ]
stopifnot(setequal(std$method, ORDER), setequal(ctl$method, ORDER))
n_pos_std <- unique(std$n_pos); n_pos_ctl <- unique(ctl$n_pos)
stopifnot(length(n_pos_std) == 1, length(n_pos_ctl) == 1)

BASE <- 0.50
W    <- 0.36

bars <- do.call(rbind, lapply(seq_along(ORDER), function(i) {
  m <- ORDER[i]
  data.frame(method = m, x = i,
             arm  = c("standard", "controlled"),
             dx   = c(-W / 2, W / 2),
             auc  = c(std$auc[match(m, std$method)], ctl$auc[match(m, ctl$method)]),
             fill = if (m == "PU full model") COL_OURS else COL_OTHER,
             stringsAsFactors = FALSE)
}))
bars$xmin <- bars$x + bars$dx - W / 2
bars$xmax <- bars$x + bars$dx + W / 2

# --- Geometry ---------------------------------------------------------------
XLIM <- c(0.5, 3.5)
YLIM <- c(BASE, 0.686)
FIG_W_CM <- 12.0
FIG_H_CM <- 7.8

PANEL_W_CM <- FIG_W_CM - 2.2
PANEL_H_CM <- FIG_H_CM - 1.9

HATCH_SLOPE <- (PANEL_W_CM / PANEL_H_CM) * (diff(YLIM) / diff(XLIM))
HATCH_GAP_CM <- 0.105
HATCH_DB     <- HATCH_GAP_CM / ((PANEL_H_CM / diff(YLIM)) * cos(pi / 4))

hatch_segments <- function(xmin, xmax, ymin, ymax, slope = HATCH_SLOPE, db = HATCH_DB) {
  bs   <- c(ymin, ymax) - slope * rep(c(xmin, xmax), each = 2)
  bseq <- seq(floor(min(bs) / db) * db, ceiling(max(bs) / db) * db, by = db)
  segs <- lapply(bseq, function(b) {
    xa <- (ymin - b) / slope
    xb <- (ymax - b) / slope
    lo <- max(xmin, min(xa, xb))
    hi <- min(xmax, max(xa, xb))
    if (hi <= lo) return(NULL)
    data.frame(x = lo, xend = hi, y = slope * lo + b, yend = slope * hi + b)
  })
  do.call(rbind, segs)
}

hatched <- bars[bars$arm == "controlled", ]
hatch <- do.call(rbind, lapply(seq_len(nrow(hatched)), function(i)
  hatch_segments(hatched$xmin[i], hatched$xmax[i], BASE, hatched$auc[i])))

# --- Significance brackets ---------------------------------------------------
IG   <- which(ORDER == "PU full model")
GPS_ <- which(ORDER == "GPS overall")
br <- data.frame(
  dx    = c(-W / 2, W / 2),
  y     = c(0.660, 0.672),
  stars = c(dl$arms$standard_taskA$stars, dl$arms$leakage_controlled$stars),
  stringsAsFactors = FALSE)
br$x1 <- IG + br$dx
br$x2 <- GPS_ + br$dx
BT <- 0.0022

# --- Hand-drawn legend --------------------------------------------------------
LEG_SW_X <- 2.36
LEG_SW_W <- 0.14
LEG_Y    <- c(0.6715, 0.6465)
LEG_H    <- 0.0055
leg <- data.frame(
  xmin = LEG_SW_X, xmax = LEG_SW_X + LEG_SW_W,
  ymin = LEG_Y - LEG_H, ymax = LEG_Y + LEG_H,
  label = c(sprintf("Standard\n(all %d trial targets)", n_pos_std),
            sprintf("Leakage-controlled\n(%d GPS-blind targets)", n_pos_ctl)),
  stringsAsFactors = FALSE)
leg_hatch <- hatch_segments(leg$xmin[2], leg$xmax[2], leg$ymin[2], leg$ymax[2])

p <- ggplot() +
  geom_rect(data = bars, aes(xmin = xmin, xmax = xmax, ymin = BASE, ymax = auc, fill = fill),
            colour = "black", linewidth = LW(HAIRLINE)) +
  geom_segment(data = hatch, aes(x = x, xend = xend, y = y, yend = yend),
               colour = "black", linewidth = LW(HATCH_LW)) +
  geom_rect(data = bars, aes(xmin = xmin, xmax = xmax, ymin = BASE, ymax = auc),
            fill = NA, colour = "black", linewidth = LW(HAIRLINE)) +
  geom_text(data = bars, aes(x = (xmin + xmax) / 2, y = auc, label = sprintf("%.3f", auc)),
            vjust = -0.55, size = GEOM_PT(DENSE)) +
  geom_segment(data = br, aes(x = x1, xend = x2, y = y, yend = y), linewidth = LW(HAIRLINE)) +
  geom_segment(data = br, aes(x = x1, xend = x1, y = y - BT, yend = y), linewidth = LW(HAIRLINE)) +
  geom_segment(data = br, aes(x = x2, xend = x2, y = y - BT, yend = y), linewidth = LW(HAIRLINE)) +
  geom_text(data = br, aes(x = (x1 + x2) / 2, y = y + 0.0009, label = stars),
            vjust = 0, size = GEOM_PT(TITLE_SIZE)) +
  geom_rect(data = leg, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "#cccccc", colour = "black", linewidth = LW(HAIRLINE)) +
  geom_segment(data = leg_hatch, aes(x = x, xend = xend, y = y, yend = yend),
               colour = "black", linewidth = LW(HATCH_LW)) +
  geom_text(data = leg, aes(x = xmax + 0.06, y = (ymin + ymax) / 2, label = label),
            hjust = 0, size = GEOM_PT(DENSE), lineheight = 0.95) +
  scale_fill_identity() +
  scale_x_continuous(breaks = seq_along(ORDER), labels = unname(XLAB[ORDER])) +
  scale_y_continuous(breaks = seq(0.500, 0.675, 0.025)) +
  coord_cartesian(xlim = XLIM, ylim = YLIM, clip = "off") +
  labs(x = NULL, y = "AUC (recovering held-out trial targets)") +
  base_theme +
  theme(plot.margin = margin(6, 6, 6, 6))

OUT     <- file.path(repo_root, "final_plots", "supplementary",
                     "leakage_controlled_comparison.png")
OUT_PDF <- file.path(repo_root, "final_plots", "pdf",
                     "leakage_controlled_comparison.pdf")
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)

save_figure(p, OUT, OUT_PDF, FIG_W_CM, FIG_H_CM)

cat(sprintf("   %d bars | brackets %s / %s | n_pos %d/%d\n",
            nrow(bars), br$stars[1], br$stars[2], n_pos_std, n_pos_ctl))
