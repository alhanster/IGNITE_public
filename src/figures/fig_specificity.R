#!/usr/bin/env Rscript
# fig_specificity.R
#
# Immune-vs-cardiac specificity control (Fig 7): one panel per recovery target-set,
# with DeLong paired-test significance brackets comparing the model against each OT score.
#
# Inputs are read from figure_data/, produced by make tables, which must run before make figures.

for (loc in c("en_US.UTF-8", "C.UTF-8", "UTF-8")) {
  if (suppressWarnings(Sys.setlocale("LC_ALL", loc)) != "") break
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir  <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()
repo_root   <- normalizePath(file.path(script_dir, "..", ".."))
stopifnot(file.exists(file.path(repo_root, "Makefile")))
source(file.path(script_dir, "vector_output.R"))
FIGDATA     <- file.path(repo_root, "figure_data")

FIGURE_DPI     <- 300

# Shared style
font_family <- "Helvetica"
update_geom_defaults("text", list(family = font_family))

AXIS_TITLE <- 7
AXIS_TEXT  <- 6
TINY       <- 6
TITLE_SIZE <- 8
TAG_SIZE   <- 8
DENSE      <- 6
STRIP_SIZE <- 7
GEOM_PT <- function(pt) pt / .pt
LW      <- function(pt) pt / .pt
HAIRLINE <- 0.8
RULE     <- 1.2

base_theme <- theme_classic(base_size = AXIS_TITLE) +
  theme(text        = element_text(family = font_family),
        axis.title  = element_text(size = AXIS_TITLE),
        axis.text   = element_text(size = AXIS_TEXT, colour = "black"),
        axis.line   = element_line(linewidth = LW(HAIRLINE)),
        axis.ticks  = element_line(linewidth = LW(HAIRLINE)),
        plot.margin = margin(6, 10, 6, 10),
        legend.position = "none",
        plot.tag    = element_text(size = TAG_SIZE, face = "bold"),
        strip.text       = element_text(size = STRIP_SIZE),
        strip.background = element_blank())

RECOVER_COLORS <- c(immune_exclusive = "#E15759", cardiac_exclusive = "#4E79A7")

SCORE_ORDER  <- c("our_model", "OT_genetic_immune", "OT_genetic_cardiac")
SCORE_LABELS <- c(our_model          = "IGNITE",
                  OT_genetic_immune  = "Immune-\nspecific\nGenetic\nScore",
                  OT_genetic_cardiac = "Cardiac-\nspecific\nGenetic\nScore")
XPOS <- setNames(seq_along(SCORE_ORDER), SCORE_ORDER)

mat_path   <- file.path(FIGDATA, "specificity_matrix.csv")
delong_path <- file.path(FIGDATA, "delong_specificity_tests.csv")
if (!file.exists(mat_path))    stop("Missing ", mat_path, "\nRun `make tables` first.")
if (!file.exists(delong_path)) stop("Missing ", delong_path, "\nRun `make tables` first.")

m <- read.csv(mat_path, stringsAsFactors = FALSE)
m$score   <- factor(m$score,   levels = SCORE_ORDER)
m$recover <- factor(m$recover, levels = names(RECOVER_COLORS))

dl <- read.csv(delong_path, stringsAsFactors = FALSE)

n_by <- read.csv(file.path(FIGDATA, "specificity_panel_n.csv"), stringsAsFactors = FALSE)
names(n_by)[names(n_by) == "n_pos"] <- "n"
facet_labs <- setNames(
  c(sprintf("Recovering immune-exclusive\ntargets (n=%d)",
            n_by$n[n_by$recover == "immune_exclusive"]),
    sprintf("Recovering cardiac-exclusive\ntargets (n=%d)",
            n_by$n[n_by$recover == "cardiac_exclusive"])),
  c("immune_exclusive", "cardiac_exclusive"))

cihi <- setNames(m$ci_hi, paste(m$recover, m$score))
bar_pairs <- list(
  list(B = "OT_genetic_immune",  bump = 0.018),
  list(B = "OT_genetic_cardiac", bump = 0.030))

brk <- do.call(rbind, lapply(c("immune_exclusive", "cardiac_exclusive"), function(rc) {
  d <- dl[dl$panel == rc & dl$A == "our_model", ]
  base_top <- max(cihi[paste(rc, c("our_model", "OT_genetic_immune"))])
  lo_y <- base_top + 0.035
  hi_y <- max(cihi[paste(rc, c("our_model", "OT_genetic_cardiac"))] + 0.035, lo_y + 0.032)
  data.frame(
    recover = rc,
    x    = c(XPOS["our_model"], XPOS["our_model"]),
    xend = c(XPOS["OT_genetic_immune"], XPOS["OT_genetic_cardiac"]),
    y    = c(lo_y, hi_y),
    label = c(d$sig[d$B == "OT_genetic_immune"], d$sig[d$B == "OT_genetic_cardiac"]),
    stringsAsFactors = FALSE)
}))
brk$recover <- factor(brk$recover, levels = names(RECOVER_COLORS))
brk$ytick   <- brk$y - 0.008
brk$xmid    <- (brk$x + brk$xend) / 2
brk$lab_y   <- brk$y + 0.006

make_panel <- function(rc, tag, show_y) {
  md <- droplevels(m[m$recover == rc, ])
  bd <- brk[brk$recover == rc, ]
  g <- ggplot(md, aes(x = score, y = auc)) +
    geom_col(width = 0.66, fill = RECOVER_COLORS[[rc]],
             colour = "white", linewidth = LW(HAIRLINE)) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.16,
                  linewidth = LW(RULE), colour = "#333333") +
    geom_text(aes(y = ci_hi, label = sprintf("%.2f", auc)),
              vjust = -0.6, size = GEOM_PT(TITLE_SIZE)) +
    geom_segment(data = bd, aes(x = x, xend = xend, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = LW(HAIRLINE)) +
    geom_segment(data = bd, aes(x = x, xend = x, y = y, yend = ytick),
                 inherit.aes = FALSE, linewidth = LW(HAIRLINE)) +
    geom_segment(data = bd, aes(x = xend, xend = xend, y = y, yend = ytick),
                 inherit.aes = FALSE, linewidth = LW(HAIRLINE)) +
    geom_text(data = bd, aes(x = xmid, y = lab_y, label = label),
              inherit.aes = FALSE, size = GEOM_PT(DENSE), vjust = 0) +
    facet_wrap(~recover, labeller = as_labeller(facet_labs)) +
    scale_x_discrete(labels = SCORE_LABELS) +
    coord_cartesian(ylim = c(0.5, 0.82)) +
    labs(x = NULL, tag = tag,
         y = if (show_y) "AUC (recovering held-out trial targets)" else NULL) +
    base_theme +
    theme(axis.title.y = element_text(size = AXIS_TITLE, margin = margin(r = 2)),
          axis.text.x  = element_text(size = TINY))
  if (!show_y)
    g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  g
}

pA <- make_panel("immune_exclusive",  "a", show_y = TRUE)
pB <- make_panel("cardiac_exclusive", "b", show_y = FALSE)

p <- pA + pB

fname     <- "specificity_immune_vs_cardiac_faceted.png"
final_dir <- file.path(repo_root, "final_plots")
dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)

OUT     <- file.path(final_dir, fname)
OUT_PDF <- file.path(repo_root, "final_plots", "pdf", "specificity_immune_vs_cardiac_faceted.pdf")

FIG_W_CM <- 12.0
FIG_H_CM <- 7.8

save_figure(p, OUT, OUT_PDF, FIG_W_CM, FIG_H_CM)
