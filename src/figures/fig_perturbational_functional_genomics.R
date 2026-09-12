#!/usr/bin/env Rscript
# fig_perturbational_functional_genomics.R
#
# Fig 3 layout: A / B / (C | D), a full-width row per evidence class followed by
# two condition-faceted cytokine panels side by side.
#
#   A  Vs-Th0 significance (-log10 adjusted p), Th1/Th2/Th17/Treg
#   B  Residualized trans-regulatory burden (expected_n_regulators_residuals)
#   C  Significantly regulated cytokines, by condition
#   D  Significantly regulated cytokine receptors

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

source(file.path(script_dir, "palette.R"))   # TARGET_COLORS (approved/in_trials/other)
source(file.path(script_dir, "vector_output.R"))  # save_figure(): the PNG + PDF pair

FIGDATA  <- file.path(repo_root, "figure_data")
OUT_DIR  <- file.path(repo_root, "final_plots")
OUT      <- file.path(OUT_DIR, "perturbational_functional_genomics_with_thde.png")
OUT_PDF  <- file.path(OUT_DIR, "pdf", "perturbational_functional_genomics_with_thde.pdf")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Shared style
font_family <- "Helvetica"
AXIS_TITLE <- 8
AXIS_TEXT  <- 7
STRIP_SIZE <- AXIS_TEXT
TAG_SIZE   <- 8
LW       <- function(pt) pt / .pt
STROKE   <- function(pt) pt / (.stroke / 2)
HAIRLINE <- 0.8

base_theme <- theme_classic(base_size = AXIS_TITLE) +
  theme(text        = element_text(family = font_family),
        axis.title  = element_text(size = AXIS_TITLE),
        axis.text   = element_text(size = AXIS_TEXT),
        axis.line   = element_line(linewidth = LW(HAIRLINE)),
        axis.ticks  = element_line(linewidth = LW(HAIRLINE)),
        plot.margin = margin(6, 10, 6, 10),
        legend.position = "none",
        plot.tag    = element_text(size = TAG_SIZE, face = "bold"))

# Group labels reflect counts computed by build_perturbational_fg_tables.py.
TARGET_LEVELS <- c("approved", "in_trials", "other")


make_target_factor <- function(df, n_tab, sep = "\n") {
  gv   <- function(g) as.integer(n_tab$n[n_tab$group == g])
  lab  <- c(approved  = sprintf("Approved%s(n = %s)",  sep, format(gv("approved"),  big.mark = ",")),
            in_trials = sprintf("In-trial%s(n = %s)", sep, format(gv("in_trials"), big.mark = ",")),
            other     = sprintf("Other%s(n = %s)",     sep, format(gv("other"),     big.mark = ",")))
  list(Target = factor(lab[df$group], levels = lab[TARGET_LEVELS]),
       fill   = setNames(unname(TARGET_COLORS[TARGET_LEVELS]), lab[TARGET_LEVELS]))
}

# Panel A: Th-subset vs-Th0 significance, faceted by target status
TH_LABELS   <- c("Th1", "Th2", "Th17", "Treg")

th_long <- read.csv(file.path(FIGDATA, "th_de_long.csv"), stringsAsFactors = FALSE)
th_long$.y     <- th_long$neglog10_adjp
th_long$subset <- factor(th_long$subset, levels = TH_LABELS)

th_ntab <- read.csv(file.path(FIGDATA, "th_de_group_n.csv"), stringsAsFactors = FALSE)
th_n    <- function(g) as.integer(th_ntab$n[th_ntab$group == g])
th_lab  <- c(approved  = sprintf("Approved (n = %s)",  format(th_n("approved"),  big.mark = ",")),
             in_trials = sprintf("In-trial (n = %s)", format(th_n("in_trials"), big.mark = ",")),
             other     = sprintf("Other (n = %s)",     format(th_n("other"),     big.mark = ",")))
th_long$Panel <- factor(th_lab[th_long$group], levels = th_lab[TARGET_LEVELS])
th_fill <- setNames(unname(TARGET_COLORS[TARGET_LEVELS]), TARGET_LEVELS)

panelA <- ggplot(th_long, aes(x = subset, y = .y, fill = group)) +
  geom_violin(trim = TRUE, colour = NA, alpha = 0.75, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
               colour = "grey30", linewidth = LW(HAIRLINE)) +
  facet_wrap(~Panel, nrow = 1) +
  coord_cartesian(ylim = c(0, 15)) +
  scale_fill_manual(values = th_fill) +
  labs(x = NULL, y = "vs-Th0 significance (-log10 adj p)", tag = "a") +
  base_theme +
  theme(axis.title.y     = element_text(size = AXIS_TITLE, margin = margin(r = 2)),
        strip.text       = element_text(size = STRIP_SIZE),
        strip.background = element_rect(fill = "grey92", colour = NA))

# Panel B: residualized trans-regulatory burden
b_data <- read.csv(file.path(FIGDATA, "regulator_residual.csv"), stringsAsFactors = FALSE)
b_data$.y <- b_data$residual
b_ntab <- read.csv(file.path(FIGDATA, "regulator_residual_group_n.csv"),
                   stringsAsFactors = FALSE)
b_tf <- make_target_factor(b_data, b_ntab); b_data$Target <- b_tf$Target

panelB <- ggplot(b_data, aes(x = Target, y = .y, fill = Target)) +
  geom_violin(trim = TRUE, colour = NA, alpha = 0.75, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", colour = "grey30",
               linewidth = LW(HAIRLINE)) +
  coord_cartesian(ylim = c(-30, 100)) +
  scale_fill_manual(values = b_tf$fill) +
  labs(x = NULL, y = "Residualized trans-regulatory burden", tag = "b") +
  base_theme +
  theme(axis.title.y = element_text(size = AXIS_TITLE, margin = margin(r = 2)))

# Shared builder for the two condition-faceted cytokine panels (C and D)
make_cyt_panel <- function(data_csv, n_csv, y_label, tag) {
  d <- read.csv(file.path(FIGDATA, data_csv), stringsAsFactors = FALSE)
  d$.y    <- d$n_sig

  d$panel <- factor(d$condition, levels = c("Rest", "Stim8hr", "Stim48hr"),
                    labels = c("Rest", "Stim 8hr", "Stim 48hr"))

  n_tab <- read.csv(file.path(FIGDATA, n_csv), stringsAsFactors = FALSE)

  tf <- make_target_factor(d, n_tab, sep = " "); d$Target <- tf$Target

  ggplot(d, aes(x = Target, y = .y, fill = Target)) +
    geom_violin(trim = TRUE, colour = NA, alpha = 0.75, scale = "width") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 1.4,
                 fill = "white", colour = "black", stroke = STROKE(HAIRLINE)) +
    facet_wrap(~panel, nrow = 1) +
    coord_cartesian(ylim = c(0, 5)) +
    scale_fill_manual(values = tf$fill) +
    labs(x = NULL, y = y_label, tag = tag) +
    base_theme +
    theme(axis.title.y     = element_text(size = AXIS_TITLE, margin = margin(r = 2)),
          axis.text.x      = element_text(size = AXIS_TEXT, angle = 45, hjust = 1),
          plot.margin      = margin(6, 10, 6, 18),
          strip.text       = element_text(size = STRIP_SIZE),
          strip.background = element_rect(fill = "grey92", colour = NA))
}

panelC <- make_cyt_panel("cytokine_counts.csv", "cytokine_group_n.csv",
                         "Sig. regulated cytokines\n(per gene)", "c")
panelD <- make_cyt_panel("cytokine_receptor_counts.csv", "cytokine_receptor_group_n.csv",
                         "Sig. regulated cytokine\nreceptors (per gene)", "d")

# Compose and save
bottom <- panelC | panelD
bottom_free <- if (utils::packageVersion("patchwork") >= "1.2.0") {
  patchwork::free(bottom, side = "l")
} else bottom
fig <- panelA / panelB / bottom_free + plot_layout(heights = c(1, 1, 1))

FIG_W_CM <- 18.0
FIG_H_CM <- 17.0

save_figure(fig, OUT, OUT_PDF, FIG_W_CM, FIG_H_CM)