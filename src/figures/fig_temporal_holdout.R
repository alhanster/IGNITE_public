#!/usr/bin/env Rscript
# fig_temporal_holdout.R
#
# Three-panel prospective temporal holdout (Fig 6): design timeline, paired
# half-violin comparison of rank percentile at T0, and T0 sensitivity sweep
# across freeze years.
#
# Inputs are written by build_temporal_holdout.py to figure_data/.
# See REPRODUCIBILITY.md, Temporal holdout (T0) labels, for the split and file-naming details.

for (loc in c("en_US.UTF-8", "C.UTF-8", "UTF-8")) {
  if (suppressWarnings(Sys.setlocale("LC_ALL", loc)) != "") break
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
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

# --- Palette (semantic) ---
C_FULL <- "#1f5fa8"
C_EMG  <- "#b5341a"
C_UNL  <- "#9aa3ad"
C_GEN  <- "#E69F00"

stars <- function(p) {
  ifelse(p < 1e-4, "****",
  ifelse(p < 1e-3, "***",
  ifelse(p < 1e-2, "**",
  ifelse(p < 0.05, "*", "ns"))))
}

font_family <- "Helvetica"
update_geom_defaults("text",  list(family = font_family))
update_geom_defaults("label", list(family = font_family))

AXIS_TITLE <- 8
AXIS_TEXT  <- 7
TITLE_SIZE <- 8
TAG_SIZE   <- 8
DENSE      <- 6.5
TINY       <- 6
GEOM_PT <- function(pt) pt / .pt
LW      <- function(pt) pt / .pt
HAIRLINE <- 0.8
RULE     <- 1.2

base_theme <- theme_classic(base_size = AXIS_TITLE) +
  theme(text        = element_text(family = font_family),
        axis.title  = element_text(size = AXIS_TITLE),
        axis.text   = element_text(size = AXIS_TEXT),
        axis.line   = element_line(linewidth = LW(HAIRLINE)),
        axis.ticks  = element_line(linewidth = LW(HAIRLINE)),
        plot.title  = element_text(size = TITLE_SIZE, hjust = 0),
        plot.tag    = element_text(size = TAG_SIZE, face = "bold"),
        plot.margin = margin(6, 10, 6, 10))

# --- Load data ---
split  <- jsonlite::fromJSON(need("t0_split.json"))
T0     <- as.integer(split$T0)
scores <- read.csv(need(sprintf("scores_at_T0_%d.csv", T0)), stringsAsFactors = FALSE)
aucd   <- read.csv(need("auc_emergent.csv"), stringsAsFactors = FALSE)
sens   <- read.csv(need("t0_sensitivity.csv"), stringsAsFactors = FALSE)

pstats <- read.csv(need("t0_panel_stats.csv"), stringsAsFactors = FALSE)
nP     <- split$n_pos_at_T0
nE     <- split$n_emergent
nKnown <- split$n_known_nonapproved_at_T0
nDrop  <- length(split$dropped_undated)
nHeld  <- nKnown + nDrop

auc_full <- aucd$auc[aucd$model == "full"]
auc_gen  <- aucd$auc[aucd$model == "genetics"]

# ==== Panel a: design timeline ====
tl <- data.frame(
  y      = c(3, 2, 1),
  x0     = c(1998, T0, 1998),
  x1     = c(T0, 2026, 2026),
  fill   = c(C_FULL, C_EMG, C_UNL),
  alpha  = c(0.9, 0.9, 0.5),
  lab    = c(sprintf("train: %d\napproved-by-%d", nP, T0),
             sprintf("test: %d\nemergent\n(first trial\n> %d)", nE, T0),
             sprintf("%d labelled targets held out of negatives\n(%d already in clinic, %d undatable)",
                     nHeld, nKnown, nDrop)),
  labx   = c((1998 + T0) / 2, (T0 + 2026) / 2, (1998 + 2026) / 2),
  labcol = c("white", "white", "#333333"),
  labsz  = c(GEOM_PT(TINY), GEOM_PT(TINY), GEOM_PT(TINY))
)
panelA <- ggplot(tl) +
  geom_rect(aes(xmin = x0, xmax = x1, ymin = y - 0.32, ymax = y + 0.32,
                fill = I(fill), alpha = I(alpha))) +
  geom_text(aes(x = labx, y = y, label = lab, colour = I(labcol), size = I(labsz)),
            lineheight = 0.9) +
  geom_vline(xintercept = T0, linetype = "dashed", linewidth = LW(HAIRLINE)) +
  annotate("text", x = T0, y = 3.6, label = sprintf("T0 = %d", T0),
           hjust = 0.5, vjust = 0, size = GEOM_PT(TINY), family = font_family) +
  scale_x_continuous(breaks = seq(2000, 2025, 5), name = "Year") +
  scale_y_continuous() +
  coord_cartesian(xlim = c(1998, 2027), ylim = c(0.5, 3.85)) +
  labs(title = sprintf("Freeze labels at %d and\ntest post-%d trial entries", T0, T0)) +
  base_theme +
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.y = element_blank())

# ==== Panel b: paired half-violins of rank percentile at T0, one per model ====
pool <- scores[scores$in_pool %in% c(TRUE, "True", "true"), ]
pool$is_emergent <- pool$is_emergent %in% c(TRUE, "True", "true")

srow <- sens[sens$T0 == T0, ]
stopifnot(nrow(srow) == 1,
          abs(srow$auc_full - auc_full) < 5e-4, abs(srow$auc_genetics - auc_gen) < 5e-4)
nHigher <- pstats$n_full_higher
nPaired <- pstats$n_emergent_paired

HALFW <- 0.34
half_violin <- function(x, xc, side, halfw = HALFW, adjust = 1.15) {
  d <- stats::density(x, from = 0, to = 1, adjust = adjust)
  w <- d$y / max(d$y) * halfw
  data.frame(x = c(xc + side * w, rep(xc, length(w))), y = c(d$x, rev(d$x)))
}

VB <- list(list(xc = 1, se = -1, col = C_FULL, e = "rank_pctile_full",
                lab = "Full Model",
                med = pstats$median_emergent,
                q1  = pstats$q1_emergent,          q3 = pstats$q3_emergent),
           list(xc = 2, se = +1, col = C_GEN,  e = "rank_pctile_genetics",
                lab = "Genetics-Only",
                med = pstats$median_emergent_genetics,
                q1  = pstats$q1_emergent_genetics, q3 = pstats$q3_emergent_genetics))

polys <- do.call(rbind, lapply(seq_along(VB), function(i) {
  v  <- VB[[i]]
  ev <- pool[[v$e]][pool$is_emergent]; ev <- ev[!is.na(ev)]
  transform(half_violin(ev, v$xc, v$se), grp = paste0("e", i), fill = v$col)
}))
boxes <- do.call(rbind, lapply(VB, function(v) data.frame(
  x = v$xc + v$se * 0.07, q1 = v$q1, q3 = v$q3, med = v$med,
  labx = v$xc + v$se * 0.13, hj = if (v$se < 0) 1 else 0, lab = sprintf("%.2f", v$med))))

pairs_df <- pool[pool$is_emergent, c("rank_pctile_full", "rank_pctile_genetics")]
pairs_df <- pairs_df[stats::complete.cases(pairs_df), ]
stopifnot(nrow(pairs_df) == nPaired,
          sum(pairs_df$rank_pctile_full > pairs_df$rank_pctile_genetics) == nHigher)
LX0 <- 1 + 0.12; LX1 <- 2 - 0.12

YB <- 1.10
ann_b <- stars(srow$delong_p_two_sided)

panelB <- ggplot() +
  geom_hline(yintercept = pstats$median_unlabeled, colour = C_UNL, linetype = "dotted",
             linewidth = LW(HAIRLINE)) +
  annotate("text", x = 2.50, y = pstats$median_unlabeled + 0.03, label = "pool\nmedian",
           hjust = 0.5, vjust = 0, size = GEOM_PT(TINY), colour = C_UNL, lineheight = 0.9,
           family = font_family) +
  geom_segment(data = pairs_df,
               aes(x = LX0, xend = LX1, y = rank_pctile_full, yend = rank_pctile_genetics),
               colour = C_UNL, linewidth = LW(0.5), alpha = 0.22) +
  geom_polygon(data = polys, aes(x = x, y = y, group = grp, fill = I(fill)), colour = NA,
               alpha = 0.88) +
  geom_segment(data = boxes, aes(x = x, xend = x, y = q1, yend = q3),
               colour = "#111111", linewidth = LW(RULE), lineend = "butt") +
  geom_point(data = boxes, aes(x = x, y = med), shape = 21, fill = "white",
             colour = "#111111", size = 1.6, stroke = LW(HAIRLINE)) +
  geom_text(data = boxes, aes(x = labx, y = med, label = lab, hjust = hj),
            size = GEOM_PT(DENSE), family = font_family) +
  annotate("segment", x = 1, xend = 2, y = YB, yend = YB, colour = "#111111",
           linewidth = LW(HAIRLINE)) +
  annotate("segment", x = c(1, 2), xend = c(1, 2), y = YB, yend = YB - 0.025,
           colour = "#111111", linewidth = LW(HAIRLINE)) +
  annotate("text", x = 1.5, y = YB + 0.045, label = ann_b, size = GEOM_PT(TITLE_SIZE),
           family = font_family) +
  scale_x_continuous(breaks = c(1, 2), labels = sapply(VB, `[[`, "lab")) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0")) +
  coord_cartesian(xlim = c(0.55, 2.62), ylim = c(-0.02, YB + 0.12)) +
  labs(x = NULL, y = sprintf("Rank Percentile at T0 = %d", T0),
       title = sprintf("Emergent genes rank higher under\nfull-model for %d of %d genes",
                       nHigher, nPaired)) +
  base_theme +
  theme(axis.text.x  = element_text(size = AXIS_TEXT, margin = margin(t = 6)),
        axis.ticks.x = element_blank())

# ==== Panel c: T0 sensitivity ====
sens <- sens[order(sens$T0), ]
sens$stars <- stars(sens$delong_p_two_sided)

sens_long <- rbind(
  data.frame(T0 = sens$T0, auc = sens$auc_full,     n_emergent = sens$n_emergent, model = "full"),
  data.frame(T0 = sens$T0, auc = sens$auc_genetics, n_emergent = sens$n_emergent, model = "genetics")
)
sens_long$model <- factor(sens_long$model, levels = c("full", "genetics"))
model_cols <- c(full = C_FULL, genetics = C_GEN)

panelC <- ggplot(sens_long, aes(x = T0, y = auc, colour = model)) +
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = "#888888", linewidth = LW(HAIRLINE)) +
  annotate("text", x = min(sens$T0), y = 0.508, label = "Chance",
           hjust = 0, size = GEOM_PT(TINY), colour = "#666666") +
  geom_line(linewidth = LW(RULE)) +
  geom_point(size = 2.2) +
  geom_text(data = sens, aes(x = T0, y = auc_full, label = sprintf("n=%d", n_emergent)),
            inherit.aes = FALSE, vjust = -1.9, size = GEOM_PT(TINY), colour = "#333333") +
  geom_text(data = sens, aes(x = T0, y = pmax(auc_full, auc_genetics) + 0.05, label = stars),
            inherit.aes = FALSE, size = GEOM_PT(TITLE_SIZE), colour = "#111111",
            family = font_family) +
  scale_colour_manual(values = model_cols,
                       labels = c(full = "Full Model", genetics = "Genetics-Only"),
                       name = NULL) +
  scale_x_continuous(breaks = sens$T0) +
  coord_cartesian(xlim = c(min(sens$T0) - 0.4, max(sens$T0) + 0.4), ylim = c(0.45, 0.80)) +
  labs(title = "Signal is stable across\nthe freeze year",
       x = "T0 (Freeze Year)", y = "AUC (Emergent vs Unlabeled)") +
  base_theme +
  theme(legend.position = c(0.74, 0.35),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = AXIS_TEXT),
        legend.margin = margin(0, 0, 0, 0))

# ==== Assemble ====
fig <- (panelA | panelB | panelC) +
  plot_layout(widths = c(1.05, 1.0, 1.0)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = TAG_SIZE, family = font_family, face = "bold"))

OUT     <- file.path(repo_root, "final_plots", "figure_prospective_temporal_holdout.png")
OUT_PDF <- file.path(repo_root, "final_plots", "pdf", "figure_prospective_temporal_holdout.pdf")
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)

FIG_W_CM <- 18.0
FIG_H_CM <- 7.8

save_figure(fig, OUT, OUT_PDF, FIG_W_CM, FIG_H_CM)
