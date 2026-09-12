#!/usr/bin/env Rscript
# fig_stat4_vignette.R
#
# Builds a four-panel STAT4 case-study vignette as native ggplot2/patchwork objects,
# not a static image.

for (loc in c("en_US.UTF-8", "C.UTF-8", "UTF-8")) {
  if (suppressWarnings(Sys.setlocale("LC_ALL", loc)) != "") break
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(jsonlite)
  library(ggtext)
})

# --- repo-root anchor (script lives in src/figures/) ------------------------
args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir  <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()
repo_root   <- normalizePath(file.path(script_dir, "..", ".."))
stopifnot(file.exists(file.path(repo_root, "Makefile")))
source(file.path(script_dir, "vector_output.R"))

DATA_DIR <- file.path(repo_root, "figure_data")
OUT_DIR  <- file.path(repo_root, "final_plots")
OUT      <- file.path(OUT_DIR, "figure_stat4_vignette.png")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

C_GEN   <- "#E69F00"; C_OBS <- "#0072B2"; C_PERTURBATIONAL <- "#CC79A7"; C_FOCAL <- "#c0392b"
C_OTH   <- "#9a9a9a"
PATH_COL <- c(cd3 = "#e07b39", tcrk = C_FOCAL, costim = "#d9a520", calci = "#7d2e8d",
              other = C_OTH)

font_family <- "Helvetica"
update_geom_defaults("text", list(family = font_family))
PT           <- 1 / 2.845
TAG_SIZE      <- 8
TITLE_SIZE <- 7
AXIS_TITLE <- 6
AXIS_TEXT  <- 6
GT      <- AXIS_TITLE * PT
LW       <- function(pt) pt / .pt
HAIRLINE <- 0.8
RULE     <- 1.2
base_theme <- theme_classic(base_size = AXIS_TITLE) +
  theme(text        = element_text(family = font_family, colour = "#222222"),
        axis.title  = element_text(size = AXIS_TITLE, colour = "#222222"),
        axis.text   = element_text(size = AXIS_TEXT, colour = "#555555"),
        axis.ticks  = element_line(colour = "#999999", linewidth = LW(HAIRLINE)),
        axis.line   = element_line(colour = "#666666", linewidth = LW(HAIRLINE)),
        plot.title  = element_text(size = TITLE_SIZE, hjust = 0, face = "plain",
                                   colour = "#111111", margin = margin(b = 5)),
        plot.margin = margin(8, 14, 8, 14),
        legend.position = "none",
        plot.tag    = element_text(size = TAG_SIZE, face = "bold"))

meta <- fromJSON(file.path(DATA_DIR, "vignette_meta.json"))

# ============ Panel a: PU-score distribution with STAT4 marked ==============
scores <- read.csv(file.path(DATA_DIR, "vignette_score_hist.csv"),
                    stringsAsFactors = FALSE)
lab_a2 <- sprintf("atop(italic('STAT4') * ': rank %d / %s', '(top %.2f%%)')",
                   meta$stat4_rank, format(meta$n_genes, big.mark = ","), meta$stat4_top_pct)
panelA <- ggplot(scores, aes(x = pu_score)) +
  geom_histogram(bins = 60, fill = "#dcdcdc", colour = "white", linewidth = LW(HAIRLINE)) +
  geom_vline(xintercept = meta$approved_median, colour = "#555555",
             linetype = "dashed", linewidth = LW(HAIRLINE)) +
  geom_vline(xintercept = meta$stat4_score, colour = C_FOCAL, linewidth = LW(RULE)) +
  annotate("text", x = meta$approved_median - 0.02, y = Inf,
           label = "Approved Median", hjust = 1, vjust = 1.15, size = GT,
           lineheight = 0.9, colour = "#555555") +
  annotate("segment", x = meta$stat4_score - 0.30, xend = meta$stat4_score - 0.01,
           y = 470, yend = 340, colour = C_FOCAL, linewidth = LW(HAIRLINE),
           arrow = arrow(length = unit(0.14, "cm"), type = "closed")) +
  annotate("text", x = meta$stat4_score - 0.31, y = 485,
           label = lab_a2, parse = TRUE,
           hjust = 1, vjust = 0.5, size = GT, lineheight = 0.95, colour = C_FOCAL) +
  scale_x_continuous(breaks = seq(0, 1, 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
  labs(x = "PU-model target score", y = "Genes",
       title = expression(italic("STAT4")*" ranks among the very top")) +
  base_theme +
  theme(axis.title.y = element_text(margin = margin(r = 0)))

# ============ Panel b: four evidence axes ===================================
pb <- read.csv(file.path(DATA_DIR, "vignette_panelB_axes.csv"), stringsAsFactors = FALSE,
               check.names = FALSE)
pb <- pb[order(pb$order), ]
pb$axis[pb$axis == "Trans-regulatory convergence\n(perturbational)"] <-
  "Trans-regulatory burden (perturbational)"
pb$axis <- factor(pb$axis, levels = rev(pb$axis))
pb$col  <- c(genetic = C_GEN, observational = C_OBS, perturbational = C_PERTURBATIONAL)[pb$evidence]
panelB <- ggplot(pb, aes(x = percentile, y = axis, fill = I(col))) +
  geom_col(width = 0.66) +
  geom_vline(xintercept = 90, colour = "#888888", linetype = "dotted", linewidth = LW(HAIRLINE)) +
  geom_text(aes(label = sprintf("%.0f", percentile)), hjust = 1.4,
            colour = "white", size = GT) +
  annotate("text", x = 90, y = length(levels(pb$axis)) + 0.55, label = "90th pctile",
           size = GT, colour = "#888888", vjust = 0) +
  scale_x_continuous(breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0.02))) +
  scale_y_discrete(expand = expansion(add = c(0.6, 0.95))) +
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = "Percentile among measured genes", y = NULL,
       title = "STAT4 percentiles relative to all measured genes") +
  base_theme + theme(axis.text.y = element_text(size = AXIS_TEXT, colour = "#333333",
                                                lineheight = 0.9))

# ============ Panel c: composition of the immune-GWAS feature ==============
pc   <- read.csv(file.path(DATA_DIR, "vignette_panelC_gwas_diseases.csv"), stringsAsFactors = FALSE)
mc   <- fromJSON(file.path(DATA_DIR, "vignette_panelC_gwas_meta.json"))
pc   <- pc[pc$contributes_to_gwas_score == 1, ]
pc   <- pc[order(-pc$gwas_credible_sets_score), ]
pc$label <- factor(pc$label, levels = rev(pc$label))
stopifnot(nrow(pc) == mc$n_contributing_diseases,
          sum(pc$n_credible_sets) == mc$n_credible_sets_total,
          max(pc$gwas_credible_sets_score) < mc$gwas_score_feature_value)
feat_c <- mc$gwas_score_feature_value
XMAX_C <- 0.84
XMIN_C <- 0.25
title_c <- sprintf(
  "Immune GWAS score of *STAT4*:<br>%d credible sets across %d immune disease terms",
  mc$n_credible_sets_total, mc$n_contributing_diseases)
panelC <- ggplot(pc, aes(x = gwas_credible_sets_score, y = label)) +
  geom_segment(aes(x = 0, xend = gwas_credible_sets_score, yend = label),
               colour = "#d0d0d0", linewidth = LW(HAIRLINE)) +
  geom_vline(xintercept = feat_c, linetype = "dashed", colour = "#666666", linewidth = LW(HAIRLINE)) +
  geom_point(colour = C_GEN, size = 1.9) +
  geom_text(aes(label = n_credible_sets), hjust = 0, nudge_x = 0.0155,
            size = AXIS_TEXT * PT, colour = "#777777") +
  annotate("text", x = feat_c - 0.022, y = nrow(pc) + 1.35, hjust = 1, vjust = 1,
           size = GT, colour = "#555555",
           label = sprintf("model feature = %.3f", feat_c)) +
  scale_x_continuous(breaks = seq(0, 0.75, 0.25), expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.5))) +
  coord_cartesian(xlim = c(XMIN_C, XMAX_C)) +
  labs(x = "Open Targets GWAS credible-set score", y = NULL, title = title_c) +
  base_theme +
  theme(axis.text.y = element_text(size = AXIS_TEXT, colour = "#333333", hjust = 1),
        plot.title = element_markdown())

# ============ Panel d: kNN nearest approved-target neighbors ==============
pd <- read.csv(file.path(DATA_DIR, "vignette_panelD_knn.csv"), stringsAsFactors = FALSE,
               check.names = FALSE)
pd <- pd[order(pd$order), ]
pd$y   <- rev(seq_len(nrow(pd)))
pd$col <- PATH_COL[pd$pathway]
title_d <- "Nearest neighbors of *STAT4* in knockdown-signature space"
x_d     <- "Cosine similarity of knockdown signature to *STAT4* (Stim 48hr)"
panelD <- ggplot(pd) +
  geom_segment(aes(x = 0, xend = cosine, y = y, yend = y), colour = "#d0d0d0",
               linewidth = LW(HAIRLINE)) +
  geom_vline(xintercept = meta$knn_p95, colour = "#999999", linetype = "dashed",
             linewidth = LW(HAIRLINE)) +
  annotate("text", x = meta$knn_p95 - 0.006, y = max(pd$y) + 0.72, label = "95th pctile",
           size = GT, colour = "#888888", vjust = 1, hjust = 1) +
  geom_point(aes(x = cosine, y = y, colour = I(col)), size = 2.6) +
  geom_point(aes(x = cosine, y = y), shape = 21, colour = "white", fill = NA,
             size = 2.6, stroke = LW(HAIRLINE)) +
  geom_text(aes(x = cosine + 0.014, y = y + 0.24, label = gene), hjust = 0,
            size = GT) +
  geom_text(aes(x = cosine + 0.014, y = y - 0.02, label = mechanism, colour = I(col)),
            hjust = 0, fontface = "italic", size = GT, family = "Arial Unicode MS") +
  geom_text(aes(x = cosine + 0.014, y = y - 0.28, label = paste0("\u2192 ", approved_drugs_label)),
            hjust = 0, colour = "#666666", size = GT, family = "Arial Unicode MS") +
  scale_x_continuous(breaks = seq(0, 0.6, 0.2), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous() +
  coord_cartesian(xlim = c(0, 0.85), ylim = c(0.35, max(pd$y) + 0.95)) +
  labs(x = x_d, y = NULL, title = title_d) +
  base_theme +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.title = element_markdown(), axis.title.x = element_markdown(),
        plot.tag.position = c(0.24, 1))

# --- assemble 2x2 panel grid -------------------------------------------
free_l <- function(p) {
  if (utils::packageVersion("patchwork") >= "1.2.0") patchwork::free(p, side = "l") else p
}
fig <- free_l(panelA | panelB) / free_l(panelC | panelD) +
  plot_layout(heights = c(0.38, 0.62)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = TAG_SIZE, family = font_family, face = "bold"))

OUT_PDF <- file.path(repo_root, "final_plots", "pdf", "figure_stat4_vignette.pdf")

FIG_W_CM <- 18.0
FIG_H_CM <- 14.0

save_figure(fig, OUT, OUT_PDF, FIG_W_CM, FIG_H_CM)
