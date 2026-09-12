#!/usr/bin/env Rscript
# fig_immune_stages.R
#
# Three-panel immune-indication trial-validation figure (Fig 5).
#
# Panel a: for the full model's top 50 non-approved PU nominations, the highest
# clinical-trial stage of an immune-indication drug against each target, encoded
# by bar length and colour; a grey stub marks targets with no immune-indication drug.
# Panel b: percentage of genes with an immune-indication drug for the full model's
# top 50, the genetics-only model's top 50, and a bootstrapped control, at two
# trial-stage thresholds (any immune trial, Ph I+; advanced, Ph III), with 95%
# bootstrap intervals, fold-enrichment and empirical p-values.
# Panel c: the any-immune-trial arm of panel b across shortlist depths 25/50/100.
#
# Inputs from figure_data/, written by build_immune_stages_data.py.

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

# --- Shared style ----------------------------------------------------------
col_top  <- "#2a7ab0"
col_gen  <- "#E69F00"
col_ctrl <- "#b0b0b0"

STAGE_NUM <- c("Phase 3" = 3, "Phase 2" = 2, "Phase 1" = 1)
STAGE_COL <- c("Phase 3" = "#6a51a3", "Phase 2" = "#9e9ac8", "Phase 1" = "#dadaeb")

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

top <- read.csv(need("immune_stages_top50.csv"), stringsAsFactors = FALSE)
top <- top[order(top$rank_nonapproved), ]

stage   <- trimws(as.character(top$furthest_imm_stage))
stage[stage == "" | stage == "NA"] <- NA
max_num <- unname(STAGE_NUM[stage]); max_num[is.na(max_num)] <- -1
top$y        <- rev(seq_len(nrow(top)))
top$plot_val <- ifelse(max_num < 0, 0.10, max_num)
top$fill     <- ifelse(max_num < 0, "#f0f0f0", unname(STAGE_COL[stage]))
top$fill[is.na(top$fill)] <- "#f0f0f0"

panelA <- ggplot(top, aes(x = plot_val, y = y)) +
  geom_col(aes(fill = fill), orientation = "y", width = 0.66,
           colour = "white", linewidth = LW(HAIRLINE)) +
  scale_fill_identity() +
  scale_x_continuous(breaks = 0:3, labels = c("none", "Phase I", "Phase II", "Phase III"),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(breaks = top$y, labels = top$gene,
                     expand = expansion(add = 0.6)) +
  coord_cartesian(xlim = c(0, 3.3)) +
  labs(x = "Highest stage of an immune-indication drug against the target", y = NULL,
       title = sprintf(
         "Top %d full model nominations excluding\napproved immune-indication targets", nrow(top)),
       tag = "a") +
  base_theme +
  theme(axis.text.y = element_text(size = TINY, face = "italic"))

# --- Panels b, c: enrichment vs bootstrapped control ------------------------
counts <- jsonlite::fromJSON(need("immune_stages_counts.json"))
tn     <- counts$N_top
sweep  <- counts$depth_sweep

grp_top  <- "Full model"
grp_gen  <- "Genetics-only model"
grp_ctrl <- sprintf("Bootstrapped controls\n(shortlist-size draws \u00d7 %s)",
                    format(counts$n_bootstrap, big.mark = ","))
GRPS   <- c(grp_top, grp_gen, grp_ctrl)
FILLS  <- setNames(c(col_top, col_gen, col_ctrl), GRPS)
C_ANN_TOP <- "#08519c"
C_ANN_GEN <- "#9a6a00"

pstr <- function(p) if (p < 1e-4) "p<0.0001" else sprintf("p=%.4f", p)
fmt_ann <- function(fold, p) sprintf("%.1f\u00d7\n%s", fold, sapply(p, pstr))

Y_MAX  <- 26
ANN_DY <- 1.6
DODGE  <- 0.26
BAR_W  <- 0.24

enrich_panel <- function(d, xlabs, title, tag, show_legend, gen_pct_hjust = 0.5) {
  stopifnot(nrow(d) == length(xlabs),
            abs(d$gen_ctrl - d$ctrl) < 0.1,
            c(d$top_pct, d$gen_pct, d$ctrl_hi) < Y_MAX)
  x <- seq_len(nrow(d))
  bars <- rbind(
    data.frame(x = x - DODGE, pct = d$top_pct, lo = NA,        hi = NA,        group = grp_top),
    data.frame(x = x,         pct = d$gen_pct, lo = NA,        hi = NA,        group = grp_gen),
    data.frame(x = x + DODGE, pct = d$ctrl,    lo = d$ctrl_lo, hi = d$ctrl_hi, group = grp_ctrl)
  )
  bars$group <- factor(bars$group, levels = GRPS)
  ann_mod <- rbind(
    data.frame(x = x - DODGE, pct = d$top_pct, colour = C_ANN_TOP, hj = 0.5, pct_hj = 0.5,
               lab = fmt_ann(d$top_fold, d$top_p)),
    data.frame(x = x,         pct = d$gen_pct, colour = C_ANN_GEN, hj = 0, pct_hj = gen_pct_hjust,
               lab = fmt_ann(d$gen_fold, d$gen_p))
  )

  p <- ggplot(bars, aes(x = x, y = pct, fill = group)) +
    geom_col(width = BAR_W) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.10, colour = "#444444",
                  linewidth = LW(RULE), na.rm = TRUE) +
    geom_text(data = ann_mod, aes(x = x, y = pct, label = sprintf("%.0f%%", pct), hjust = pct_hj),
              inherit.aes = FALSE, vjust = -0.4, size = GEOM_PT(TITLE_SIZE)) +
    geom_text(data = ann_mod,
              aes(x = x, y = pct + ANN_DY, label = lab, colour = colour, hjust = hj),
              inherit.aes = FALSE, size = GEOM_PT(DENSE), vjust = 0, show.legend = FALSE) +
    scale_colour_identity() +
    scale_fill_manual(values = FILLS, name = NULL) +
    scale_x_continuous(breaks = x, labels = xlabs) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(xlim = c(0.5, nrow(d) + 0.5), ylim = c(0, Y_MAX)) +
    labs(x = NULL, y = "% of genes with an in-trial immune-indication drug",
         title = title, tag = tag) +
    base_theme +
    theme(axis.text.x = element_text(size = AXIS_TEXT))

  if (show_legend) {
    p + theme(legend.position        = "inside",
              legend.position.inside = c(0.98, 0.98),
              legend.justification   = c(1, 1),
              legend.text            = element_text(size = TINY, family = font_family),
              legend.key.size        = unit(9, "pt"))
  } else {
    p + theme(legend.position = "none")
  }
}

dB <- data.frame(
  top_pct  = c(counts$top_phI_pct,           counts$top_ph3_pct),
  gen_pct  = c(counts$gen_top_phI_pct,       counts$gen_top_ph3_pct),
  ctrl     = c(counts$ctrl_phI_mean_pct,     counts$ctrl_ph3_mean_pct),
  ctrl_lo  = c(counts$ctrl_phI_ci[1],        counts$ctrl_ph3_ci[1]),
  ctrl_hi  = c(counts$ctrl_phI_ci[2],        counts$ctrl_ph3_ci[2]),
  gen_ctrl = c(counts$gen_ctrl_phI_mean_pct, counts$gen_ctrl_ph3_mean_pct),
  top_fold = c(counts$fold_phI,              counts$fold_ph3),
  gen_fold = c(counts$gen_fold_phI,          counts$gen_fold_ph3),
  top_p    = c(counts$bootstrap_p_phI_bh,     counts$bootstrap_p_ph3_bh),
  gen_p    = c(counts$gen_bootstrap_p_phI_bh, counts$gen_bootstrap_p_ph3_bh)
)
panelB <- enrich_panel(
  dB, c("Any immune trial\n(Phase I+)", "Advanced immune\ntrial (Phase III)"),
  sprintf("Top %d nominations enriched for\nimmune clinical trials", tn), "b",
  show_legend = TRUE)

sw_f <- sweep[sweep$model == "full", ]
sw_g <- sweep[sweep$model == "genetics_only", ]
sw_f <- sw_f[order(sw_f$depth), ]
sw_g <- sw_g[order(sw_g$depth), ]
stopifnot(identical(sw_f$depth, sw_g$depth), setequal(sw_f$depth, counts$panelC_depths))
dC <- data.frame(
  top_pct  = sw_f$top_phI_pct,        gen_pct  = sw_g$top_phI_pct,
  ctrl     = sw_f$ctrl_phI_mean_pct,  ctrl_lo  = sw_f$ctrl_phI_ci_lo,
  ctrl_hi  = sw_f$ctrl_phI_ci_hi,     gen_ctrl = sw_g$ctrl_phI_mean_pct,
  top_fold = sw_f$fold_phI,           gen_fold = sw_g$fold_phI,
  top_p    = sw_f$bootstrap_p_phI_bh, gen_p    = sw_g$bootstrap_p_phI_bh
)
i50 <- which(sw_f$depth == tn)
stopifnot(length(i50) == 1,
          isTRUE(all.equal(dC$top_pct[i50], counts$top_phI_pct)),
          isTRUE(all.equal(dC$gen_pct[i50], counts$gen_top_phI_pct)),
          isTRUE(all.equal(dC$top_fold[i50], counts$fold_phI)))
panelC <- enrich_panel(
  dC, sprintf("Top %d", sw_f$depth),
  "Any immune trial (Phase I+)\nacross shortlist depth", "c",
  show_legend = FALSE, gen_pct_hjust = 0.41)

# --- Compose and save --------------------------------------------------------
fig <- panelA + (panelB / panelC) + plot_layout(widths = c(1.4, 1))

OUT     <- file.path(repo_root, "final_plots", "figure_trial_validation_immune_stages_R.png")
OUT_PDF <- file.path(repo_root, "final_plots", "pdf", "figure_trial_validation_immune_stages_R.pdf")
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)

FIG_W_CM <- 18.0
FIG_H_CM <- 19.0

save_figure(fig, OUT, OUT_PDF, FIG_W_CM, FIG_H_CM)
