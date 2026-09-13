#!/usr/bin/env Rscript
# See REPRODUCIBILITY.md.

for (loc in c("en_US.UTF-8", "C.UTF-8", "UTF-8")) {
  if (suppressWarnings(Sys.setlocale("LC_ALL", loc)) != "") break
}
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(jsonlite) })

args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir  <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()
repo_root   <- normalizePath(file.path(script_dir, "..", ".."))
stopifnot(file.exists(file.path(repo_root, "Makefile")))
source(file.path(script_dir, "vector_output.R"))
need <- function(f) {
  p <- file.path(repo_root, "figure_data", f)
  if (!file.exists(p))
    stop("Missing ", p, "\nThis file is committed; restore it from git, or rebuild with ",
         "`make permutation-null` (~2 h). See figure_data/PROVENANCE.md.")
  p
}

null <- read.csv(need("label_permutation_null.csv"), stringsAsFactors = FALSE, check.names = FALSE)
pv   <- read.csv(need("label_permutation_pvalues.csv"), stringsAsFactors = FALSE, check.names = FALSE)

# Staleness guard for panel a.
ladder <- read.csv(file.path(repo_root, "figure_data", "panelA_auc_ladder.csv"),
                   stringsAsFactors = FALSE)
ladder_key <- c("genetic only" = "genetic", "+ observational" = "genetic+observational",
                "+ perturbational" = "genetic+perturbational", "full" = "full")
lad <- setNames(ladder$auc, ladder_key[ladder$model])
lev <- pv[pv$kind == "level", ]
drift <- abs(lad[lev$metric] - lev$observed)
if (any(is.na(drift)) || max(drift) > 1e-12) {
  bad <- lev$metric[is.na(drift) | drift > 1e-12]
  check_or_warn(paste0("label_permutation_pvalues.csv is stale against panelA_auc_ladder.csv.\n",
       "  disagreeing rungs: ", paste(bad, collapse = ", "), "\n",
       "  max |difference| : ", format(max(drift, na.rm = TRUE)), "\n",
       "Run `make permutation-null-finalize` to recompute the p-values against the ",
       "current ladder (seconds -- it reuses the committed null table)."))
}

scr_json <- jsonlite::fromJSON(need("scrambled_feature_control.json"))
scr_vals <- as.numeric(unlist(scr_json$scrambled))
scr_real <- as.numeric(scr_json$real)
scr_sum  <- scr_json$summary

# Staleness guard for panel b.
inc_obs <- pv$observed[pv$metric == "full"] - pv$observed[pv$metric == "genetic"]
if (length(inc_obs) != 1 || abs(inc_obs - scr_real) > 1e-12) {
  stop("scrambled_feature_control.json is stale against label_permutation_pvalues.csv.\n",
       "  real (scrambled control)      : ", format(scr_real, digits = 17), "\n",
       "  full - genetic (pvalues)      : ", format(inc_obs, digits = 17), "\n",
       "  |difference|                  : ", format(abs(inc_obs - scr_real)), "\n",
       "These are the same quantity. Rebuild with `make permutation-null` (~1.6 h for the ",
       "scrambled control), deleting outputs/permutation/ first so a stale checkpoint is not ",
       "resumed. See figure_data/PROVENANCE.md.")
}

col_map <- c(genetic = "#E69F00", `genetic+observational` = "#0072B2",
             `genetic+perturbational` = "#CC79A7", full = "#66C2A5")
lab_map <- c(genetic = "Genetics-Only", `genetic+observational` = "Genetics + Observational Block",
             `genetic+perturbational` = "Genetics + Perturbational Block", full = "Full Model")

font_family <- "Helvetica"
update_geom_defaults("text", list(family = font_family))

AXIS_TITLE <- 8
AXIS_TEXT  <- 7
TITLE_SIZE <- 8
TAG_SIZE   <- 8
STRIP_SIZE <- AXIS_TEXT
DENSE      <- 6.5
HAIRLINE   <- 0.8
RULE       <- 1.2

GEOM_PT <- function(pt) pt / .pt
LW      <- function(pt) pt / .pt

base_theme <- theme_classic(base_size = AXIS_TITLE) +
  theme(text = element_text(family = font_family),
        axis.title = element_text(size = AXIS_TITLE), axis.text = element_text(size = AXIS_TEXT),
        axis.line   = element_line(linewidth = LW(HAIRLINE)),
        axis.ticks  = element_line(linewidth = LW(HAIRLINE)),
        plot.title = element_text(size = TITLE_SIZE, hjust = 0),
        plot.tag = element_text(size = TAG_SIZE, face = "bold"),
        legend.position = "none", plot.margin = margin(6, 10, 6, 10),
        strip.text = element_text(size = STRIP_SIZE),
        strip.background = element_blank())

pfmt <- function(p) ifelse(p < 1/ (pv$n_perm[1] + 1) + 1e-12,
                           sprintf("p < %.1g", 1 / (pv$n_perm[1] + 1)),
                           sprintf("p = %.3f", p))

rungs <- c("genetic", "genetic+observational", "genetic+perturbational", "full")
la <- do.call(rbind, lapply(rungs, function(r)
  data.frame(rung = r, auc = null[[r]], stringsAsFactors = FALSE)))
la$rung <- factor(la$rung, levels = rungs, labels = lab_map[rungs])
obsA <- pv[pv$kind == "level", ]
obsA$rung <- factor(obsA$metric, levels = rungs, labels = lab_map[rungs])
obsA$col  <- col_map[obsA$metric]
obsA$lab  <- sprintf("Obs = %.3f\n%s", obsA$observed, pfmt(obsA$emp_p))

panelA <- ggplot(la, aes(x = auc)) +
  geom_histogram(bins = 40, fill = "#bbbbbb", colour = NA) +
  geom_vline(data = obsA, aes(xintercept = observed, colour = I(col)),
             linewidth = LW(RULE)) +
  geom_text(data = obsA, aes(x = observed, y = Inf, label = lab, colour = I(col)),
            hjust = 1.06, vjust = 1.3, size = GEOM_PT(DENSE), family = font_family, lineheight = 0.9) +
  facet_wrap(~ rung, ncol = 1, scales = "free_y") +
  coord_cartesian(xlim = c(0.40, 0.82)) +
  labs(title = sprintf("Label-permutation null (%d permutations)", pv$n_perm[1]),
       x = "5-fold CV ROC-AUC", y = "Permutations") +
  base_theme

# Scrambled-feature control: correct test for the FG increment, not the label-permutation null.
fg_green <- "#66C2A5"
sb <- data.frame(val = sort(scr_vals))
sb$rank <- seq_len(nrow(sb))
n_scr   <- nrow(sb)
n_ge_real <- as.integer(scr_sum$n_ge_real)
stopifnot(length(n_ge_real) == 1, !is.na(n_ge_real),
          n_ge_real == sum(scr_vals >= scr_real))
p_lab   <- sprintf("empirical p = %.3f  (%d of %d scrambles \u2265 real)",
                   scr_sum$emp_p, n_ge_real, n_scr)
real_lab <- sprintf("Obs Full Model Gain = +%.4f", scr_real)
scr_mean_lab <- sprintf("scrambled: %+.4f \u00b1 %.4f", scr_sum$scrambled_mean, scr_sum$scrambled_sd)

panelB <- ggplot(sb, aes(x = val, y = rank)) +
  geom_vline(xintercept = 0, colour = "#999999", linetype = "dotted", linewidth = LW(HAIRLINE)) +
  geom_vline(xintercept = scr_sum$scrambled_mean, colour = "#888888", linewidth = LW(HAIRLINE)) +
  geom_point(colour = "#666666", fill = "#cccccc", shape = 21, size = 1.7, stroke = LW(HAIRLINE)) +
  geom_vline(xintercept = scr_real, colour = fg_green, linewidth = LW(RULE)) +
  # vjust clears the top few draws: the label runs left from the rule, and the highest-ranked
  # points sit just under it. The draws fall away to the left as rank drops, so nudging the
  # label down moves it out of their path rather than into it.
  annotate("text", x = scr_real, y = n_scr, label = real_lab, colour = fg_green,
           hjust = 1.06, vjust = 3.0, size = GEOM_PT(DENSE), family = font_family) +
  annotate("text", x = scr_sum$scrambled_mean, y = 1, label = scr_mean_lab, colour = "#555555",
           hjust = -0.05, vjust = -0.4, size = GEOM_PT(DENSE), family = font_family) +
  annotate("text", x = Inf, y = -Inf, label = p_lab, colour = "#222222",
           hjust = 1.03, vjust = -0.8, size = GEOM_PT(DENSE), family = font_family) +
  scale_y_continuous(breaks = NULL, expand = expansion(mult = c(0.10, 0.10))) +
  labs(title = sprintf("Scrambled-feature control (%d draws)", length(scr_vals)),
       x = "\u0394 CV ROC-AUC (Genetics-Only \u2192 Full Model)",
       y = "Scrambled-feature draws (ranked)") +
  base_theme +
  theme(axis.title.x = element_text(family = "Arial Unicode MS"))

fig <- (panelA | panelB) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = TAG_SIZE, family = font_family, face = "bold"))


OUT     <- file.path(repo_root, "final_plots", "figure_label_permutation.png")
OUT_PDF <- file.path(repo_root, "final_plots", "pdf", "figure_label_permutation.pdf")

FIG_W_CM <- 18.0
FIG_H_CM <- FIG_W_CM * 7.4 / 9.6

save_figure(fig, OUT, OUT_PDF, FIG_W_CM, FIG_H_CM)
