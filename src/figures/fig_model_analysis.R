#!/usr/bin/env Rscript
# fig_model_analysis.R
#
# Six-panel composite figure (ggplot2/patchwork) for the model-analysis section.

for (loc in c("en_US.UTF-8", "C.UTF-8", "UTF-8")) {
  if (suppressWarnings(Sys.setlocale("LC_ALL", loc)) != "") break
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})

# Resolve paths relative to the script location, not the working directory.
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

# Shared style
col_genetic <- "#E69F00"
col_perturbational  <- "#CC79A7"
col_obs     <- "#0072B2"
col_full    <- "#66C2A5"
col_weak    <- "#F0E442"
col_strong  <- "#D55E00"

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

LW       <- function(pt) pt / .pt
HAIRLINE <- 0.8
RULE     <- 1.2

base_theme <- theme_classic(base_size = AXIS_TITLE) +
  theme(text        = element_text(family = font_family),
        axis.title  = element_text(size = AXIS_TITLE),
        axis.text   = element_text(size = AXIS_TEXT),
        axis.line   = element_line(linewidth = LW(HAIRLINE)),
        axis.ticks  = element_line(linewidth = LW(HAIRLINE)),
        plot.title  = element_text(size = TITLE_SIZE, hjust = 0),
        plot.title.position = "plot",
        plot.tag    = element_text(size = TAG_SIZE, face = "bold"),
        plot.margin = margin(6, 10, 6, 10))

# Display names for panels A, B and C, per the panel-letter mapping in the file header.
LADDER_LABELS <- c("genetic only"    = "Genetics-\nOnly",
                   "+ observational" = "Genetics +\nObservational\nBlock",
                   "+ perturbational" = "Genetics +\nPerturbational\nBlock",
                   "full"            = "Full\nModel")
MARGINAL_LABELS <- c(perturbational = "Genetics + Perturbational Block\nvs Genetics-Only",
                     observational  = "Genetics + Observational Block\nvs Genetics-Only",
                     full           = "Full Model\nvs Genetics-Only")
SHAPLEY_LABELS  <- c(genetic        = "Genetics Block",
                     perturbational = "Perturbational Block",
                     observational  = "Observational Block")
DEPTH_LABELS    <- c("top-50" = "Top 50", "top-100" = "Top 100",
                     "top-200" = "Top 200", "top-500" = "Top 500")
BLOCK_LABELS <- c("Trans-regulatory residual"           = "Residualized trans-regulatory burden",
                  "Regulatory burden"                   = "Regulator burden correlation",
                  "Cytokine regulation"                 = "Cytokine regulation",
                  "Cytokine-receptor regulation"        = "Cytokine-receptor regulation",
                  "Polarization knockdown-coefficients" = "Ranked polarization coefficient",
                  "Coverage flags"                      = "Coverage flags")

RECOVERY_LEVELS <- c("Genetics-Only", "Full Model")
STRATUM_LEVELS  <- c("Weak Genetic Prior", "Strong Genetic Prior")

pfmt <- function(p) {
  ifelse(p < 1e-4, sprintf("p=%.0e", p), sprintf("p=%.4f", p))
}

# Panel A: AUC ladder
A <- read.csv(need("panelA_auc_ladder.csv"), stringsAsFactors = FALSE)
A$label <- LADDER_LABELS[A$model]
stopifnot(!any(is.na(A$label)))
A$x     <- seq_len(nrow(A))
A$col   <- c(col_genetic, col_obs, col_perturbational, col_full)

Asig <- read.csv(need("panelA_significance.csv"), stringsAsFactors = FALSE)
xpos <- c(genetic = 1, "genetic+observational" = 2, "genetic+perturbational" = 3, full = 4)
vs_gen <- Asig[Asig$model_a == "genetic" &
               Asig$model_b %in% c("genetic+observational", "genetic+perturbational", "full"), ]
vs_gen <- vs_gen[match(c("genetic+observational", "genetic+perturbational", "full"), vs_gen$model_b), ]
vs_gen$x1 <- xpos[["genetic"]]
vs_gen$x2 <- xpos[vs_gen$model_b]
vs_gen$y  <- c(0.797, 0.803, 0.809)
bt <- 0.0015

panelA <- ggplot(A, aes(x = x, y = auc)) +
  geom_line(colour = "#999999", linewidth = LW(RULE)) +
  geom_point(aes(colour = col), size = 2.6) +
  geom_text(aes(label = sprintf("%.3f", auc)), vjust = -1.1, size = GEOM_PT(DENSE)) +
  geom_segment(data = vs_gen, aes(x = x1, xend = x2, y = y, yend = y),
               inherit.aes = FALSE, linewidth = LW(HAIRLINE)) +
  geom_segment(data = vs_gen, aes(x = x1, xend = x1, y = y - bt, yend = y),
               inherit.aes = FALSE, linewidth = LW(HAIRLINE)) +
  geom_segment(data = vs_gen, aes(x = x2, xend = x2, y = y - bt, yend = y),
               inherit.aes = FALSE, linewidth = LW(HAIRLINE)) +
  geom_text(data = vs_gen, aes(x = (x1 + x2) / 2, y = y + 0.0007, label = stars),
            inherit.aes = FALSE, size = GEOM_PT(TITLE_SIZE)) +
  scale_colour_identity() +
  scale_x_continuous(breaks = A$x, labels = A$label, expand = expansion(add = 0.15)) +
  scale_y_continuous() +
  coord_cartesian(xlim = c(0.6, 4.5), ylim = c(0.755, 0.813)) +
  labs(x = NULL, y = "CV ROC-AUC\n(approved vs unlabeled)",
       title = "Adding functional genomics to genetic priors", tag = "a") +
  base_theme +
  theme(axis.text.x = element_text(size = DENSE))

# Panel B: Shapley shares (stacked)
B <- read.csv(need("panelB_shapley_shares.csv"), stringsAsFactors = FALSE)
meta <- jsonlite::fromJSON(need("panels_ABC_meta.json"))
B$evidence <- factor(B$evidence, levels = c("observational", "perturbational", "genetic"))
B <- B[order(B$evidence), ]
B$col <- c(observational = col_obs, perturbational = col_perturbational, genetic = col_genetic)[as.character(B$evidence)]
B <- B[order(match(B$evidence, c("genetic", "perturbational", "observational"))), ]
B$ymax <- cumsum(B$shapley_auc)
B$ymin <- B$ymax - B$shapley_auc
B$ymid <- (B$ymin + B$ymax) / 2
B$lab  <- sprintf("%s\n%.0f%%", SHAPLEY_LABELS[as.character(B$evidence)], B$share_pct)
B$txtcol <- ifelse(B$evidence == "genetic", "black", "white")

panelB <- ggplot(B) +
  geom_rect(aes(xmin = 0.6, xmax = 1.4, ymin = ymin, ymax = ymax, fill = col),
            colour = "white", linewidth = LW(RULE), alpha = 0.8) +
  geom_text(aes(x = 1, y = ymid, label = lab, colour = txtcol),
            size = GEOM_PT(DENSE), lineheight = 0.9) +
  annotate("text", x = 1, y = meta$full_above_chance + 0.006,
           label = sprintf("total %.3f", meta$full_above_chance),
           size = GEOM_PT(DENSE), colour = "#555") +
  scale_fill_identity() + scale_colour_identity() +
  scale_x_continuous() +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  coord_cartesian(xlim = c(0.5, 1.5), ylim = c(0, meta$full_above_chance + 0.012)) +
  labs(x = NULL, y = "Shapley value (AUC above chance)",
       title = "Share of model performance by evidence type", tag = "b") +
  base_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

# Panel C: marginal contribution (paired)
Cc <- read.csv(need("panelC_marginal.csv"), stringsAsFactors = FALSE)
Cc$comparison <- MARGINAL_LABELS[Cc$evidence]
stopifnot(!any(is.na(Cc$comparison)))
# Rows read top-to-bottom in panel A's ladder order (observational, then perturbational).
# Keyed on evidence, not CSV row order, so a reordered table cannot silently relabel the axis.
C_ORDER <- c(observational = 3, perturbational = 2, full = 1)
Cc$y   <- C_ORDER[Cc$evidence]
stopifnot(!any(is.na(Cc$y)))
Cc$col <- c(perturbational = col_perturbational, observational = col_obs, full = col_full)[Cc$evidence]

panelC <- ggplot(Cc, aes(x = delta, y = y)) +
  geom_col(aes(fill = col), width = 0.6, orientation = "y") +
  geom_errorbar(aes(xmin = delta - sem, xmax = delta + sem), orientation = "y",
                width = 0.18, colour = "#333333", linewidth = LW(RULE)) +
  geom_text(aes(x = delta + sem + 0.0012, label = pfmt(nb_p)),
            hjust = 0, size = GEOM_PT(DENSE), colour = "#555") +
  scale_fill_identity() +
  scale_y_continuous(breaks = Cc$y, labels = Cc$comparison) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(0, 0.042)) +
  labs(x = "\u0394AUC over Genetics-Only (SEM)", y = NULL,
       title = "Marginal contribution (Nadeau-Bengio corrected)", tag = "c") +
  base_theme +
  theme(axis.text.y = element_text(size = DENSE))

# Panel D: score distribution by target class (violin+box)
SG   <- read.csv(need("score_by_group.csv"), stringsAsFactors = FALSE)
SGst <- jsonlite::fromJSON(need("score_by_group_stats.json"))
sg_lev  <- c("non-target", "in-trial", "approved")
SG$grp  <- factor(SG$grp, levels = sg_lev)
sg_cols <- c("non-target" = "#999999", "in-trial" = "#009E73", "approved" = "#56B4E9")
sg_xi   <- setNames(seq_along(sg_lev), sg_lev)
sg_pw   <- SGst$pairwise_mannwhitney
sg_br <- data.frame(
  x1  = c(sg_xi[["non-target"]], sg_xi[["in-trial"]], sg_xi[["non-target"]]),
  x2  = c(sg_xi[["in-trial"]],   sg_xi[["approved"]], sg_xi[["approved"]]),
  y   = c(1.06, 1.06, 1.16),
  lab = c(sg_pw$stars[sg_pw$a == "non-target" & sg_pw$b == "in-trial"],
          sg_pw$stars[sg_pw$a == "in-trial"   & sg_pw$b == "approved"],
          sg_pw$stars[sg_pw$a == "non-target" & sg_pw$b == "approved"]))
sg_tick <- 0.015

panelD <- ggplot(SG, aes(grp, rank_pctile_cv, fill = grp)) +
  geom_violin(scale = "width", width = 0.85, colour = NA, alpha = 0.55, trim = TRUE) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", colour = "#333333",
               linewidth = LW(HAIRLINE)) +
  scale_fill_manual(values = sg_cols) +
  geom_segment(data = sg_br, aes(x = x1, xend = x2, y = y, yend = y),
               inherit.aes = FALSE, linewidth = LW(HAIRLINE)) +
  geom_segment(data = sg_br, aes(x = x1, xend = x1, y = y - sg_tick, yend = y),
               inherit.aes = FALSE, linewidth = LW(HAIRLINE)) +
  geom_segment(data = sg_br, aes(x = x2, xend = x2, y = y - sg_tick, yend = y),
               inherit.aes = FALSE, linewidth = LW(HAIRLINE)) +
  geom_text(data = sg_br, aes(x = (x1 + x2) / 2, y = y + 0.012, label = lab),
            inherit.aes = FALSE, colour = "#111111", size = GEOM_PT(TITLE_SIZE)) +
  annotate("text", x = 0.6, y = 1.30, hjust = 0, size = GEOM_PT(DENSE), colour = "#555",
           label = sprintf("In-trial vs Other: p = %.1e",
                           sg_pw$p_holm[sg_pw$a == "non-target" & sg_pw$b == "in-trial"])) +
  scale_x_discrete(labels = c(
    sprintf("Other\n(n = %s)",     format(SGst$groups$`non-target`$n, big.mark = ",")),
    sprintf("In-trial\n(n = %d)", SGst$groups$`in-trial`$n),
    sprintf("Approved\n(n = %d)\ntraining positives",  SGst$groups$approved$n))) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.01, 0))) +
  coord_cartesian(ylim = c(0, 1.36)) +
  labs(x = NULL, y = "Cross-Fitted Model Score\n(Rank Percentile)",
       title = "Model score rises with clinical validation", tag = "d") +
  base_theme +
  theme(legend.position = "none", axis.text.x = element_text(size = AXIS_TEXT))

# Panel E: target recovery vs depth
D <- read.csv(need("panelD_recovery.csv"), stringsAsFactors = FALSE)
Dl <- rbind(
  data.frame(depth = D$depth, x = seq_len(nrow(D)) - 0.19, n = D$genetics,
             group = RECOVERY_LEVELS[1], col = col_genetic),
  data.frame(depth = D$depth, x = seq_len(nrow(D)) + 0.19, n = D$genetics_perturbseq,
             group = RECOVERY_LEVELS[2], col = col_full)
)
Dl$group <- factor(Dl$group, levels = RECOVERY_LEVELS)
gain_lab <- data.frame(x = seq_len(nrow(D)) + 0.19, n = D$genetics_perturbseq,
                       lab = sprintf("+%d", D$gain))

panelE <- ggplot(Dl, aes(x = x, y = n, fill = group)) +
  geom_col(width = 0.36) +
  geom_text(aes(label = n), vjust = -0.5, size = GEOM_PT(DENSE), show.legend = FALSE) +
  geom_text(data = gain_lab, aes(x = x, y = n, label = lab),
            inherit.aes = FALSE, hjust = 0.5, vjust = -1.9, size = GEOM_PT(DENSE),
            colour = col_full) +
  scale_fill_manual(values = setNames(c(col_genetic, col_full), RECOVERY_LEVELS),
                    name = NULL) +
  scale_x_continuous(breaks = seq_len(nrow(D)), labels = DEPTH_LABELS[D$depth]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(ylim = c(0, max(D$genetics_perturbseq) * 1.26)) +
  labs(x = NULL, y = "Number of held-out immune\ntrial targets recovered",
       title = "Full model recovers more immune\ntargets at top ranks", tag = "e") +
  base_theme +
  theme(legend.position        = "inside",
        legend.position.inside = c(0.02, 0.98),
        legend.justification    = c(0, 1),
        legend.text            = element_text(size = DENSE, family = font_family),
        legend.key.size        = unit(9, "pt"),
        axis.text.x            = element_text(size = DENSE))

# Panel F: which perturbational block adds incremental signal
E <- read.csv(need("panelE_block_signal.csv"), stringsAsFactors = FALSE)
stopifnot(all(E$block %in% names(BLOCK_LABELS)))
E$block <- factor(BLOCK_LABELS[E$block], levels = rev(BLOCK_LABELS[E$block]))
El <- rbind(
  data.frame(block = E$block, lift = E$weak_lift,   stratum = STRATUM_LEVELS[1]),
  data.frame(block = E$block, lift = E$strong_lift, stratum = STRATUM_LEVELS[2])
)
El$stratum <- factor(El$stratum, levels = STRATUM_LEVELS)

panelF <- ggplot(El, aes(x = lift, y = block, fill = stratum)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_vline(xintercept = 0, colour = "black", linewidth = LW(RULE)) +
  scale_fill_manual(values = setNames(c(col_weak, col_strong), STRATUM_LEVELS),
                    name = NULL, guide = guide_legend(nrow = 2)) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.08))) +
  labs(x = "\u0394AUC vs Genetic Prior-Only", y = NULL,
       title = "Perturb-seq features that add signal", tag = "f") +
  base_theme +
  theme(axis.text.y     = element_text(size = DENSE),
        legend.position        = "inside",
        legend.position.inside = c(0.98, 0.58),
        legend.justification   = c(1, 0.5),
        legend.text     = element_text(size = DENSE, family = font_family),
        legend.key.size = unit(9, "pt"))

# Compose and save
free_l <- function(p) {
  if (utils::packageVersion("patchwork") >= "1.2.0") patchwork::free(p, side = "l") else p
}
panelC_free <- patchwork::free(panelC, type = "label", side = "b")

fig <- free_l(panelA | panelB) / free_l(panelC_free | panelD) / free_l(panelE | panelF) +
  plot_layout(heights = c(1, 1, 1))

OUT     <- file.path(repo_root, "final_plots", "figure_model_analysis_composite.png")
OUT_PDF <- file.path(repo_root, "final_plots", "pdf", "figure_model_analysis_composite.pdf")
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)

FIG_W_CM <- 18.0
FIG_H_CM <- 20.0

save_figure(fig, OUT, OUT_PDF, FIG_W_CM, FIG_H_CM)
