#!/usr/bin/env Rscript
# See REPRODUCIBILITY.md.

for (loc in c("en_US.UTF-8", "C.UTF-8", "UTF-8")) {
  if (suppressWarnings(Sys.setlocale("LC_ALL", loc)) != "") break
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# Repo-root anchor: paths are computed relative to this script's location, src/figures/.
args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir  <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()
repo_root   <- normalizePath(file.path(script_dir, "..", ".."))
stopifnot(file.exists(file.path(repo_root, "Makefile")))

source(file.path(script_dir, "palette.R"))
source(file.path(script_dir, "vector_output.R"))

FIGDATA  <- file.path(repo_root, "figure_data")
OUT_DIR  <- file.path(repo_root, "final_plots")
OUT      <- file.path(OUT_DIR, "genetic_data.png")
OUT_PDF  <- file.path(OUT_DIR, "pdf", "genetic_data.pdf")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Shared style
font_family <- "Helvetica"
update_geom_defaults("text", list(family = font_family))

AXIS_TITLE <- 8
AXIS_TEXT  <- 7
TAG_SIZE   <- 8
GEOM_PT    <- AXIS_TEXT / .pt
PCT_PT     <- 6.5 / .pt

LW       <- function(pt) pt / .pt
STROKE   <- function(pt) pt / (.stroke / 2)
HAIRLINE <- 0.8
RULE     <- 1.2

base_theme <- theme_classic(base_size = AXIS_TITLE) +
  theme(text        = element_text(family = font_family),
        axis.title  = element_text(size = AXIS_TITLE),
        axis.text   = element_text(size = AXIS_TEXT),
        axis.line   = element_line(linewidth = LW(HAIRLINE)),
        axis.ticks  = element_line(linewidth = LW(HAIRLINE)),
        plot.title  = element_text(size = AXIS_TITLE, hjust = 0),
        plot.tag    = element_text(size = TAG_SIZE, face = "bold"),
        plot.margin = margin(6, 10, 6, 10))

# Panel a: IEI-gene enrichment forest
csv_path <- file.path(FIGDATA, "iei_enrichment_forest.csv")
stopifnot(file.exists(csv_path))
res <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)

axis_n <- read.csv(file.path(FIGDATA, "iei_enrichment_axis.csv"), stringsAsFactors = FALSE)
GROUP_KEY <- c("Targets of approved immune drugs"  = "approved",
               "Targets of immune drugs in trials" = "in_trials",
               "All other genes"                   = "other",
               "Druggable genome"                  = "other")
BASE_Y <- c(approved = 3, in_trials = 2, other = 1)
DODGE  <- 0.15
lab_prim <- "Primary: vs all other genes"
lab_sens <- "Sensitivity: vs druggable genome (Finan et al., 2017)"

pa <- res
pa$ckey <- GROUP_KEY[pa$group]
pa$hue  <- TARGET_COLORS[pa$ckey]
pa$bg   <- factor(ifelse(pa$background == "all_genes", lab_prim, lab_sens),
                  levels = c(lab_prim, lab_sens))
pa$dy   <- ifelse(pa$background == "all_genes", DODGE, -DODGE)
pa$y    <- BASE_Y[pa$ckey] + pa$dy
pa$fill_col <- ifelse(pa$bg == lab_prim, pa$hue, "white")
pa$grp  <- factor(pa$ckey, levels = c("approved", "in_trials", "other"))
grp_values <- c(approved = unname(TARGET_COLORS["approved"]),
                in_trials = unname(TARGET_COLORS["in_trials"]),
                other = unname(TARGET_COLORS["other"]))
fmt_n <- function(x) formatC(x, format = "d", big.mark = ",")
pa$txt <- pa$label_text
pa$txt_x <- ifelse(is.na(pa$OR_hi), pa$OR, pa$OR_hi) + 0.3
ylabels <- c(
  sprintf("Approved\n(n = %s)",  fmt_n(axis_n$n[axis_n$ckey == "approved"])),
  sprintf("In-trial\n(n = %s)", fmt_n(axis_n$n[axis_n$ckey == "in_trials"])),
  "Reference\n(OR = 1)")
xlim_hi <- max(pa$OR_hi, na.rm = TRUE) + 6.5

panelA <- ggplot(pa, aes(x = OR, y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#c0392b", linewidth = LW(RULE)) +
  annotate("text", x = 1.1, y = 3.7, label = "No enrichment (OR = 1)",
           colour = "#c0392b", size = GEOM_PT, hjust = 0, family = font_family) +
  geom_errorbarh(data = pa[!is.na(pa$OR_lo), ],
                 aes(xmin = OR_lo, xmax = OR_hi, colour = grp), height = 0, linewidth = LW(RULE)) +
  geom_point(aes(fill = bg), shape = 21, alpha = 0) +
  geom_point(aes(colour = grp), shape = 21, size = 2.9, stroke = STROKE(RULE),
             fill = pa$fill_col) +
  geom_text(aes(x = txt_x, label = txt), hjust = 0, size = GEOM_PT) +
  scale_colour_manual(values = grp_values, name = NULL, guide = "none") +
  scale_fill_manual(values = c("grey30", "white"), breaks = c(lab_prim, lab_sens), name = NULL,
                    guide = guide_legend(order = 2,
                                         override.aes = list(shape = 21, colour = "grey30", alpha = 1))) +
  scale_x_continuous(breaks = seq(0, floor(max(pa$OR_hi, na.rm = TRUE)), by = 2.5),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(breaks = c(3, 2, 1), labels = ylabels) +
  coord_cartesian(xlim = c(0, xlim_hi), ylim = c(0.4, 4.0), clip = "off") +
  labs(x = "Odds ratio of IEI-gene enrichment (95% CI)", y = NULL, tag = "a") +
  base_theme +
  theme(axis.line.y            = element_blank(),
        axis.ticks.y           = element_blank(),
        legend.position        = "inside",
        legend.position.inside = c(0.46, 0.30),
        legend.justification   = c(0, 1),
        legend.background      = element_blank(),
        legend.key             = element_blank(),
        legend.text            = element_text(size = AXIS_TEXT),
        legend.key.size        = unit(10, "pt"))

# Panel b: immune-GWAS score violin
plot_data_pos <- read.csv(file.path(FIGDATA, "gwas_by_drug_status.csv"),
                          stringsAsFactors = FALSE)
grp_n_tab     <- read.csv(file.path(FIGDATA, "genetic_data_group_n.csv"),
                          stringsAsFactors = FALSE)
zero_inf      <- read.csv(file.path(FIGDATA, "gwas_zero_inflation.csv"),
                          stringsAsFactors = FALSE)

n_gwas <- function(g) as.integer(grp_n_tab$n[grp_n_tab$panel == "gwas" &
                                             grp_n_tab$group == g])
lab_approved  <- sprintf("Approved\n(n = %s)",  format(n_gwas("approved"),  big.mark = ","))
lab_in_trials <- sprintf("In-trial\n(n = %s)", format(n_gwas("in_trials"), big.mark = ","))
lab_other     <- sprintf("Other\n(n = %s)",     format(n_gwas("other"),     big.mark = ","))
group_labels  <- c(approved = lab_approved, in_trials = lab_in_trials, other = lab_other)
lev_b         <- c(lab_approved, lab_in_trials, lab_other)

plot_data_pos$Target <- factor(group_labels[plot_data_pos$group], levels = lev_b)
fill_palette <- setNames(unname(TARGET_COLORS[c("approved", "in_trials", "other")]), lev_b)

pct_pos <- data.frame(
  Target = factor(group_labels[zero_inf$group], levels = lev_b),
  lab    = sprintf("%.0f%% > 0", 100 * zero_inf$pct_gt0),
  stringsAsFactors = FALSE)

panelB <- ggplot(plot_data_pos, aes(x = Target, y = gwas_score, fill = Target)) +
  geom_violin(trim = TRUE, colour = NA, alpha = 0.75, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
               alpha = 0.6, colour = "grey30", linewidth = LW(HAIRLINE)) +
  geom_text(data = pct_pos, aes(x = Target, y = 1.02, label = lab),
            inherit.aes = FALSE, size = PCT_PT, colour = "grey30", vjust = 0) +
  coord_cartesian(ylim = c(0, 1.08)) +
  scale_fill_manual(values = fill_palette) +
  labs(x = NULL, y = "Immune-GWAS score, genes with signal (>0)", tag = "b") +
  base_theme +
  theme(legend.position = "none",
        axis.title.y    = element_text(size = AXIS_TITLE, margin = margin(r = 2)))

# Panel c: missense-constraint violin (mis.z_score).
mis_data <- read.csv(file.path(FIGDATA, "mis_z_by_drug_status.csv"),
                     stringsAsFactors = FALSE)
n_mis <- function(g) as.integer(grp_n_tab$n[grp_n_tab$panel == "mis_z" &
                                            grp_n_tab$group == g])
lab_approved_c  <- sprintf("Approved\n(n = %s)",  format(n_mis("approved"),  big.mark = ","))
lab_in_trials_c <- sprintf("In-trial\n(n = %s)", format(n_mis("in_trials"), big.mark = ","))
lab_other_c     <- sprintf("Other\n(n = %s)",     format(n_mis("other"),     big.mark = ","))
group_labels_c  <- c(approved = lab_approved_c, in_trials = lab_in_trials_c, other = lab_other_c)
lev_c           <- c(lab_approved_c, lab_in_trials_c, lab_other_c)

mis_data$Target <- factor(group_labels_c[mis_data$group], levels = lev_c)
fill_palette_c <- setNames(unname(TARGET_COLORS[c("approved", "in_trials", "other")]), lev_c)

panelC <- ggplot(mis_data, aes(x = Target, y = mis_z_score, fill = Target)) +
  geom_violin(trim = TRUE, colour = NA, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50", linewidth = LW(HAIRLINE)) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", colour = "grey30",
               linewidth = LW(HAIRLINE)) +
  coord_cartesian(ylim = c(-5, 8)) +
  scale_fill_manual(values = fill_palette_c) +
  labs(x = NULL, y = "Missense Z-score", tag = "c") +
  base_theme +
  theme(legend.position = "none",
        axis.title.y    = element_text(size = AXIS_TITLE, margin = margin(r = 2)))

# Compose and save
bottom <- panelB | panelC
bottom_free <- if (utils::packageVersion("patchwork") >= "1.2.0") {
  patchwork::free(bottom, side = "l")
} else bottom
fig <- panelA / bottom_free + plot_layout(heights = c(1.15, 1))

FIG_W_CM <- 18.0
FIG_H_CM <- 16.3

save_figure(fig, OUT, OUT_PDF, FIG_W_CM, FIG_H_CM)
