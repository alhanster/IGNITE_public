#!/usr/bin/env Rscript
# fig_discordance.R
#
# Produces the genetics-vs-full-model discordance figure corresponding to Supplementary S2.

for (loc in c("en_US.UTF-8", "C.UTF-8", "UTF-8")) {
  if (suppressWarnings(Sys.setlocale("LC_ALL", loc)) != "") break
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(jsonlite)
  library(ggrepel)
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
need_cache <- function(f) {
  p <- file.path(repo_root, "figure_data", f)
  if (!file.exists(p))
    stop("Missing ", p,
         "\nThis is a COMMITTED table from the opt-in STRING/GO cache, not a `make tables`",
         "\noutput. Restore it from git (`git checkout -- figure_data`), or rebuild it with",
         "\n  make discordance-network   (~10 min, one-time 135 MB download)",
         "\nthen re-run `make tables` to fold it back into discordance_stats.json.")
  p
}
OUT <- file.path(repo_root, "final_plots", "supplementary",
                 "figure_discordance_genetics_vs_full.png")

COL_FULL <- "#66C2A5"
COL_GEN  <- "#E69F00"
GREY     <- "#666666"
LIGHT    <- "#CFCFCF"

font_family <- "Helvetica"
update_geom_defaults("text", list(family = font_family))

AXIS_TITLE <- 6
AXIS_TEXT  <- 6
TITLE_SIZE <- 7
TAG_SIZE   <- 8
DENSE      <- 6
TINY       <- 6
HAIRLINE   <- 0.8
RULE       <- 1.2

GEOM_PT <- function(pt) pt / .pt
LW      <- function(pt) pt / .pt

theme_set(theme_classic(base_size = AXIS_TITLE) +
  theme(text         = element_text(family = font_family),
        plot.title   = element_text(size = TITLE_SIZE, hjust = 0.5),
        axis.title   = element_text(size = AXIS_TITLE),
        axis.text    = element_text(size = AXIS_TEXT, colour = "black"),
        legend.text  = element_text(size = TINY),
        legend.title = element_blank(),
        legend.key.size = unit(7, "pt"),
        legend.background = element_blank(),
        plot.tag     = element_text(size = TAG_SIZE, face = "bold")))

d <- read.csv(need("discordance_scatter.csv"), check.names = FALSE)
S <- fromJSON(need("discordance_stats.json"))

# Key gate runs before the first key value is read; downstream logic assumes a valid key.
REQUIRED_KEYS <- c("hi_cut", "lo_cut",
                   "n_features_full", "n_features_fg", "spearman_arms",
                   "n_core_promoted", "n_core_demoted",
                   "net_promoted_edges", "net_promoted_null_mean", "net_promoted_p",
                   "net_promoted_connected", "net_demoted_edges", "net_demoted_p",
                   "n_null_network")
missing_keys <- REQUIRED_KEYS[!REQUIRED_KEYS %in% names(S) |
                              vapply(REQUIRED_KEYS,
                                     function(k) length(S[[k]]) == 0L, logical(1))]
if (length(missing_keys)) {
  stop("discordance_stats.json is missing key(s): ",
       paste(missing_keys, collapse = ", "),
       "\n  Refresh it:  make tables  ->  make discordance-network  ->  make tables")
}

prom <- d %>% filter(core == "promoted") %>% arrange(desc(rank_shift))
dem  <- d %>% filter(core == "demoted")  %>% arrange(rank_shift)

# Panel a: the two gene rankings, with core genes labeled by name.
rest <- d %>% filter(core == "none")
GRP_LAB <- c(sprintf("promoted (n=%d)", S$n_core_promoted),
             sprintf("demoted (n=%d)",  S$n_core_demoted))
core_pts <- bind_rows(prom, dem) %>%
  mutate(grp = factor(core, levels = c("promoted", "demoted"), labels = GRP_LAB))
LABEL_GENES <- c("TAF11", "LSM14A", "EDC3", "IER5", "TMEM273", "SLC26A4")
lab_a <- bind_rows(prom, dem) %>%
  filter(gene %in% LABEL_GENES) %>%
  mutate(grp = factor(core, levels = c("promoted", "demoted"), labels = GRP_LAB))
lab_missing <- setdiff(LABEL_GENES, lab_a$gene)
if (length(lab_missing)) {
  stop("panel a's label list names gene(s) that are not in the core set: ",
       paste(lab_missing, collapse = ", "),
       "\n  These six are hand-picked. If the discordance set moved, re-pick them against",
       "\n  figure_data/discordance_core_genes.csv rather than deleting the name here.")
}

PROM_LAB_FLOOR <- max(prom$pctile_full) * 100 + 2.5

pa <- ggplot() +
  geom_point(data = rest, aes(pctile_genetics * 100, pctile_full * 100),
             colour = LIGHT, size = 0.22, stroke = 0) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              colour = GREY, linewidth = LW(HAIRLINE)) +
  geom_vline(xintercept = c(S$lo_cut, S$hi_cut) * 100, linetype = "dashed",
             colour = GREY, linewidth = LW(HAIRLINE)) +
  geom_hline(yintercept = c(S$lo_cut, S$hi_cut) * 100, linetype = "dashed",
             colour = GREY, linewidth = LW(HAIRLINE)) +
  geom_point(data = core_pts, aes(pctile_genetics * 100, pctile_full * 100,
                                  colour = grp), size = 1.25, stroke = 0) +
  geom_text_repel(data = filter(lab_a, core == "promoted"),
                  aes(pctile_genetics * 100, pctile_full * 100,
                      label = gene, colour = grp),
                  family = font_family,
                  size = GEOM_PT(DENSE), fontface = "italic", segment.size = 0.2,
                  segment.alpha = 0.85, min.segment.length = 0,
                  box.padding = 0.42, point.padding = 0.25,
                  nudge_x = -14, ylim = c(PROM_LAB_FLOOR, NA),
                  direction = "y", force = 6, force_pull = 0.3,
                  max.overlaps = Inf, seed = 4, show.legend = FALSE) +
  geom_text_repel(data = filter(lab_a, core == "demoted"),
                  aes(pctile_genetics * 100, pctile_full * 100,
                      label = gene, colour = grp),
                  family = font_family,
                  size = GEOM_PT(DENSE), fontface = "italic", segment.size = 0.2,
                  segment.alpha = 0.85, min.segment.length = 0,
                  box.padding = 0.42, point.padding = 0.25,
                  nudge_x = 14, ylim = c(NA, S$lo_cut * 100 - 1),
                  direction = "y", force = 6, force_pull = 0.3,
                  max.overlaps = Inf, seed = 4, show.legend = FALSE) +
  scale_colour_manual(values = setNames(c(COL_FULL, COL_GEN), GRP_LAB)) +
  coord_fixed(ratio = 1, xlim = c(0, 100), ylim = c(0, 100), clip = "off") +
  labs(x = "Genetics-Only PU Percentile", y = "Full Model PU Percentile",
       title = "Discordant genes between Full Model and Genetics-Only") +
  theme(plot.title = element_text(hjust = 0),
        legend.position = "inside",
        legend.position.inside = c(0.98, 0.98),
        legend.justification.inside = c(1, 1),
        legend.background = element_rect(fill = scales::alpha("white", 0.8), colour = NA),
        legend.margin = margin(2, 4, 2, 4))

# Panel b: tests whether the genes marked as discordant are biologically similar to one another.
net_all  <- read.csv(need_cache("discordance_network.csv"))
edg_all  <- read.csv(need_cache("discordance_network_edges.csv"))
nod_all  <- read.csv(need_cache("discordance_network_nodes.csv"))
expect_n <- c(promoted = nrow(prom), demoted = nrow(dem))
got_n    <- setNames(net_all$n_genes, net_all$direction)[names(expect_n)]

problems <- character(0)
if (!identical(as.integer(c(S$n_core_promoted, S$n_core_demoted)),
               as.integer(c(nrow(prom), nrow(dem))))) {
  problems <- c(problems, sprintf(
    "discordance_stats.json says %d/%d core genes (promoted/demoted), discordance_scatter.csv has %d/%d",
    S$n_core_promoted, S$n_core_demoted, nrow(prom), nrow(dem)))
}
if (!identical(as.integer(got_n), as.integer(expect_n))) {
  problems <- c(problems, sprintf(
    "discordance_network.csv describes %s gene(s) (promoted/demoted), current core sets are %d/%d",
    paste(got_n, collapse = "/"), expect_n[["promoted"]], expect_n[["demoted"]]))
}
orphan <- setdiff(nod_all$gene, d$gene[d$core != "none"])
if (length(orphan)) {
  problems <- c(problems, sprintf(
    "%d cached network node(s) are no longer core genes: %s",
    length(orphan), paste(utils::head(orphan, 8), collapse = ", ")))
}
if (length(problems)) {
  check_or_warn(paste0("STRING/GO cache is stale:\n  - ", paste(problems, collapse = "\n  - "),
       "\n  Refresh it:  make tables  ->  make discordance-network  ->  make tables"))
}

edg <- edg_all %>% filter(direction == "promoted")
nod <- nod_all %>%
  filter(direction == "promoted", degree > 0) %>%
  arrange(desc(degree), gene)

ang <- seq(pi / 2, -3 * pi / 2, length.out = nrow(nod) + 1)[seq_len(nrow(nod))]
nod$x <- cos(ang); nod$y <- sin(ang)
edg2 <- edg %>%
  inner_join(nod %>% select(a = gene, xa = x, ya = y), by = "a") %>%
  inner_join(nod %>% select(b = gene, xb = x, yb = y), by = "b")

pb <- ggplot() +
  geom_segment(data = edg2, aes(xa, ya, xend = xb, yend = yb, linewidth = score),
               colour = COL_FULL, alpha = 0.55, lineend = "round") +
  geom_point(data = nod, aes(x, y, size = degree), colour = COL_FULL) +
  geom_text(data = nod, aes(x * 1.42, y * 1.42, label = gene,
                            hjust = ifelse(x > 0.15, 0, ifelse(x < -0.15, 1, 0.5))),
            size = GEOM_PT(DENSE), fontface = "italic", colour = "black") +
  scale_linewidth_continuous(range = c(LW(HAIRLINE), LW(RULE)), guide = "none") +
  scale_size_continuous(range = c(1.2, 3.4), guide = "none") +
  coord_fixed(ratio = 1, xlim = c(-2.15, 2.15), ylim = c(-2.15, 2.15), clip = "off") +
  labs(x = NULL, y = NULL,
       title = sprintf("Promoted genes interact: %d edges vs %.1f\nexpected (%s)",
                       S$net_promoted_edges, S$net_promoted_null_mean,
                       ifelse(S$net_promoted_p < 0.001, "p < 0.001",
                              sprintf("p = %.3f", S$net_promoted_p))),
       subtitle = sprintf(
         "%d of %d connected; %d isolated.\n%d null draws. Demoted: %d edge, p = %.2f",
         S$net_promoted_connected, S$n_core_promoted,
         S$n_core_promoted - S$net_promoted_connected, S$n_null_network,
         S$net_demoted_edges, S$net_demoted_p)) +
  theme_void(base_size = AXIS_TITLE) +
  theme(text = element_text(family = font_family),
        plot.margin = margin(2, 2, 2, 2),
        plot.title = element_text(size = TITLE_SIZE, hjust = 0.5),
        plot.subtitle = element_text(size = DENSE, hjust = 0.5, colour = GREY),
        plot.tag = element_text(size = TAG_SIZE, face = "bold"))

# Panel c: identifies which features drive the discordance shift.
FEAT_LAB <- c(
  expected_n_regulators_residuals = "Regulatory burden (residual)",
  reg_burden_sig_Stim8hr = "Regulator burden (8 h stim)",
  reg_burden_sig_Stim48hr = "Regulator burden (48 h stim)",
  zscore_Th17 = "Th17 DE z-score", zscore_Th2 = "Th2 DE z-score",
  zscore_Th1 = "Th1 DE z-score", zscore_Treg = "Treg DE z-score",
  polar_rank_range = "Polarization rank range",
  polar_coef_rank_Rest = "Polarization coef (rest)",
  polar_coef_rank_Stim48hr = "Polarization coef (48 h)",
  n_sig_regulated_cytokine_receptors_Rest = "Cytokine receptors reg. (rest)",
  n_sig_regulated_cytokines_Rest = "Cytokines regulated (rest)",
  has_cytokine = "Cytokine measured",
  mis.z_score = "Missense constraint", lof.oe_ci.upper = "LoF tolerance",
  gwas_score = "GWAS score", IEI = "IEI gene",
  gene_burden_score = "Rare-variant burden")

ft <- read.csv(need("discordance_features.csv")) %>%
  filter(fdr < 0.05) %>%
  mutate(lab = ifelse(feature %in% names(FEAT_LAB), FEAT_LAB[feature], feature),
         dir_lab = factor(direction, levels = c("promoted", "demoted"),
                          labels = c("promoted", "demoted")),
         lab = paste0(lab, "  (", substr(direction, 1, 4), ")")) %>%
  arrange(direction == "promoted", rank_biserial) %>%
  mutate(lab = factor(lab, levels = lab))

pc <- ggplot(ft, aes(rank_biserial, lab, colour = dir_lab)) +
  geom_vline(xintercept = 0, colour = "black", linewidth = LW(HAIRLINE)) +
  geom_segment(aes(x = 0, xend = rank_biserial, yend = lab), linewidth = LW(HAIRLINE)) +
  geom_point(size = 1.3) +
  scale_colour_manual(values = c(promoted = COL_FULL, demoted = COL_GEN)) +
  scale_x_continuous(breaks = seq(-1, 1, 0.5)) +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(x = "Rank-Biserial Correlation\n(vs matched controls)", y = NULL,
       title = sprintf(
         "Every discriminating feature is\nfunctional-genomics (%d of %d, FDR<0.05)",
         sum(ft$block == "functional genomics"), nrow(ft))) +
  theme(axis.text.y = element_text(size = TINY),
        legend.position = "none")

# Panels a and b keep the cells they had in the former four-panel layout (8.4767 and 8.5614 cm, read from its
# gtable); panel c keeps its former 4.9378 cm panel width, centred between equal spacers.
bottom <- plot_spacer() + pc + plot_spacer() +
  plot_layout(widths = unit(c(1, 4.9378, 1), c("null", "cm", "null")))
fig <- (free(pa, side = "blr") + free(pb, side = "blr") + bottom) +
  plot_layout(design = c(area(1, 1), area(1, 2), area(2, 1, 2, 2)),
              widths = unit(c(8.4767, 8.5614), "cm"), heights = c(9.5, 8.0)) &
  theme(plot.tag = element_text(size = TAG_SIZE, face = "bold"))
fig <- fig + plot_annotation(tag_levels = "a")


OUT_PDF <- file.path(repo_root, "final_plots", "pdf",
                     "figure_discordance_genetics_vs_full.pdf")

FIG_W_CM <- 18.0
FIG_H_CM <- 17.5

save_figure(fig, OUT, OUT_PDF, FIG_W_CM, FIG_H_CM)

cat("wrote", file.path("final_plots", "supplementary", basename(OUT)), "\n")
cat(sprintf("   marked %d promoted / %d demoted | network %d edges p=%.3f | %d features\n",
            nrow(prom), nrow(dem), S$net_promoted_edges, S$net_promoted_p,
            nrow(ft)))
