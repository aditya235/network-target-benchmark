# ============================================================
# 03_figures.R -- publication figures.
#
# Two claims drive the figure set:
#   (1) hub ranking is degree       -> Fig 2, Fig 3
#   (2) most of the signal is general druggability, not disease
#                                   -> Fig 4, Fig 5
#
# Outputs: results/figures/*.pdf and *.png
# ============================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

RES <- "results"; FIG <- file.path(RES, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

# OUP asks for high-resolution .eps or .tif. PDF and PNG are kept for
# drafting; TIFF (LZW-compressed, 600 dpi) and EPS are the submission
# formats.
save_fig <- function(p, name, w, h) {
  for (ext in c("pdf", "png"))
    ggsave(file.path(FIG, paste0(name, ".", ext)), p, width = w, height = h, dpi = 300)
  ggsave(file.path(FIG, paste0(name, ".tiff")), p, width = w, height = h,
         dpi = 600, device = "tiff", compression = "lzw")
  ok <- tryCatch({
    ggsave(file.path(FIG, paste0(name, ".eps")), p, width = w, height = h,
           device = cairo_ps, fallback_resolution = 600); TRUE
  }, error = function(e) FALSE)
  msg("wrote ", name, if (!ok) "  (eps failed; tiff written)" else "")
}

TOPO_M <- c("Degree", "MNC (CytoHubba)", "Betweenness")
GRAPH_M <- c("spectral + XGBoost", "RWR", "GAT")
pal_fam <- c("Topology ranking" = "#2471A3",
             "Graph learning"   = "#C0392B",
             "No graph"         = "#E67E22")
fam_of <- function(m) fifelse(m %in% TOPO_M, "Topology ranking",
                       fifelse(m == "XGBoost (no graph)", "No graph", "Graph learning"))

base_thm <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 11))

# ============================================================
# Figure 2 -- the degree control
# ============================================================
bm <- fread(file.path(RES, "benchmark_raw.csv"))
bm[, family := fam_of(method)]

# Never hard-code the disease count: one pre-specified disease was dropped
# when its identifier returned no data, and captions saying "34" survived
# the change silently.
ND <- uniqueN(bm$disease)
msg("figures built on ", ND, " diseases")

# per disease x method x evaluation, averaged over folds
per_dis <- bm[, .(AUROC = mean(AUROC)), by = .(disease, method, evaluation, family)]
ord <- per_dis[evaluation == "degree_matched", .(m = mean(AUROC)),
               by = method][order(-m), method]
per_dis[, method := factor(method, levels = ord)]
per_dis[, evaluation := factor(evaluation, c("all_nodes", "degree_matched"),
                               c("All nodes", "Degree-matched"))]

pA <- ggplot(per_dis, aes(method, AUROC, fill = evaluation)) +
  geom_hline(yintercept = .5, linetype = "dashed", colour = "grey40") +
  geom_boxplot(outlier.size = .5, alpha = .9, width = .7,
               position = position_dodge(.8)) +
  scale_fill_manual(values = c("All nodes" = "#BFC9CA",
                               "Degree-matched" = "#C0392B")) +
  coord_flip(ylim = c(0.45, 1)) +
  labs(x = NULL, y = "AUROC", fill = NULL,
       title = "A  Performance before and after controlling for degree",
       subtitle = sprintf("Each point is one of %d diseases", ND)) +
  base_thm

# share of above-chance performance explained by degree
eff <- per_dis[, .(AUROC = mean(AUROC)), by = .(method, evaluation)]
eff <- dcast(eff, method ~ evaluation, value.var = "AUROC")
setnames(eff, c("All nodes", "Degree-matched"), c("all_nodes", "matched"))
eff[, pct_degree := 100 * (all_nodes - matched) / (all_nodes - 0.5)]
eff[, family := fam_of(as.character(method))]
setorder(eff, pct_degree)
eff[, method := factor(method, levels = method)]

pB <- ggplot(eff, aes(method, pct_degree, fill = family)) +
  geom_col(width = .68) +
  geom_text(aes(label = sprintf("%.0f%%", pct_degree)), hjust = -0.15, size = 3) +
  scale_fill_manual(values = pal_fam) +
  coord_flip(ylim = c(0, 115)) +
  labs(x = NULL, y = "% of above-chance performance explained by degree",
       fill = NULL, title = "B  How much was connectivity?") +
  base_thm

save_fig(pA, "fig2a_degree_control", 7.5, 4.5)
save_fig(pB, "fig2b_degree_share", 7, 4)
save_fig(pA / pB + plot_layout(heights = c(1.15, 1)), "fig2_degree", 8, 8.5)

# ============================================================
# Figure 3 -- is the collapse consistent across diseases?
# ============================================================
dm <- per_dis[evaluation == "Degree-matched"]
dm[, disease := factor(disease, levels = dm[method == ord[1]][order(AUROC), disease])]

pF3 <- ggplot(dm, aes(method, disease, fill = AUROC)) +
  geom_tile(colour = "white", linewidth = .3) +
  scale_fill_gradient2(low = "#2471A3", mid = "#F7F7F7", high = "#C0392B",
                       midpoint = .5, limits = c(.4, 1), oob = squish) +
  labs(x = NULL, y = NULL, fill = "AUROC",
       title = sprintf("Degree-matched performance, all %d diseases", ND),
       subtitle = sprintf("Hub-ranking methods sit at chance in all %d", ND)) +
  base_thm +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = element_text(size = 7))
save_fig(pF3, "fig3_per_disease", 7.5, 8)

# ============================================================
# Figure 4 -- how much of the signal is disease-specific?
# ============================================================
ag_file <- file.path(RES, "agnostic_comparison.csv")
if (file.exists(ag_file)) {
  ag <- fread(ag_file)
  long <- melt(ag, id.vars = c("disease", "n_labels", "pct_labels_also_approved"),
               measure.vars = c("own_cv_auroc", "agnostic_auroc", "other_disease_auroc"),
               variable.name = "model", value.name = "AUROC")
  long[, model := factor(model,
      levels = c("own_cv_auroc", "agnostic_auroc", "other_disease_auroc"),
      labels = c("Own disease", "No disease info", "Another disease"))]

  pA4 <- ggplot(long, aes(model, AUROC, group = disease)) +
    geom_line(alpha = .25, colour = "grey45") +
    geom_point(aes(colour = model), size = 1.8, alpha = .85) +
    stat_summary(aes(group = 1), fun = median, geom = "line",
                 colour = "#C0392B", linewidth = 1.1) +
    stat_summary(aes(group = 1), fun = median, geom = "point",
                 colour = "#C0392B", size = 3.2) +
    scale_colour_manual(values = c("Own disease" = "#1E8449",
                                   "No disease info" = "#7F8C8D",
                                   "Another disease" = "#2471A3"),
                        guide = "none") +
    labs(x = NULL, y = "Degree-matched AUROC",
         title = "A  What does knowing the disease buy you?",
         subtitle = "One line per disease; red is the median") +
    base_thm

  pB4 <- ggplot(ag, aes(pct_labels_also_approved)) +
    geom_histogram(bins = 18, fill = "#2471A3", colour = "white", alpha = .9) +
    geom_vline(xintercept = median(ag$pct_labels_also_approved),
               linetype = "dashed", colour = "#C0392B", linewidth = .8) +
    annotate("text", x = median(ag$pct_labels_also_approved) - 1,
             y = Inf, vjust = 1.8, hjust = 1, size = 3, colour = "#C0392B",
             label = sprintf("median %.1f%%", median(ag$pct_labels_also_approved))) +
    labs(x = "% of a disease's drug-target labels that are\nalso general approved-drug targets",
         y = "Diseases",
         title = "B  Why: the labels are mostly not disease-specific") +
    base_thm

  save_fig(pA4 + pB4 + plot_layout(widths = c(1.2, 1)), "fig4_disease_specificity", 11, 4.8)

  ag[, adv := own_cv_auroc - agnostic_auroc]
  setorder(ag, adv)
  ag[, disease := factor(disease, levels = disease)]
  pC4 <- ggplot(ag, aes(disease, adv, fill = adv > 0.02)) +
    geom_col(width = .7) +
    geom_hline(yintercept = 0, colour = "grey30") +
    scale_fill_manual(values = c("TRUE" = "#1E8449", "FALSE" = "#BFC9CA"),
                      labels = c("= 0.02", "> 0.02"), name = "Gain") +
    coord_flip() +
    labs(x = NULL, y = "AUROC gain from disease-specific information",
         title = "Disease information helps, but modestly",
         subtitle = sprintf("median +%.3f; %d of %d diseases gain > 0.02",
                            median(ag$adv), sum(ag$adv > 0.02), nrow(ag))) +
    base_thm + theme(axis.text.y = element_text(size = 7))
  save_fig(pC4, "fig5_disease_advantage", 7, 7)
}

# ============================================================
# Figure 6 -- cross-disease prediction similarity
# ============================================================
cc <- file.path(RES, "pred_correlation_spectral.csv")
if (file.exists(cc)) {
  M <- as.matrix(fread(cc), rownames = 1)
  dd <- as.data.table(as.table(M)); setnames(dd, c("d1", "d2", "rho"))
  o <- rownames(M)[hclust(as.dist(1 - M))$order]
  dd[, `:=`(d1 = factor(d1, o), d2 = factor(d2, o))]
  pF6 <- ggplot(dd, aes(d1, d2, fill = rho)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2471A3", mid = "white", high = "#C0392B",
                         midpoint = .5, limits = c(0, 1)) +
    labs(x = NULL, y = NULL, fill = "Spearman",
         title = "Prediction similarity between diseases",
         subtitle = "Rankings are correlated but not interchangeable") +
    base_thm +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6))
  save_fig(pF6, "fig6_cross_disease", 8.5, 8)
}

msg("03_figures.R complete -> ", FIG)
