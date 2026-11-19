# ============================================================
# 05_figure6_subnetwork.R -- Figure 6, centrality inside
# disease-specific subnetworks.
#
# Reads results/subnetwork_centrality.csv, written by
# python/subnetwork_centrality.py.
#
# Panel A: per-disease matched AUROC by method and arm
# Panel B: share of above-chance performance taken by connectivity
#
# Output: results/figures/fig6_subnetwork.{pdf,png,tiff,eps}
# ============================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

RES <- "results"; FIG <- file.path(RES, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

r <- fread(file.path(RES, "subnetwork_centrality.csv"))
r <- r[method != "Global degree"]
r[, arm := factor(fifelse(arm == "seeds", "Seed genes only", "Plus first neighbours"),
                  levels = c("Seed genes only", "Plus first neighbours"))]

# Order methods by matched performance in the seed-only arm, so the
# reading order matches the argument rather than the alphabet.
ord <- r[arm == "Seed genes only",
         .(m = median(AUROC_matched, na.rm = TRUE)), by = method][order(-m), method]
r[, method := factor(method, levels = ord)]

# Arms do not share a denominator: seed-only subnetworks are evaluable
# only where enough labelled positives survive the gene selection, so a
# single n would overstate the sample behind the green boxes.
narm <- r[!is.na(AUROC_matched), .(n = uniqueN(disease)), by = arm]
nlab <- paste(sprintf("%s n = %d", narm$arm, narm$n), collapse = "; ")
pal <- c("Seed genes only" = "#1D9E75", "Plus first neighbours" = "#D85A30")

base_thm <- theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 10))

pA <- ggplot(r[!is.na(AUROC_matched)],
             aes(method, AUROC_matched, colour = arm, fill = arm)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey35", linewidth = .4) +
  geom_boxplot(width = .6, alpha = .25, outlier.shape = NA,
               position = position_dodge(preserve = "single")) +
  geom_point(position = position_jitterdodge(jitter.width = .12, dodge.width = .6),
             size = .9, alpha = .6) +
  scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
  coord_flip() +
  labs(x = NULL, y = "Degree-matched AUROC", title = "A  What survives the control",
       subtitle = sprintf("one point per disease; dashed line is chance (%s)", nlab)) +
  base_thm

# Share is a pooled quantity, so one bar per method and arm rather than a
# distribution. Capped at 100 for display: above-100 values mean matched
# performance dipped below chance, which is noise around zero, not more
# than total attribution. The cap is stated in the caption.
sh <- r[, .(all = mean(AUROC_subnetwork, na.rm = TRUE),
            mat = mean(AUROC_matched, na.rm = TRUE)), by = .(arm, method)]
sh[, share := 100 * (all - mat) / (all - 0.5)]
sh[, capped := share > 100]
sh[, share_plot := pmin(share, 100)]

pB <- ggplot(sh, aes(method, share_plot, fill = arm)) +
  geom_col(width = .62, position = position_dodge(preserve = "single"), alpha = .9) +
  geom_text(aes(label = sprintf("%.0f%%%s", share, fifelse(capped, "*", ""))),
            position = position_dodge(width = .62), hjust = 1.18,
            size = 2.5, colour = "white", fontface = "bold") +
  scale_fill_manual(values = pal) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, .02))) +
  coord_flip() +
  labs(x = NULL, y = "% of above-chance performance explained by degree",
       title = "B  How much was connectivity",
       subtitle = "* bar capped at 100%; matched value fell just below chance") +
  base_thm

p <- pA / pB + plot_layout(heights = c(1.15, 1))

save_fig <- function(name, plot, w, h) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), plot, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIG, paste0(name, ".png")), plot, width = w, height = h, dpi = 300)
  ggsave(file.path(FIG, paste0(name, ".tiff")), plot, width = w, height = h,
         dpi = 600, device = "tiff", compression = "lzw")
  tryCatch(ggsave(file.path(FIG, paste0(name, ".eps")), plot, width = w, height = h,
                  device = cairo_ps, fallback_resolution = 600),
           error = function(e) msg("eps failed for ", name, "; tiff written"))
  msg("wrote ", name)
}

save_fig("fig6_subnetwork", p, 7.2, 7.4)
msg("05_figure6_subnetwork.R complete -> ", FIG)
