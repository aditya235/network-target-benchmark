# ============================================================
# 00_figure1_design.R -- Figure 1, study design schematic.
#
# Standalone: needs no pipeline output, so it can be run before or
# after the benchmark. Counts are read from data/ when those files
# exist and fall back to the values reported in the manuscript.
#
# Output: results/figures/fig1_design.{pdf,png,tiff,eps}
# ============================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

RES <- "results"; FIG <- file.path(RES, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

# ---- counts, read from disk where available ----------------
n_nodes <- 15882L; n_edges <- 236712L; n_dis <- 33L; n_appr <- 946L
if (file.exists("data/nodes_topology.csv"))
  n_nodes <- nrow(fread("data/nodes_topology.csv"))
if (file.exists("data/edges.csv"))
  n_edges <- nrow(fread("data/edges.csv"))
if (file.exists("data/disease_summary.csv"))
  n_dis <- nrow(fread("data/disease_summary.csv"))
if (file.exists("data/tractability.csv")) {
  tr <- fread("data/tractability.csv")
  if ("approved_drug" %in% names(tr)) n_appr <- sum(tr$approved_drug == TRUE)
}
msg("figure 1 on ", n_nodes, " nodes / ", n_edges, " edges / ", n_dis, " diseases")

# ---- boxes -------------------------------------------------
# One row per box. Keep every vector the same length; the schematic is
# edited by hand often enough that a mismatched column is the likeliest
# way to break this script.
boxes <- data.table(
  id = c("string", "ot", "graph", "topo", "labels",
         "cent", "diff", "ml", "evals", "agno", "res"),
  x  = c(1.0, 1.0, 3.3, 3.3, 3.3, 5.7, 5.7, 5.7, 8.1, 8.1, 10.4),
  y  = c(7.6, 3.2, 7.6, 5.2, 2.8, 8.0, 6.3, 4.2, 6.8, 2.8, 5.0),
  w  = c(1.6, 1.6, 1.8, 1.8, 1.8, 1.9, 1.9, 1.9, 2.0, 2.0, 1.9),
  h  = c(1.5, 1.4, 1.6, 1.4, 2.4, 1.2, 1.2, 1.8, 2.8, 2.4, 2.8),
  kind = c("source", "source", "asset", "asset", "asset",
           "method", "method", "method", "eval", "eval", "out"),
  lab = c(
    # Plain ASCII throughout. The pdf and eps devices here convert
    # multibyte characters through mbcsToSbcs and silently substitute a
    # dot, so an earlier version rendered "combined score . 700".
    "STRING v12.0\ncombined score\n700 or higher",
    "Open Targets\nPlatform",
    sprintf("Interactome\n%s proteins\n%s edges",
            format(n_nodes, big.mark = ","), format(n_edges, big.mark = ",")),
    "9 topology\nfeatures\ncomputed once",
    sprintf("%d diseases\n\nlabels:\nclinical evidence\n\nfeatures: genetic,\nanimal, somatic", n_dis),
    "Centrality\ndegree, MNC,\nbetweenness",
    "Diffusion\nrandom walk\nwith restart",
    "Machine learning\nXGBoost,\nspectral + XGBoost,\ngraph attention",
    "Every method\nevaluated twice\n\n1. All nodes\n    (conventional)\n\n2. Degree-matched\n    3 negatives per\n    positive, same\n    degree stratum",
    sprintf("Disease-agnostic\nbaseline\n\n%d approved-drug\ntargets\n\ntopology features\nonly, no disease\ninput of any kind", n_appr),
    "Per disease:\nAUROC under\nboth schemes\n\nhow much of the\nsignal is degree\n\nhow much is\ndisease-specific")
)
stopifnot(length(unique(sapply(boxes, length))) == 1L)
boxes[, `:=`(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2)]
boxes[, light := kind %in% c("method", "eval", "out")]

# ---- arrows ------------------------------------------------
# Topology is derived from the graph, not from STRING directly, and
# Open Targets feeds the labels only. An earlier draft drew Open Targets
# into the interactome, which it does not touch.
# No topo -> agno arrow. Routed straight it passes behind the machine
# learning box and reads as though the models feed the agnostic baseline,
# which is the one input that baseline must not have. The agno box states
# its own inputs instead.
arr <- data.table(
  from = c("string", "graph", "ot",
           "topo", "graph", "labels", "topo", "labels", "graph",
           "cent", "diff", "ml",
           "evals", "agno"),
  to   = c("graph", "topo", "labels",
           "cent", "diff", "diff", "ml", "ml", "ml",
           "evals", "evals", "evals",
           "res", "res")
)
# labels -> diff is not decoration: random walk with restart is seeded
# from the training positives, so it consumes labels. Centrality is the
# only method family here that uses none, and the figure has to show
# that difference because it is what the paper measures.
stopifnot(nrow(arr) == length(arr$from), all(c(arr$from, arr$to) %in% boxes$id))

# Boxes in the same column connect edge to edge vertically. Attaching
# every arrow right-to-left draws same-column links as a backwards loop.
seg <- merge(arr, boxes[, .(from = id, fx = x, fy = y,
                            fxmax = xmax, fymin = ymin, fymax = ymax)],
             by = "from", sort = FALSE)
seg <- merge(seg, boxes[, .(to = id, tx = x, ty = y,
                            txmin = xmin, tymin = ymin, tymax = ymax)],
             by = "to", sort = FALSE)
seg[, same_col := abs(fx - tx) < 0.5]
seg[, `:=`(
  x1 = fifelse(same_col, fx, fxmax),
  y1 = fifelse(same_col, fifelse(fy > ty, fymin, fymax), fy),
  x2 = fifelse(same_col, tx, txmin),
  y2 = fifelse(same_col, fifelse(fy > ty, tymax, tymin), ty)
)]

note <- data.table(
  x = 3.3, y = 1.1,
  lab = "34 diseases pre-specified; multiple sclerosis excluded\n(no generic Open Targets entry carrying association data)")

p <- ggplot() +
  geom_segment(data = seg, aes(x = x1, y = y1, xend = x2, yend = y2),
               colour = "grey45", linewidth = .4,
               arrow = arrow(length = unit(0.14, "cm"), type = "closed")) +
  geom_rect(data = boxes,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = kind),
            colour = "grey25", linewidth = .4, alpha = .95) +
  geom_text(data = boxes[light == TRUE],  aes(x = x, y = y, label = lab),
            size = 2.6, lineheight = 1.1, colour = "white") +
  geom_text(data = boxes[light == FALSE], aes(x = x, y = y, label = lab),
            size = 2.6, lineheight = 1.1, colour = "grey10") +
  geom_text(data = note, aes(x = x, y = y, label = lab),
            size = 2.2, colour = "grey30", fontface = "italic", lineheight = 1.1) +
  scale_fill_manual(values = c(source = "#D5DBDB", asset = "#AEB6BF",
                               method = "#2471A3", eval = "#C0392B",
                               out = "#1E8449"), guide = "none") +
  annotate("text", x = c(1.0, 3.3, 5.7, 8.1, 10.4), y = 9.4,
           label = c("Sources", "Shared assets", "Methods", "Evaluation", "Outputs"),
           fontface = "bold", size = 3.1, colour = "grey20") +
  coord_cartesian(xlim = c(0, 11.5), ylim = c(0.3, 9.7)) +
  theme_void() +
  theme(plot.margin = margin(6, 6, 6, 6))

# cairo_pdf rather than the default pdf device: it embeds fonts properly
# and does not route text through mbcsToSbcs.
ggsave(file.path(FIG, "fig1_design.pdf"), p, width = 10.5, height = 6.4,
       device = cairo_pdf)
ggsave(file.path(FIG, "fig1_design.png"), p, width = 10.5, height = 6.4, dpi = 300)
ggsave(file.path(FIG, "fig1_design.tiff"), p, width = 10.5, height = 6.4,
       dpi = 600, device = "tiff", compression = "lzw")
tryCatch(
  ggsave(file.path(FIG, "fig1_design.eps"), p, width = 10.5, height = 6.4,
         device = cairo_ps, fallback_resolution = 600),
  error = function(e) msg("eps failed; tiff written"))

msg("00_figure1_design.R complete -> ", FIG)
