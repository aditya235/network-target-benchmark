# ============================================================
# 04_supplementary.R -- supplementary tables, submission-ready.
#
# S1  the 34 diseases, pre-specified, with resolved identifiers and
#     label counts (the table a reviewer checks first)
# S2  full per-method per-evaluation results
# S3  per-disease degree effect
# S4  disease-specificity comparison
# S5  cross-disease prediction similarity
#
# Written as CSV (editable, per OUP's requirement that tables not be
# images) plus a single Markdown file for pasting into the supplement.
# ============================================================

suppressPackageStartupMessages({ library(data.table) })

RES <- "results"; SUP <- file.path(RES, "supplementary")
dir.create(SUP, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

md <- c("# Supplementary Material",
        "",
        "**Network methods find druggable proteins, not disease genes:",
        "a controlled benchmark across 34 diseases**",
        "")

add <- function(dt, id, title, note = NULL) {
  fwrite(dt, file.path(SUP, sprintf("Table_%s.csv", id)))
  md <<- c(md, sprintf("## Supplementary Table %s. %s", id, title), "")
  if (!is.null(note)) md <<- c(md, note, "")
  hdr <- paste(names(dt), collapse = " | ")
  sep <- paste(rep("---", ncol(dt)), collapse = " | ")
  rows <- apply(dt, 1, function(r) paste(r, collapse = " | "))
  md <<- c(md, paste0("| ", hdr, " |"), paste0("| ", sep, " |"),
           paste0("| ", rows, " |"), "")
  msg("Table ", id, ": ", nrow(dt), " rows")
}

# ---- S1. the disease panel ----------------------------------------------
feas <- fread(file.path("data", "disease_feasibility.csv"))
dsum_f <- file.path("data", "disease_summary.csv")
if (file.exists(dsum_f)) {
  dsum <- fread(dsum_f)
  feas <- merge(feas, dsum, by.x = "query", by.y = "disease", all.x = TRUE,
                suffixes = c("", ".g"))
}
s1 <- feas[, .(
  Disease            = query,
  `Open Targets ID`  = efo_id,
  `Matched term`     = matched_name,
  `Associated targets` = n_associated,
  `In network`       = if ("n_in_graph" %in% names(feas)) n_in_graph else NA_integer_,
  `Positive labels`  = if ("n_clinical.g" %in% names(feas)) n_clinical.g else n_clinical,
  `With genetic evidence` = if ("n_genetic.g" %in% names(feas)) n_genetic.g else n_genetic
)]
setorder(s1, Disease)
add(s1, "S1", "Disease panel",
    paste("Thirty-four diseases, pre-specified before any analysis.",
          "Identifiers were resolved programmatically through the Open Targets",
          "search API. Positive labels are proteins with non-zero `clinical`",
          "evidence for that disease that are present in the interaction network."))

# ---- S2. full benchmark results ------------------------------------------
bm <- fread(file.path(RES, "benchmark_raw.csv"))
mcols <- intersect(c("AUROC", "AUPRC"), names(bm))
pcol <- grep("^Precision@", names(bm), value = TRUE)
s2 <- bm[, c(lapply(.SD, function(x) round(mean(x), 4)),
             lapply(.SD, function(x) round(sd(x), 4))),
         by = .(Evaluation = evaluation, Method = method),
         .SDcols = c(mcols, pcol)]
setnames(s2, 3:ncol(s2),
         c(paste0(c(mcols, pcol), "_mean"), paste0(c(mcols, pcol), "_sd")))
setorder(s2, Evaluation, -AUROC_mean)
add(s2, "S2", "Benchmark results by method and evaluation scheme",
    "Mean and standard deviation across 34 diseases x 3 repeats x 5 folds.")

# ---- S3. per-disease degree effect ----------------------------------------
pd <- bm[, .(AUROC = mean(AUROC)), by = .(disease, method, evaluation)]
w <- dcast(pd, disease + method ~ evaluation, value.var = "AUROC")
if (all(c("all_nodes", "degree_matched") %in% names(w))) {
  w[, degree_share := round(100 * (all_nodes - degree_matched) /
                              pmax(all_nodes - 0.5, 1e-9), 1)]
  s3 <- w[, .(Disease = disease, Method = method,
              `AUROC all nodes` = round(all_nodes, 4),
              `AUROC degree-matched` = round(degree_matched, 4),
              `% explained by degree` = degree_share)]
  setorder(s3, Disease, Method)
  add(s3, "S3", "Degree effect, every disease and method",
      paste("Percentage of above-chance performance removed by degree",
            "matching, calculated as (AUROC_all - AUROC_matched) /",
            "(AUROC_all - 0.5)."))
}

# ---- S4. disease specificity -----------------------------------------------
ag_f <- file.path(RES, "agnostic_comparison.csv")
if (file.exists(ag_f)) {
  ag <- fread(ag_f)
  s4 <- ag[, .(
    Disease = disease,
    `Positive labels` = n_labels,
    `% labels also general drug targets` = pct_labels_also_approved,
    `Own disease model` = own_cv_auroc,
    `No disease information` = agnostic_auroc,
    `Another disease's model` = other_disease_auroc,
    `Gain from disease information` = disease_advantage
  )]
  setorder(s4, -`Gain from disease information`)
  add(s4, "S4", "Contribution of disease-specific information",
      paste("All values are degree-matched AUROC. The disease-agnostic model",
            "was trained only to recognise targets of approved drugs, using",
            "topology features, with no disease input of any kind."))
}

# ---- S5. cross-disease similarity -------------------------------------------
cs_f <- file.path(RES, "cross_disease_summary.json")
if (file.exists(cs_f)) {
  js <- jsonlite::fromJSON(cs_f)
  s5 <- rbindlist(lapply(names(js), function(m) {
    v <- js[[m]]
    data.table(Model = toupper(m),
               `Median rank correlation` = round(v$pred_rho_median, 3),
               `Top-100 prediction overlap` = round(v[[grep("^top", names(v))[1]]], 3),
               `Label set overlap` = round(v$label_jaccard_median, 3),
               `Cross-prediction AUROC` = round(v$cross_auroc_other, 3))
  }))
  add(s5, "S5", "Similarity of predictions between diseases",
      paste("Overlap is the Jaccard index. Prediction overlap should be",
            "compared with label overlap, not with the rank correlation,",
            "which is a different scale."))
}

writeLines(md, file.path(SUP, "Supplementary_Material.md"))
msg("wrote ", file.path(SUP, "Supplementary_Material.md"))
msg("04_supplementary.R complete")
