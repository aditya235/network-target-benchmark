# ============================================================
# 02_build_assets.R -- everything that does not depend on the disease.
#
# The interactome, its topology, and the tractability labels are shared
# across all 34 diseases, so they are built once here. Only the
# association scores and drug labels vary per disease.
#
# Reuses files already downloaded for the earlier project when present,
# rather than re-fetching ~100 MB.
#
# Outputs: data/edges.csv, data/nodes_topology.csv, data/graph.rds,
#          data/tractability.csv, data/disease_assoc.csv
# ============================================================

suppressPackageStartupMessages({
  library(data.table); library(igraph); library(httr2); library(jsonlite)
  library(org.Hs.eg.db); library(AnnotationDbi)
})

OUT <- "data"; RAW <- "data/raw"
dir.create(RAW, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
if (!exists("%||%")) `%||%` <- function(x, y) if (is.null(x)) y else x   # base R < 4.4

OLD <- path.expand("~/t2d_gnn/data/raw")     # reuse earlier downloads
STRING_VERSION <- "12.0"; TAXON <- 9606; SCORE_MIN <- 700

OT <- "https://api.platform.opentargets.org/api/v4/graphql"
gql <- function(q, v = NULL, tries = 3) {
  for (k in seq_len(tries)) {
    r <- tryCatch(request(OT) |> req_body_json(list(query = q, variables = v)) |>
                    req_timeout(60) |> req_error(is_error = function(x) FALSE) |>
                    req_perform(), error = function(e) NULL)
    if (!is.null(r)) {
      out <- tryCatch(resp_body_json(r, simplifyVector = FALSE), error = function(e) NULL)
      if (!is.null(out$data)) return(out$data)
    }
    Sys.sleep(2 * k)
  }
  NULL
}

# ---- 1. STRING interactome ---------------------------------------------
get_string <- function(kind) {
  fn <- sprintf("%d.protein.%s.v%s.txt.gz", TAXON, kind, STRING_VERSION)
  dest <- file.path(RAW, fn)
  if (file.exists(dest)) return(dest)
  old <- file.path(OLD, fn)
  if (file.exists(old)) { msg("reusing ", fn); file.copy(old, dest); return(dest) }
  url <- sprintf("https://stringdb-downloads.org/download/protein.%s.v%s/%s",
                 kind, STRING_VERSION, fn)
  msg("downloading ", fn, " ...")
  utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  dest
}

if (!file.exists(file.path(OUT, "graph.rds"))) {
  links <- fread(get_string("links.detailed"))
  info  <- fread(get_string("info"))
  setnames(info, c("#string_protein_id", "preferred_name"), c("string_id", "symbol"),
           skip_absent = TRUE)
  links <- links[combined_score >= SCORE_MIN]
  map <- setNames(info$symbol, info$string_id)
  links[, `:=`(a = map[protein1], b = map[protein2])]
  links <- links[!is.na(a) & !is.na(b) & a != b]
  links[, `:=`(s1 = pmin(a, b), s2 = pmax(a, b))]
  edges <- unique(links[, .(source = s1, target = s2, score = combined_score)],
                  by = c("source", "target"))

  g <- graph_from_data_frame(edges[, .(source, target)], directed = FALSE)
  E(g)$weight <- edges$score / 1000
  comp <- components(g)
  g <- induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
  msg("graph: ", vcount(g), " nodes, ", ecount(g), " edges")

  msg("computing topology ...")
  nodes <- data.table(
    gene = V(g)$name, degree = degree(g), strength = strength(g),
    betweenness = betweenness(g, normalized = TRUE),
    closeness = closeness(g, normalized = TRUE),
    eigen = eigen_centrality(g)$vector, pagerank = page_rank(g)$vector,
    kcore = coreness(g),
    clustering = transitivity(g, type = "local", isolates = "zero"))
  adj <- as_adj_list(g)
  nodes[, mnc := vapply(seq_len(vcount(g)), function(i) {
    nb <- adj[[i]]; if (length(nb) < 2) return(length(nb))
    max(components(induced_subgraph(g, nb))$csize) }, numeric(1))]

  edges <- edges[source %in% V(g)$name & target %in% V(g)$name]
  fwrite(edges, file.path(OUT, "edges.csv"))
  fwrite(nodes, file.path(OUT, "nodes_topology.csv"))
  saveRDS(g, file.path(OUT, "graph.rds"))
} else {
  msg("graph already built; skipping")
  nodes <- fread(file.path(OUT, "nodes_topology.csv"))
}

# ---- 2. tractability: the disease-agnostic label set --------------------
# Curated, per-protein druggability. Unlike Open Targets' `clinical`
# association score it does not propagate a drug across its whole
# mechanism, so it is the clean label definition. Both are carried
# forward: comparing them across diseases is one of the study's aims.
tr_path <- file.path(OUT, "tractability.csv")
old_tr  <- path.expand("~/t2d_gnn/data/processed/tractability_raw.csv")
if (!file.exists(tr_path) && file.exists(old_tr)) {
  msg("reusing tractability from the earlier project")
  file.copy(old_tr, tr_path)
}
if (!file.exists(tr_path)) {
  msg("fetching tractability (this takes a few minutes) ...")
  map <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = unique(nodes$gene), keytype = "SYMBOL", columns = "ENSEMBL"))
  setDT(map); map <- unique(map[!is.na(ENSEMBL)], by = "SYMBOL")
  chunks <- split(map$ENSEMBL, ceiling(seq_len(nrow(map)) / 100))
  res <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    d <- gql('query($ids:[String!]!){ targets(ensemblIds:$ids){
               approvedSymbol tractability { label modality value } } }',
             list(ids = as.list(chunks[[i]])))
    tg <- d$targets; if (is.null(tg)) next
    res[[i]] <- rbindlist(lapply(tg, function(t) {
      tr <- t$tractability; if (is.null(tr) || !length(tr)) return(NULL)
      data.table(gene = t$approvedSymbol,
                 label_tr = vapply(tr, function(x) x$label %||% NA_character_, character(1)),
                 modality = vapply(tr, function(x) x$modality %||% NA_character_, character(1)),
                 value = vapply(tr, function(x) isTRUE(x$value), logical(1)))
    }), fill = TRUE)
    if (i %% 25 == 0) msg("  ", i, "/", length(chunks))
  }
  fwrite(rbindlist(res, fill = TRUE), tr_path)
}
tract <- fread(tr_path)
approved <- tract[label_tr == "Approved Drug" & value == TRUE, unique(gene)]
msg(length(approved), " proteins are approved-drug targets (disease-agnostic)")

# ---- 3. per-disease association and drug evidence ------------------------
feas <- fread(file.path(OUT, "disease_feasibility.csv"))
msg("fetching association data for ", nrow(feas), " diseases ...")

q_assoc <- '
query($id: String!, $idx: Int!, $size: Int!) {
  disease(efoId: $id) {
    associatedTargets(page: {index: $idx, size: $size}) {
      rows { target { approvedSymbol } datatypeScores { id score } } } } }'

fetch_disease <- function(efo) {
  acc <- list()
  for (i in 0:24) {
    d <- gql(q_assoc, list(id = efo, idx = i, size = 500L))
    rows <- d$disease$associatedTargets$rows
    if (is.null(rows) || !length(rows)) break
    acc[[length(acc) + 1L]] <- rbindlist(lapply(rows, function(r) {
      sc <- setNames(vapply(r$datatypeScores, function(x) x$score, numeric(1)),
                     vapply(r$datatypeScores, function(x) x$id, character(1)))
      g <- function(k) unname(if (k %in% names(sc)) sc[[k]] else 0)
      data.table(gene = r$target$approvedSymbol,
                 genetic = g("genetic_association"),
                 clinical = g("clinical"),
                 animal = g("animal_model"),
                 somatic = g("somatic_mutation"),
                 literature = g("literature"))
    }))
    if (length(rows) < 500L) break
  }
  if (!length(acc)) return(NULL)
  unique(rbindlist(acc, fill = TRUE), by = "gene")
}

all_dis <- rbindlist(lapply(seq_len(nrow(feas)), function(i) {
  tb <- fetch_disease(feas$efo_id[i])
  if (is.null(tb)) return(NULL)
  tb[, `:=`(disease = feas$query[i], efo_id = feas$efo_id[i])]
  cat(sprintf("  %-38s %5d genes\n", substr(feas$query[i], 1, 37), nrow(tb)))
  tb
}), fill = TRUE)

fwrite(all_dis, file.path(OUT, "disease_assoc.csv"))
msg(uniqueN(all_dis$disease), " diseases, ",
    nrow(all_dis), " disease-gene rows written")

# ---- 4. summary ----------------------------------------------------------
in_graph <- all_dis[gene %in% nodes$gene]
summ <- in_graph[, .(
  n_in_graph   = .N,
  n_clinical   = sum(clinical > 0),
  n_genetic    = sum(genetic > 0),
  n_approved   = sum(gene %in% approved)
), by = disease][order(-n_clinical)]
fwrite(summ, file.path(OUT, "disease_summary.csv"))
print(summ)

msg("02_build_assets.R complete")
