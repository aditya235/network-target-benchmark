# ============================================================
# 01_select_diseases.R -- which diseases can this study actually use?
#
# Gating step. Everything downstream depends on how many diseases have
# enough drug-target labels to train and evaluate on. Run this before
# writing any more pipeline.
#
# Disease IDs are resolved from NAMES via the Open Targets search API
# rather than hard-coded, because IDs recalled from memory are exactly
# how the previous project lost an afternoon.
#
# Outputs: data/disease_feasibility.csv
# ============================================================

suppressPackageStartupMessages({
  library(data.table); library(httr2); library(jsonlite)
})

OT   <- "https://api.platform.opentargets.org/api/v4/graphql"
OUT  <- "data"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
msg  <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

gql <- function(q, v = NULL, quiet = TRUE) {
  r <- tryCatch(
    request(OT) |> req_body_json(list(query = q, variables = v)) |>
      req_timeout(60) |> req_error(is_error = function(x) FALSE) |> req_perform(),
    error = function(e) NULL)
  if (is.null(r)) return(NULL)
  out <- tryCatch(resp_body_json(r, simplifyVector = FALSE), error = function(e) NULL)
  if (!is.null(out$errors) && !quiet)
    message("OT: ", out$errors[[1]]$message)
  out$data
}

# ---- 0. confirm the search field's shape before relying on it ---------
msg("checking the search API shape ...")
probe_shapes <- list(
  'query($q:String!){ search(queryString:$q, entityNames:["disease"], page:{index:0,size:3}){
     hits { id name entity } } }',
  'query($q:String!){ search(queryString:$q, entityNames:["disease"]){
     hits { id name entity } } }',
  'query($q:String!){ search(queryString:$q){ hits { id name entity } } }'
)
SEARCH_Q <- NULL
for (q in probe_shapes) {
  d <- gql(q, list(q = "type 2 diabetes"))
  if (!is.null(d$search$hits)) { SEARCH_Q <- q; break }
}
if (is.null(SEARCH_Q))
  stop("Could not find a working search query shape. Introspect `Query.search` ",
       "and update probe_shapes.")
msg("search API OK")

resolve_disease <- function(name) {
  d <- gql(SEARCH_Q, list(q = name))
  hits <- d$search$hits
  if (is.null(hits) || !length(hits)) return(NULL)
  h <- Filter(function(x) identical(x$entity, "disease"), hits)
  if (!length(h)) h <- hits
  list(id = h[[1]]$id, matched_name = h[[1]]$name)
}

# ---- 1. candidate diseases ---------------------------------------------
# Common, non-cancerous conditions spanning therapeutic areas. Deliberately
# broader than needed: the label counts below decide what survives.
# PRE-SPECIFY this list and keep it under version control -- selecting
# diseases after seeing results is the fastest way to lose a reviewer.
candidates <- c(
  "type 2 diabetes mellitus", "type 1 diabetes mellitus", "obesity",
  "hypertension", "coronary artery disease", "heart failure", "stroke",
  "atrial fibrillation", "hypercholesterolemia",
  "asthma", "chronic obstructive pulmonary disease",
  "rheumatoid arthritis", "systemic lupus erythematosus",
  "inflammatory bowel disease", "Crohn's disease", "ulcerative colitis",
  "psoriasis", "multiple sclerosis",
  "Alzheimer's disease", "Parkinson's disease", "epilepsy", "migraine",
  "major depressive disorder", "schizophrenia", "bipolar disorder",
  "osteoporosis", "osteoarthritis", "gout",
  "chronic kidney disease", "non-alcoholic fatty liver disease",
  "glaucoma", "age-related macular degeneration",
  "atopic dermatitis", "allergic rhinitis"
)

msg(length(candidates), " candidate diseases")

# ---- 2. resolve IDs ------------------------------------------------------
msg("resolving disease IDs ...")
res <- rbindlist(lapply(candidates, function(nm) {
  r <- resolve_disease(nm)
  Sys.sleep(0.15)
  if (is.null(r)) {
    message("  unresolved: ", nm)
    return(data.table(query = nm, efo_id = NA_character_, matched_name = NA_character_))
  }
  data.table(query = nm, efo_id = r$id, matched_name = r$matched_name)
}), fill = TRUE)
res <- res[!is.na(efo_id)]
msg(nrow(res), " resolved")

# ---- 3. count usable labels per disease ---------------------------------
# Two label definitions, deliberately compared:
#   clinical  -- Open Targets drug evidence. Propagates a drug across its
#                whole mechanism (metformin -> every complex I subunit),
#                so it inflates counts and corrupts the positive class.
#   strong    -- the same, thresholded, as a crude contamination check.
q_assoc <- '
query($id: String!, $idx: Int!, $size: Int!) {
  disease(efoId: $id) {
    associatedTargets(page: {index: $idx, size: $size}) {
      count
      rows { target { approvedSymbol } datatypeScores { id score } }
    }
  }
}'

count_labels <- function(efo) {
  rows_all <- list()
  for (i in 0:19) {
    d <- gql(q_assoc, list(id = efo, idx = i, size = 500L))
    rows <- d$disease$associatedTargets$rows
    if (is.null(rows) || !length(rows)) break
    rows_all[[length(rows_all) + 1L]] <- rbindlist(lapply(rows, function(r) {
      sc <- setNames(vapply(r$datatypeScores, function(x) x$score, numeric(1)),
                     vapply(r$datatypeScores, function(x) x$id, character(1)))
      g <- function(k) unname(if (k %in% names(sc)) sc[[k]] else 0)
      data.table(gene = r$target$approvedSymbol,
                 genetic = g("genetic_association"),
                 clinical = g("clinical"))
    }))
    if (length(rows) < 500L) break
  }
  if (!length(rows_all)) return(NULL)
  unique(rbindlist(rows_all, fill = TRUE), by = "gene")
}

msg("counting labels per disease (this is the slow part) ...")
feas <- rbindlist(lapply(seq_len(nrow(res)), function(i) {
  tb <- count_labels(res$efo_id[i])
  if (is.null(tb)) return(NULL)
  dt <- data.table(
    query        = res$query[i],
    efo_id       = res$efo_id[i],
    matched_name = res$matched_name[i],
    n_associated = nrow(tb),
    n_genetic    = sum(tb$genetic > 0),
    n_clinical   = sum(tb$clinical > 0),
    n_clin_0.5   = sum(tb$clinical > 0.5),
    n_clin_0.9   = sum(tb$clinical > 0.9)
  )
  cat(sprintf("  %-38s assoc %5d  genetic %4d  clinical %4d\n",
              substr(dt$query, 1, 37), dt$n_associated, dt$n_genetic, dt$n_clinical))
  dt
}), fill = TRUE)

setorder(feas, -n_clinical)
fwrite(feas, file.path(OUT, "disease_feasibility.csv"))

cat("\n=== feasibility ===\n")
print(feas)

MIN <- 30L
usable <- feas[n_clinical >= MIN]
cat("\n", nrow(usable), " diseases with >= ", MIN,
    " drug-evidence genes\n", sep = "")
cat("median associated targets: ", median(feas$n_associated), "\n", sep = "")

if (nrow(usable) < 10)
  cat("\nFewer than 10 usable diseases. Lower the threshold, or reconsider\n",
      "the design before building the rest of the pipeline.\n", sep = "")

msg("01_select_diseases.R complete -> ", file.path(OUT, "disease_feasibility.csv"))
