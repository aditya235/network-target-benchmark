# ============================================================
# 01c_test_ms_ids.R -- which identifier does Open Targets accept for
# multiple sclerosis?
#
# MONDO_0005301 is the correct MONDO term for multiple sclerosis, but the
# Open Targets API returned nothing for it and the disease was silently
# dropped from the benchmark. Rather than guess again, this queries each
# candidate directly and reports which returns data.
#
# Generalised: give it any disease name and it will test every candidate
# the search API offers, plus any IDs supplied by hand.
# ============================================================

suppressPackageStartupMessages({ library(httr2); library(jsonlite); library(data.table) })

if (!exists("%||%")) `%||%` <- function(x, y) if (is.null(x)) y else x
OT <- "https://api.platform.opentargets.org/api/v4/graphql"
gql <- function(q, v = NULL) {
  r <- tryCatch(request(OT) |> req_body_json(list(query = q, variables = v)) |>
                  req_timeout(45) |> req_error(is_error = function(x) FALSE) |>
                  req_perform(), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  tryCatch(resp_body_json(r, simplifyVector = FALSE), error = function(e) NULL)
}

# does this identifier resolve, and how much evidence does it carry?
probe <- function(id) {
  d <- gql('query($id:String!){ disease(efoId:$id){ id name
             associatedTargets(page:{index:0,size:1}){ count } } }',
           list(id = id))
  dis <- d$data$disease
  if (is.null(dis)) {
    st <- "no data"
    if (!is.null(d$errors)) st <- substr(d$errors[[1]]$message, 1, 60)
    return(data.table(id = id, name = NA_character_,
                      n_targets = NA_integer_, status = st))
  }
  data.table(id = id, name = dis$name %||% NA_character_,
             n_targets = dis$associatedTargets$count %||% NA_integer_,
             status = "OK")
}
DISEASE <- "multiple sclerosis"

# candidates from the search API
hits <- gql('query($q:String!){ search(queryString:$q, entityNames:["disease"],
              page:{index:0,size:15}){ hits { id name } } }',
            list(q = DISEASE))$data$search$hits
from_search <- data.table(id = character(), search_name = character())
if (length(hits)) {
  from_search <- data.table(
    id = vapply(hits, function(h) h$id, character(1)),
    search_name = vapply(hits, function(h) h$name %||% "", character(1)))
}

# plus identifiers worth testing explicitly: the canonical MONDO term and
# the EFO term, neither of which the search returned at the top
manual <- data.table(
  id = c("MONDO_0005301", "EFO_0003885", "MONDO_0005314", "MONDO_0005284"),
  search_name = c("MONDO: multiple sclerosis (canonical)",
                  "EFO: multiple sclerosis (if it exists)",
                  "relapsing-remitting MS", "chronic progressive MS"))

cand <- unique(rbind(manual, from_search), by = "id")
cat("testing", nrow(cand), "candidate identifiers\n\n")

res <- rbindlist(lapply(seq_len(nrow(cand)), function(i) {
  r <- probe(cand$id[i]); Sys.sleep(0.25)
  r[, search_name := cand$search_name[i]][]
}), fill = TRUE)

setorder(res, -n_targets, na.last = TRUE)
print(res[, .(id, name, n_targets, status, search_name)])

ok <- res[status == "OK" & !is.na(n_targets) & n_targets > 1000]
cat("\n")
if (nrow(ok)) {
  cat("USABLE (>1000 associated targets):\n")
  print(ok[, .(id, name, n_targets)])
  cat("\nPick the one whose name is the disease itself, not a subtype or\n")
  cat("a measurement, and put it in data/disease_feasibility.csv.\n")
} else {
  cat("No candidate returned substantial data. Open Targets may not hold a\n")
  cat("generic multiple sclerosis entry; consider using the most common\n")
  cat("subtype (relapsing-remitting) and saying so in the methods, or\n")
  cat("dropping the disease and reporting 33.\n")
}
