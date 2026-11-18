# ============================================================
# 01b_verify_disease_ids.R -- check what each disease name resolved to.
#
# WHY: automated ontology resolution fails silently. "multiple sclerosis"
# resolved to EFO_0803536, "multiple sclerosis symptom measurement" -- a
# GWAS measurement trait, not the disease. The entity filter did not catch
# it because Open Targets classifies such traits as diseases too.
#
# This prints the top candidates per name, flags the ones that look wrong,
# and writes a table to confirm by eye. Thirty-four lookups is worth an
# hour; a wrong identifier is invisible in every downstream result.
#
# Outputs: data/disease_id_candidates.csv
# ============================================================

suppressPackageStartupMessages({
  library(data.table); library(httr2); library(jsonlite)
})

OT <- "https://api.platform.opentargets.org/api/v4/graphql"
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
if (!exists("%||%")) `%||%` <- function(x, y) if (is.null(x)) y else x

gql <- function(q, v = NULL) {
  r <- tryCatch(request(OT) |> req_body_json(list(query = q, variables = v)) |>
                  req_timeout(45) |> req_error(is_error = function(x) FALSE) |>
                  req_perform(), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  tryCatch(resp_body_json(r, simplifyVector = FALSE)$data, error = function(e) NULL)
}

Q <- 'query($q:String!){ search(queryString:$q, entityNames:["disease"],
        page:{index:0,size:6}){ hits { id name description } } }'

feas <- fread(file.path("data", "disease_feasibility.csv"))
names_to_check <- feas$query

# Terms that mark a hit as a measurement, phenotype or procedure rather
# than a disease. These are the silent failures.
SUSPECT <- paste(c("measurement", "symptom", "score", "scale", "level",
                   "biomarker", "procedure", "medication", "self-reported",
                   "family history", "risk", "treatment", "response"),
                 collapse = "|")

res <- rbindlist(lapply(names_to_check, function(nm) {
  d <- gql(Q, list(q = nm)); Sys.sleep(0.2)
  hits <- d$search$hits
  if (is.null(hits) || !length(hits))
    return(data.table(query = nm, rank = NA_integer_, id = NA_character_,
                      name = NA_character_, suspect = NA, is_mondo = NA))
  rbindlist(lapply(seq_along(hits), function(i) {
    h <- hits[[i]]
    data.table(query = nm, rank = i, id = h$id,
               name = h$name %||% NA_character_,
               suspect = grepl(SUSPECT, h$name %||% "", ignore.case = TRUE),
               is_mondo = grepl("^MONDO", h$id))
  }))
}), fill = TRUE)

fwrite(res, file.path("data", "disease_id_candidates.csv"))

# ---- what was actually used, and does it look right? --------------------
used <- merge(feas[, .(query, used_id = efo_id, used_name = matched_name)],
              res[rank == 1, .(query, top_id = id, top_name = name)],
              by = "query", all.x = TRUE)
used[, suspect := grepl(SUSPECT, used_name, ignore.case = TRUE)]
used[, is_mondo := grepl("^MONDO", used_id)]
# does the resolved name resemble the query?
used[, name_match := mapply(function(a, b) {
  a <- tolower(gsub("[^a-z ]", "", tolower(a)))
  b <- tolower(gsub("[^a-z ]", "", tolower(b)))
  grepl(substr(a, 1, min(nchar(a), 12)), b, fixed = TRUE) ||
    grepl(substr(b, 1, min(nchar(b), 12)), a, fixed = TRUE)
}, query, used_name)]

cat("\n================ FLAGGED ================\n")
bad <- used[suspect == TRUE | is_mondo == FALSE | name_match == FALSE]
if (nrow(bad)) {
  print(bad[, .(query, used_id, used_name, suspect, is_mondo, name_match)])
  cat("\nAlternatives for the flagged entries:\n")
  for (q in bad$query) {
    cat("\n", q, "\n", sep = "")
    print(res[query == q, .(rank, id, name, suspect)])
  }
} else {
  cat("none\n")
}

cat("\n================ SUMMARY ================\n")
cat("  total diseases        : ", nrow(used), "\n", sep = "")
cat("  MONDO identifiers     : ", sum(used$is_mondo, na.rm = TRUE), "\n", sep = "")
cat("  flagged for review    : ", nrow(bad), "\n", sep = "")
cat("\nReview every flagged row by hand, correct the identifier in\n")
cat("data/disease_feasibility.csv, then rerun 02 and the benchmark for\n")
cat("the affected diseases.\n")

msg("01b_verify_disease_ids.R complete -> data/disease_id_candidates.csv")
