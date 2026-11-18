# network-target-benchmark

Code and results for *Network methods find druggable proteins, not disease genes: a controlled benchmark across 33 diseases*.

## What this is

Network-based methods are widely used to nominate drug targets, most often by ranking proteins on centrality measures such as degree or maximum neighborhood component. This benchmark asks two questions that are rarely tested:

1. **How much of a centrality ranking's apparent performance is just connectivity?** Every method is evaluated twice — conventionally, and against negatives matched on degree.
2. **How much of the remaining signal is actually disease-specific?** A model given no disease information at all, trained only to recognise approved-drug targets, is compared against disease-specific models.

Headline results: centrality ranking loses 88–100% of its above-chance performance under degree matching, while methods that learn from network structure lose 4–20%. A model with no disease information reaches AUROC 0.897 against individual diseases' labels, against 0.933 for disease-specific models — roughly 92% of the signal requires no disease information.

## Requirements

- R 4.3+ with `data.table`, `igraph`, `ggplot2`, `patchwork`, `httr`, `jsonlite`
- Python 3.12 with `torch`, `torch_geometric`, `xgboost`, `scikit-learn`, `scipy`, `numpy`, `pandas`
- A CUDA GPU is used for the graph attention network but is not required

Full benchmark runtime is approximately 100 minutes on an NVIDIA GB10.

## Data sources

Both are public and downloaded by the pipeline:

- **STRING v12.0** (taxon 9606), combined score ≥ 700
- **Open Targets Platform**, via the GraphQL API

No access credentials are needed.

## Pipeline

Run in order from the repository root.

| Step | Script | Produces |
|---|---|---|
| 1 | `R/01_select_diseases.R` | `data/disease_feasibility.csv` — resolves disease names to ontology identifiers and counts labels |
| 1b | `R/01b_verify_disease_ids.R` | flags suspect identifier resolutions for manual review |
| 1c | `R/01c_test_ms_ids.R` | probes candidate multiple sclerosis identifiers (see note below) |
| 2 | `R/02_build_assets.R` | `data/edges.csv`, `data/nodes_topology.csv`, `data/graph.rds`, `data/tractability.csv`, `data/disease_assoc.csv` |
| 3 | `python/benchmark.py` | `results/benchmark_raw.csv`, `results/benchmark_summary.csv`, `results/degree_effect.csv`, `results/final_predictions.npz` |
| 4 | `python/agnostic_baseline.py` | disease-agnostic control results |
| 5 | `python/cross_disease.py` | cross-disease correlation and prediction matrices |
| 6 | `R/03_figures.R` | Figures 2–5 as PDF, PNG, TIFF (600 dpi) and EPS |
| 7 | `R/04_supplementary.R` | Supplementary Tables S1–S5 |

`python/benchmark.py --ablation-only` reruns just the ablation without the full benchmark.

## Design notes

**Degree-matched evaluation.** For each positive, up to three negatives are drawn from the same degree stratum using 20 quantile bins over the positives' degree range. Connectivity therefore carries no information about class membership inside the labelled set. The full network is retained for message passing and diffusion; only the labelled node set is restricted.

**Excluded features.** `clinical` evidence is the label. `literature` and `genetic_literature` evidence measure how much has been published about a gene–disease pair and would import the publication bias this study sets out to measure. An assertion in `benchmark.py` halts execution if any of these reaches the feature matrix.

**Label caveat.** Positives derive from Open Targets clinical evidence, which propagates a drug's evidence across its mechanism. Metformin inhibits mitochondrial complex I, so all 41 complex I subunits in the network carry strong clinical evidence for type 2 diabetes. This inflates connectivity-correlated methods under conventional evaluation — which is the confound the matched control exists to measure. The definition is retained for comparability with prior benchmarks.

**Multiple sclerosis.** One of the 34 pre-specified diseases was dropped. Open Targets holds no generic multiple sclerosis entry carrying association data; the only term with substantial evidence is a symptom-measurement trait. Substituting a subtype after the pre-specified term failed would have been a post-hoc choice. 33 diseases are analysed. `R/01c_test_ms_ids.R` documents the identifiers tested.

## Reproducibility

Seeds are fixed throughout. Cross-validation is three repeats of five-fold stratified splits. Neural models use a stratified 20% validation split carved from each training fold for early stopping; validation data never enters training and test folds are never used for model selection.

## Citation

Manuscript under review. Citation details will be added on acceptance.

## License

MIT. See `LICENSE`.

Note that the data this code retrieves carry their own terms: STRING is CC BY 4.0, and Open Targets data are CC0 with per-source licences on some evidence streams. This licence covers the analysis code only.
