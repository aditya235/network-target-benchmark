"""
subnetwork_centrality.py -- does the degree confound survive inside a
disease-specific subnetwork?

The main benchmark computes centrality once on the full interactome. That
is not the applied workflow. The standard CytoHubba procedure builds a
subnetwork from disease-associated genes first, then ranks hubs *within
that subnetwork*. Degree in an induced subgraph is not degree in the full
graph, so the main result does not transfer automatically. This script
rebuilds that workflow per disease and reruns the same control.

Two arms, because practice splits:

  seeds       induced subgraph on the disease's own gene set
  expanded    that gene set plus its first neighbours in STRING, which is
              the "add interactors" step in STRING and Cytoscape

Methods are the CytoHubba measures that are affordable at this scale:
Degree, MNC, MCC, Bottleneck and Betweenness. MCC is the measure the
CytoHubba authors report as best performing, so it is the one that
matters most for the paper's claim.

Seed genes are chosen on non-label evidence only (genetic association,
animal model, somatic mutation). Clinical evidence is the label and must
never enter selection. Literature evidence is excluded for the same
reason as in the main benchmark: it encodes how much has been published
about the pair, which is the bias under test.

Outputs
  results/subnetwork_centrality.csv   per disease, arm, method
  results/subnetwork_summary.csv      pooled across diseases
  results/subnetwork_degree_corr.csv  subgraph degree vs global degree
  results/subnetwork_meta.json        parameters and everything skipped
"""

import argparse, json, time
from collections import deque
from math import lgamma
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx
from scipy.special import logsumexp
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

DATA, RES = Path("data"), Path("results")
RES.mkdir(exist_ok=True)

SELECT = ["genetic", "animal", "somatic"]       # never clinical, never literature
FORBIDDEN = {"clinical", "literature"}
MIN_POS = 10
SAMPLES = 300           # sampled roots for betweenness and bottleneck on big graphs
CLIQUE_BUDGET = 500_000  # MCC is abandoned past this; enumeration is super-polynomial


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


# ------------------------------------------------------------- centrality

def mnc(H):
    """Maximum neighborhood component: size of the largest connected
    component of the subgraph induced on a node's neighbours."""
    out = {}
    for v in H:
        nb = list(H[v])
        if len(nb) < 2:
            out[v] = len(nb); continue
        sub = H.subgraph(nb)
        out[v] = max((len(c) for c in nx.connected_components(sub)), default=0)
    return out


def mcc(H, budget=CLIQUE_BUDGET):
    """Maximal Clique Centrality: sum over maximal cliques C containing v
    of (|C|-1)!. Accumulated in log space -- (|C|-1)! overflows float64
    above |C| ~ 170 and STRING complexes get large. Log is monotonic, and
    only the ranking is used, so this changes nothing about the result.

    Returns None if enumeration exceeds the budget. Maximal-clique
    enumeration is super-polynomial; on a neighbour-expanded STRING
    subgraph it does not terminate in any usable time."""
    terms = {v: [] for v in H}
    n = 0
    for c in nx.find_cliques(H):
        n += 1
        if n > budget:
            return None
        w = lgamma(len(c))          # log((|C|-1)!)
        for v in c:
            terms[v].append(w)
    return {v: (logsumexp(t) if t else -np.inf) for v, t in terms.items()}


def bottleneck(H, roots):
    """CytoHubba's BottleNeck. For each root w a shortest-path tree is
    built; v scores a point for every w where more than a quarter of the
    paths out of w pass through v. Exact over `roots`; sampled roots give
    a proportional estimate, which is all the ranking needs."""
    bn = dict.fromkeys(H, 0)
    for w in roots:
        parent, order = {w: None}, [w]
        q = deque([w])
        while q:
            u = q.popleft()
            for nb in H[u]:
                if nb not in parent:
                    parent[nb] = u; order.append(nb); q.append(nb)
        size = dict.fromkeys(order, 1)
        thresh = len(order) / 4.0
        for u in reversed(order):            # children precede parents reversed
            p = parent[u]
            if p is not None:
                size[p] += size[u]
            if u != w and size[u] > thresh:
                bn[u] += 1
    return bn


def centralities(H, seed, allow_mcc):
    n = H.number_of_nodes()
    rng = np.random.default_rng(seed)
    nodes = list(H)
    roots = nodes if n <= SAMPLES else list(rng.choice(nodes, SAMPLES, replace=False))

    out = {"Degree": dict(H.degree()),
           "MNC": mnc(H),
           "Bottleneck": bottleneck(H, roots)}
    out["Betweenness"] = (nx.betweenness_centrality(H) if n <= SAMPLES
                          else nx.betweenness_centrality(H, k=SAMPLES, seed=seed))
    out["MCC"] = mcc(H) if allow_mcc else None
    return out


# --------------------------------------------------------------- control

def degree_matched(degree, y, ratio=3, n_bins=20, seed=0):
    """Positives plus up to `ratio` negatives from the same degree stratum.
    Inside the mask, degree carries no information about class. Identical
    scheme to the main benchmark."""
    rng = np.random.default_rng(seed)
    pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
    if len(pos) == 0 or len(neg) == 0:
        return np.zeros(len(y), bool)
    edges = np.unique(np.quantile(degree[pos], np.linspace(0, 1, n_bins + 1)))
    if len(edges) < 2:
        return np.zeros(len(y), bool)
    keep = list(pos)
    pb, nb = np.digitize(degree[pos], edges[1:-1]), np.digitize(degree[neg], edges[1:-1])
    for b in np.unique(pb):
        pool = neg[nb == b]
        if len(pool):
            want = min(int((pb == b).sum() * ratio), len(pool))
            keep += list(rng.choice(pool, want, replace=False))
    m = np.zeros(len(y), bool); m[np.asarray(keep, int)] = True
    return m


# ------------------------------------------------------------------ main

def run(n_seeds, ratio, seed, diseases=None):
    e = pd.read_csv(DATA / "edges.csv")
    G = nx.from_pandas_edgelist(e, "source", "target")
    assoc = pd.read_csv(DATA / "disease_assoc.csv")
    assert not FORBIDDEN & set(SELECT), "label evidence must not drive selection"
    gdeg, nodes = dict(G.degree()), set(G)

    rows, corr, skipped = [], [], []
    groups = sorted(assoc.groupby("disease"))
    if diseases:
        groups = [(d, g) for d, g in groups if d in diseases]

    for i, (dis, g) in enumerate(groups, 1):
        gg = g[g.gene.isin(nodes)].copy()
        gg["sc"] = gg[SELECT].max(axis=1)
        sd = set(gg[gg.sc > 0].nlargest(n_seeds, "sc").gene)
        if not sd:
            skipped.append((dis, "all", "no seed genes")); continue
        lab = dict(zip(g.gene, g.clinical))

        for arm, expand in (("seeds", False), ("expanded", True)):
            V = set(sd)
            if expand:
                for v in sd:
                    V |= set(G[v])
            H = G.subgraph(V).copy()
            if H.number_of_nodes() < 50:
                skipped.append((dis, arm, "subnetwork under 50 nodes")); continue
            V = list(H)
            y = np.array([1 if lab.get(v, 0) > 0 else 0 for v in V])
            if y.sum() < MIN_POS:
                skipped.append((dis, arm, f"{int(y.sum())} positives")); continue

            t0 = time.time()
            cen = centralities(H, seed, allow_mcc=not expand)
            sdeg = np.array([cen["Degree"][v] for v in V], float)
            g_deg = np.array([gdeg[v] for v in V], float)
            corr.append(dict(disease=dis, arm=arm, n_nodes=len(V),
                             n_edges=H.number_of_edges(), n_pos=int(y.sum()),
                             spearman_sub_vs_global=spearmanr(sdeg, g_deg).statistic))

            mask = degree_matched(sdeg, y, ratio, seed=seed)
            usable = mask.sum() > 0 and 0 < y[mask].sum() < mask.sum()

            scored = {k: v for k, v in cen.items() if v is not None}
            scored["Global degree"] = dict(zip(V, g_deg))
            if cen["MCC"] is None:
                skipped.append((dis, arm, "MCC: clique budget exceeded"))

            for m, d in scored.items():
                s = np.array([d[v] for v in V], float)
                s = np.where(np.isfinite(s), s, np.nanmin(s[np.isfinite(s)]) - 1
                             if np.isfinite(s).any() else 0.0)
                rows.append(dict(
                    disease=dis, arm=arm, method=m, n_nodes=len(V),
                    n_edges=H.number_of_edges(), n_pos=int(y.sum()),
                    AUROC_subnetwork=roc_auc_score(y, s),
                    AUROC_matched=roc_auc_score(y[mask], s[mask]) if usable else np.nan,
                    n_matched=int(mask.sum()) if usable else 0))
            log(f"{i:>2} {dis:<42} {arm:<9} {len(V):>6} nodes  {time.time()-t0:5.1f}s")

    raw = pd.DataFrame(rows)
    raw.to_csv(RES / "subnetwork_centrality.csv", index=False)
    pd.DataFrame(corr).to_csv(RES / "subnetwork_degree_corr.csv", index=False)

    summ = (raw.groupby(["arm", "method"])
               .agg(n_diseases=("disease", "nunique"),
                    AUROC_subnetwork=("AUROC_subnetwork", "mean"),
                    AUROC_subnetwork_sd=("AUROC_subnetwork", "std"),
                    AUROC_matched=("AUROC_matched", "mean"),
                    AUROC_matched_sd=("AUROC_matched", "std")).reset_index())
    summ["degree_share_pct"] = (100 * (summ.AUROC_subnetwork - summ.AUROC_matched)
                                / (summ.AUROC_subnetwork - 0.5)).round(1)
    summ = summ.round(4)
    summ.to_csv(RES / "subnetwork_summary.csv", index=False)

    (RES / "subnetwork_meta.json").write_text(json.dumps(dict(
        n_seeds=n_seeds, ratio=ratio, seed=seed, select_evidence=SELECT,
        min_positives=MIN_POS, samples=SAMPLES, clique_budget=CLIQUE_BUDGET,
        skipped=skipped), indent=2))

    print(); print(summ.to_string(index=False))
    if skipped:
        print(f"\nnot evaluable ({len(skipped)}):")
        for d, a, why in skipped[:12]:
            print(f"   {d} [{a}] -- {why}")
    return summ


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--n-seeds", type=int, default=500)
    p.add_argument("--ratio", type=int, default=3)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--diseases", nargs="*", default=None,
                   help="restrict to these diseases (smoke test)")
    a = p.parse_args()
    log(f"subnetwork centrality: {a.n_seeds} seeds, {a.ratio}:1 matching")
    run(a.n_seeds, a.ratio, a.seed, a.diseases)
    log("done")
