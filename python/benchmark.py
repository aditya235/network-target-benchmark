"""
benchmark.py -- target prediction across 34 diseases, with degree control.

Design notes that matter:

* The interactome, its topology features and the spectral embedding are
  disease-independent, so they are computed ONCE and reused. Only the
  association scores and labels change per disease. This is what makes
  34 diseases affordable.

* Every disease is evaluated twice: on all nodes (comparable to the
  published benchmarks) and on a degree-matched subset. The gap between
  them is how much of the reported performance was connectivity.

* Per-disease prediction vectors are saved so that rankings can be
  correlated BETWEEN diseases afterwards. If disease-specific models
  produce near-identical rankings, "disease-specific prediction" is a
  relabelling exercise.

* `clinical` is the label and must never be a feature; `literature` and
  `genetic_literature` encode prior publication and are excluded too.
  A guard enforces both.

Usage:
    python python/benchmark.py                 # all diseases
    python python/benchmark.py --quick         # 3 diseases, 1 repeat
    python python/benchmark.py --max-diseases 10
"""

import argparse, json, os, time, warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
import torch.nn.functional as F
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.preprocessing import StandardScaler
from torch_geometric.data import Data
from torch_geometric.nn import GATConv, SAGEConv

warnings.filterwarnings("ignore")
DATA, RESULTS = "data", "results"
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

TOPO = ["degree", "strength", "betweenness", "closeness", "eigen",
        "pagerank", "kcore", "clustering", "mnc"]
EVID = ["genetic", "animal", "somatic"]          # NOT clinical, NOT literature
FORBIDDEN = {"clinical", "literature", "label"}


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


# ----------------------------------------------------------------------
# shared assets
# ----------------------------------------------------------------------
def load_shared():
    nodes = pd.read_csv(f"{DATA}/nodes_topology.csv")
    edges = pd.read_csv(f"{DATA}/edges.csv")
    genes = nodes["gene"].tolist()
    idx = {g: i for i, g in enumerate(genes)}

    e = edges[edges.source.isin(idx) & edges.target.isin(idx)]
    src, dst = e.source.map(idx).to_numpy(), e.target.map(idx).to_numpy()
    edge_index = torch.tensor(np.vstack([np.concatenate([src, dst]),
                                         np.concatenate([dst, src])]), dtype=torch.long)
    edge_weight = torch.tensor(np.concatenate([e.score.to_numpy()] * 2) / 1000.0,
                               dtype=torch.float)
    log(f"graph: {len(genes)} nodes, {len(e)} edges")
    return nodes, genes, idx, edge_index, edge_weight


def spectral_embed(edge_index, edge_weight, n, dim=64):
    """Unsupervised topology embedding. Disease-independent -> computed once."""
    from scipy.sparse.linalg import eigsh
    ei = edge_index.numpy()
    A = sp.csr_matrix((edge_weight.numpy(), (ei[0], ei[1])), shape=(n, n))
    d = np.asarray(A.sum(1)).ravel(); d[d == 0] = 1
    D = sp.diags(1 / np.sqrt(d))
    L = sp.identity(n) - D @ A @ D
    vals, vecs = eigsh(L, k=min(dim + 1, n - 2), sigma=0, which="LM")
    return vecs[:, 1:dim + 1]


# ----------------------------------------------------------------------
# models
# ----------------------------------------------------------------------
class GNN(torch.nn.Module):
    def __init__(self, kind, in_dim, hidden=64, heads=4, dropout=0.3):
        super().__init__()
        self.kind, self.dropout = kind, dropout
        if kind == "GAT":
            self.c1 = GATConv(in_dim, hidden // heads, heads=heads, dropout=dropout)
            self.c2 = GATConv(hidden, hidden // 2, heads=1, dropout=dropout)
        else:
            self.c1, self.c2 = SAGEConv(in_dim, hidden), SAGEConv(hidden, hidden // 2)
        self.out = torch.nn.Linear(hidden // 2, 2)

    def forward(self, x, ei):
        x = F.dropout(F.elu(self.c1(x, ei)), p=self.dropout, training=self.training)
        return self.out(F.elu(self.c2(x, ei)))


def train_gnn(x, ei, y, train_mask, kind="GAT", epochs=200, lr=0.05, seed=0):
    """Early stopping on a held-out validation split, never on training loss."""
    torch.manual_seed(seed)
    tr_i = np.where(train_mask)[0]
    if y[tr_i].sum() >= 10 and (y[tr_i] == 0).sum() >= 10:
        tr, va = train_test_split(tr_i, test_size=0.2, stratify=y[tr_i], random_state=seed)
    else:
        tr = va = tr_i
    tm = np.zeros(len(y), bool); tm[tr] = True
    vm = np.zeros(len(y), bool); vm[va] = True

    model = GNN(kind, x.shape[1]).to(DEVICE)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=5e-4)
    xt, eit = x.to(DEVICE), ei.to(DEVICE)
    yt = torch.tensor(y, dtype=torch.long, device=DEVICE)
    npos = int(y[tm].sum())
    w = torch.tensor([1.0, max((tm.sum() - npos) / max(npos, 1), 1.0)],
                     dtype=torch.float, device=DEVICE)
    tmt = torch.tensor(tm, device=DEVICE)

    best, best_state, patience = -1.0, None, 0
    for ep in range(epochs):
        model.train(); opt.zero_grad()
        F.cross_entropy(model(xt, eit)[tmt], yt[tmt], weight=w).backward(); opt.step()
        if ep % 5 == 4:
            model.eval()
            with torch.no_grad():
                p = F.softmax(model(xt, eit), 1)[:, 1].cpu().numpy()
            if y[vm].sum() < 2:
                continue
            ap = average_precision_score(y[vm], p[vm])
            if ap > best:
                best, patience = ap, 0
                best_state = {k: v.clone() for k, v in model.state_dict().items()}
            else:
                patience += 1
                if patience >= 8:
                    break
    if best_state:
        model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        return F.softmax(model(xt, eit), 1)[:, 1].cpu().numpy()


def xgb_scores(X, y, train_mask, seed=0):
    from xgboost import XGBClassifier
    npos = int(y[train_mask].sum())
    clf = XGBClassifier(n_estimators=300, max_depth=5, learning_rate=0.1,
                        subsample=0.8, colsample_bytree=0.8,
                        scale_pos_weight=max((train_mask.sum() - npos) / max(npos, 1), 1),
                        eval_metric="aucpr", random_state=seed, n_jobs=8,
                        tree_method="hist")
    clf.fit(X[train_mask], y[train_mask])
    return clf.predict_proba(X)[:, 1]


def rwr(A_norm, seed_mask, restart=0.3, iters=100):
    """Random walk with restart -- the family of methods the prior
    benchmarks tested, included here so the replication is like-for-like."""
    p0 = seed_mask.astype(float); s = p0.sum()
    if s == 0:
        return np.zeros(len(seed_mask))
    p0 /= s; p = p0.copy()
    for _ in range(iters):
        pn = (1 - restart) * (A_norm.T @ p) + restart * p0
        if np.abs(pn - p).sum() < 1e-10:
            return pn
        p = pn
    return p


# ----------------------------------------------------------------------
# evaluation
# ----------------------------------------------------------------------
def degree_matched_subset(degree, y, ratio=3, n_bins=20, seed=0):
    rng = np.random.default_rng(seed)
    pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
    if len(pos) < 10:
        return None
    edges_ = np.unique(np.quantile(degree[pos], np.linspace(0, 1, n_bins + 1)))
    pbin = np.digitize(degree[pos], edges_[1:-1])
    nbin = np.digitize(degree[neg], edges_[1:-1])
    keep = [pos]
    for b in np.unique(pbin):
        pool = neg[nbin == b]
        if len(pool) == 0:
            continue
        want = int(ratio * np.sum(pbin == b))
        keep.append(rng.choice(pool, size=min(want, len(pool)), replace=False))
    m = np.zeros(len(y), bool); m[np.concatenate(keep)] = True
    return m


def evaluate(y_true, score, k=50):
    order = np.argsort(-score); k = min(k, len(order))
    return {"AUROC": roc_auc_score(y_true, score),
            "AUPRC": average_precision_score(y_true, score),
            f"Precision@{k}": float(y_true[order[:k]].mean())}


def run_disease(disease, nodes, genes, idx, edge_index, edge_weight, emb,
                A_norm, assoc, folds=5, repeats=3, seed=2026):
    """One disease: all methods, both evaluations. Returns (rows, preds)."""
    n = len(genes)
    sub = assoc[assoc.disease == disease]
    if not len(sub):
        return [], None

    # align disease evidence onto the shared node order
    m = sub.set_index("gene").reindex(genes)
    y = (m["clinical"].fillna(0) > 0).astype(int).to_numpy()
    if y.sum() < 20:
        log(f"  {disease}: only {y.sum()} positives, skipped")
        return [], None

    feats = nodes[TOPO].copy()
    for c in EVID:
        feats[c] = m[c].fillna(0).to_numpy()
    bad = FORBIDDEN & set(feats.columns)
    if bad:
        raise RuntimeError(f"LEAKAGE: {bad} present as features")

    X = StandardScaler().fit_transform(feats.to_numpy(dtype=float))
    Xe = np.hstack([X, emb])
    xt = torch.tensor(X, dtype=torch.float)
    degree = nodes["degree"].to_numpy()

    dm = degree_matched_subset(degree, y, seed=seed)
    evals = {"all_nodes": np.ones(n, bool)}
    if dm is not None:
        evals["degree_matched"] = dm

    rows, store = [], {}
    for tag, subset in evals.items():
        ids = np.where(subset)[0]
        for rep in range(repeats):
            skf = StratifiedKFold(folds, shuffle=True, random_state=seed + rep)
            for fold, (tr_i, te_i) in enumerate(skf.split(ids, y[ids])):
                tm = np.zeros(n, bool); tm[ids[tr_i]] = True
                te = np.zeros(n, bool); te[ids[te_i]] = True
                if y[te].sum() < 2:
                    continue
                preds = {
                    "GAT": train_gnn(xt, edge_index, y, tm, seed=seed + rep),
                    "spectral + XGBoost": xgb_scores(Xe, y, tm, seed + rep),
                    "XGBoost (no graph)": xgb_scores(X, y, tm, seed + rep),
                    "RWR": rwr(A_norm, y * tm),
                    "Degree": degree.astype(float),
                    "MNC (CytoHubba)": nodes["mnc"].to_numpy(float),
                    "Betweenness": nodes["betweenness"].to_numpy(float),
                }
                for name, p in preds.items():
                    rows.append({"disease": disease, "evaluation": tag,
                                 "method": name, "repeat": rep, "fold": fold,
                                 "n_pos": int(y.sum()),
                                 **evaluate(y[te], p[te])})
                    if tag == "degree_matched" and rep == 0 and fold == 0:
                        store[name] = p
    # one full-data model per disease for the cross-disease comparison
    allm = np.ones(n, bool)
    store["_final_GAT"] = train_gnn(xt, edge_index, y, allm, seed=seed)
    store["_final_spectral"] = xgb_scores(Xe, y, allm, seed)
    store["_y"] = y
    return rows, store


# ----------------------------------------------------------------------
def main(folds=5, repeats=3, seed=2026, max_diseases=None, quick=False):
    os.makedirs(RESULTS, exist_ok=True)
    log(f"device: {DEVICE}")
    nodes, genes, idx, edge_index, edge_weight = load_shared()
    n = len(genes)

    log("computing shared spectral embedding ...")
    emb = spectral_embed(edge_index, edge_weight, n)

    ei = edge_index.numpy()
    A = sp.csr_matrix((edge_weight.numpy(), (ei[0], ei[1])), shape=(n, n))
    d = np.asarray(A.sum(1)).ravel(); d[d == 0] = 1
    A_norm = sp.diags(1 / d) @ A

    assoc = pd.read_csv(f"{DATA}/disease_assoc.csv")
    diseases = list(pd.unique(assoc["disease"]))
    if quick:
        diseases, repeats = diseases[:3], 1
    if max_diseases:
        diseases = diseases[:max_diseases]
    log(f"{len(diseases)} diseases, {repeats}x{folds} folds each")

    all_rows, final_preds = [], {}
    for i, dis in enumerate(diseases, 1):
        t0 = time.time()
        rows, store = run_disease(dis, nodes, genes, idx, edge_index, edge_weight,
                                  emb, A_norm, assoc, folds, repeats, seed)
        all_rows.extend(rows)
        if store is not None:
            final_preds[dis] = {"GAT": store["_final_GAT"],
                                "spectral": store["_final_spectral"],
                                "y": store["_y"]}
        log(f"  [{i}/{len(diseases)}] {dis}  ({time.time()-t0:.0f}s)")
        pd.DataFrame(all_rows).to_csv(f"{RESULTS}/benchmark_raw.csv", index=False)

    df = pd.DataFrame(all_rows)
    df.to_csv(f"{RESULTS}/benchmark_raw.csv", index=False)

    # per-disease prediction matrices, for the cross-disease analysis
    np.savez_compressed(
        f"{RESULTS}/final_predictions.npz",
        genes=np.array(genes),
        diseases=np.array(list(final_preds.keys())),
        gat=np.vstack([v["GAT"] for v in final_preds.values()]),
        spectral=np.vstack([v["spectral"] for v in final_preds.values()]),
        labels=np.vstack([v["y"] for v in final_preds.values()]))

    summ = (df.groupby(["evaluation", "method"])[["AUROC", "AUPRC"]]
              .agg(["mean", "std"]).round(4))
    summ.columns = ["_".join(c) for c in summ.columns]
    summ = summ.reset_index()
    summ.to_csv(f"{RESULTS}/benchmark_summary.csv", index=False)
    print("\n", summ.to_string(index=False), flush=True)

    if {"all_nodes", "degree_matched"} <= set(df.evaluation.unique()):
        piv = df.groupby(["method", "evaluation"])["AUROC"].mean().unstack()
        piv["attributable_to_degree"] = piv["all_nodes"] - piv["degree_matched"]
        piv = piv.sort_values("degree_matched", ascending=False).round(4)
        piv.to_csv(f"{RESULTS}/degree_effect.csv")
        print("\n--- how much was degree? ---")
        print(piv.to_string(), flush=True)

    json.dump({"n_diseases": len(final_preds), "n_nodes": n,
               "folds": folds, "repeats": repeats, "device": str(DEVICE)},
              open(f"{RESULTS}/run_metadata.json", "w"), indent=2)
    log("benchmark complete")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--max-diseases", type=int, default=None)
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    main(a.folds, a.repeats, a.seed, a.max_diseases, a.quick)
    log(f"total {(time.time()-t0)/60:.1f} min")
