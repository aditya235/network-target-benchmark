"""
agnostic_baseline.py -- does disease information add anything at all?

The decisive control. Train one model with NO disease information:
labels are the 946 proteins that are targets of an approved drug for any
indication; features are network topology only. Then evaluate it on each
disease's own labels.

Three numbers to compare, per disease:
    own       -- that disease's own cross-validated model
    other     -- another disease's model, median across the other 33
    agnostic  -- this model, which has never seen any disease data

If `agnostic` approaches `own`, disease-specific target prediction is
mostly recovering the druggable proteome, and the disease-specific
machinery is doing little work.

Evaluation is degree-matched throughout, since the uncontrolled numbers
are inflated by connectivity.

Usage: python python/agnostic_baseline.py
"""

import json
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

from benchmark import (load_shared, spectral_embed, xgb_scores,
                       degree_matched_subset, TOPO, log)

DATA, RESULTS = "data", "results"


def main(folds=5, repeats=3, seed=2026):
    nodes, genes, idx, edge_index, edge_weight = load_shared()
    n = len(genes)

    # ---- disease-agnostic labels -------------------------------------
    tract = pd.read_csv(f"{DATA}/tractability.csv")
    approved = set(tract.loc[(tract.label_tr == "Approved Drug") &
                             (tract.value == True), "gene"])
    y_ag = np.array([1 if g in approved else 0 for g in genes])
    log(f"disease-agnostic positives: {y_ag.sum()}")

    log("spectral embedding ...")
    emb = spectral_embed(edge_index, edge_weight, n)
    X = StandardScaler().fit_transform(nodes[TOPO].to_numpy(dtype=float))
    Xe = np.hstack([X, emb])
    degree = nodes["degree"].to_numpy()

    # ---- train the agnostic model, cross-validated -------------------
    # Averaged over folds so no single split drives the result, and the
    # score for every protein comes from a model that did not train on it.
    sub = degree_matched_subset(degree, y_ag, seed=seed)
    ids = np.where(sub)[0]

    # Strictly out-of-fold. Nodes inside the training subset get the
    # average of only those folds that held them out; nodes outside it
    # were never trained on at all, so every fold is out-of-sample for
    # them and all are averaged. Accumulating the full prediction vector
    # each fold -- which an earlier version did -- makes the score
    # in-sample and yields a meaningless AUROC of 1.000.
    oof = np.zeros(n); cnt = np.zeros(n)
    outside = ~sub
    for rep in range(repeats):
        skf = StratifiedKFold(folds, shuffle=True, random_state=seed + rep)
        for tr_i, te_i in skf.split(ids, y_ag[ids]):
            tm = np.zeros(n, bool); tm[ids[tr_i]] = True
            p = xgb_scores(Xe, y_ag, tm, seed + rep)
            held = np.zeros(n, bool); held[ids[te_i]] = True
            held |= outside                       # never seen in training
            oof[held] += p[held]; cnt[held] += 1
    if (cnt == 0).any():
        log(f"warning: {int((cnt==0).sum())} nodes never scored")
        cnt[cnt == 0] = 1
    agnostic_score = oof / cnt
    log(f"agnostic model AUROC on its own labels, out-of-fold "
        f"(degree-matched): {roc_auc_score(y_ag[sub], agnostic_score[sub]):.4f}")

    # ---- compare, per disease ----------------------------------------
    z = np.load(f"{RESULTS}/final_predictions.npz", allow_pickle=True)
    diseases = list(z["diseases"]); labels = z["labels"]
    cross = pd.read_csv(f"{RESULTS}/cross_prediction_auroc_spectral.csv", index_col=0)
    bench = pd.read_csv(f"{RESULTS}/benchmark_raw.csv")
    own_cv = (bench[(bench.evaluation == "degree_matched") &
                    (bench.method == "spectral + XGBoost")]
              .groupby("disease")["AUROC"].mean())

    rows = []
    for j, dis in enumerate(diseases):
        y = labels[j]
        if y.sum() < 20:
            continue
        # degree-matched evaluation set for THIS disease
        m = degree_matched_subset(degree, y, seed=seed)
        if m is None:
            continue
        ag = roc_auc_score(y[m], agnostic_score[m])
        others = np.delete(cross.iloc[:, j].to_numpy(), j)
        # How disease-specific are the labels themselves? If a disease's
        # "drug targets" are mostly just general approved-drug targets,
        # then a disease-agnostic model predicting them well is no
        # surprise -- and the real finding is about the labels, not the
        # models.
        pos = np.where(y == 1)[0]
        frac_in_approved = float(y_ag[pos].mean()) if len(pos) else np.nan
        rows.append({
            "disease": dis,
            "n_labels": int(y.sum()),
            "pct_labels_also_approved": round(100 * frac_in_approved, 1),
            "own_cv_auroc": round(float(own_cv.get(dis, np.nan)), 4),
            "other_disease_auroc": round(float(np.nanmedian(others)), 4),
            "agnostic_auroc": round(float(ag), 4),
        })

    df = pd.DataFrame(rows)
    df["disease_advantage"] = (df.own_cv_auroc - df.agnostic_auroc).round(4)
    df = df.sort_values("disease_advantage")
    df.to_csv(f"{RESULTS}/agnostic_comparison.csv", index=False)

    print("\n", df.to_string(index=False), flush=True)
    print("\n=== summary (degree-matched AUROC, median across diseases) ===")
    print(f"  own disease model      : {df.own_cv_auroc.median():.4f}")
    print(f"  another disease's model: {df.other_disease_auroc.median():.4f}")
    print(f"  no disease info at all : {df.agnostic_auroc.median():.4f}")
    adv = df.disease_advantage.median()
    print(f"\n  advantage of knowing the disease: {adv:+.4f} AUROC")
    print(f"  diseases where it helps by >0.02 : "
          f"{int((df.disease_advantage > 0.02).sum())}/{len(df)}")
    print(f"\n  median % of a disease's drug-target labels that are also "
          f"general approved-drug targets: {df.pct_labels_also_approved.median():.1f}%")
    print("\nIf the advantage is near zero, disease-specific target")
    print("prediction largely recovers the druggable proteome.", flush=True)

    json.dump({"own_median": float(df.own_cv_auroc.median()),
               "other_median": float(df.other_disease_auroc.median()),
               "agnostic_median": float(df.agnostic_auroc.median()),
               "advantage_median": float(adv),
               "n_diseases": len(df)},
              open(f"{RESULTS}/agnostic_summary.json", "w"), indent=2)
    log("agnostic_baseline.py complete")


if __name__ == "__main__":
    main()
