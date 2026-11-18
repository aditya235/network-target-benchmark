"""
cross_disease.py -- are "disease-specific" predictions actually specific?

Three analyses, in increasing order of how damaging they could be:

1. Rank correlation between every pair of diseases' prediction vectors.
   If hypertension and Alzheimer's rankings correlate at rho = 0.95,
   the disease-specific part is largely decorative.

2. The same for the LABEL sets. This is the crucial control: if the
   labels themselves overlap heavily, correlated predictions are
   expected and unremarkable. The finding is only meaningful if
   predictions correlate MORE than labels do.

3. Cross-prediction: score disease j's labels using disease i's model.
   If disease A's model identifies disease B's targets about as well as
   B's own model does, the models are interchangeable.

Note on the diagonal: each final model was fit on all of its own
disease's labels, so within-disease values are optimistic. The honest
comparison is off-diagonal performance against the cross-validated
numbers in benchmark_summary.csv, not against the diagonal.

Usage: python python/cross_disease.py
"""

import json, os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

RESULTS = "results"


def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / max(len(a | b), 1)


def main(top_k=100):
    z = np.load(f"{RESULTS}/final_predictions.npz", allow_pickle=True)
    genes = z["genes"]; diseases = list(z["diseases"])
    labels = z["labels"]                       # (D, N) 0/1
    D, N = labels.shape
    print(f"{D} diseases, {N} proteins\n", flush=True)

    out = {}
    for model in ("gat", "spectral"):
        P = z[model]                            # (D, N) scores

        # ---- 1. prediction correlation --------------------------------
        rho = np.ones((D, D))
        for i in range(D):
            for j in range(i + 1, D):
                r = spearmanr(P[i], P[j]).statistic
                rho[i, j] = rho[j, i] = r
        off = rho[np.triu_indices(D, 1)]

        # ---- 2. label overlap, the control ----------------------------
        lab_j = np.ones((D, D))
        top_j = np.ones((D, D))
        for i in range(D):
            ti = np.argsort(-P[i])[:top_k]
            for j in range(i + 1, D):
                lab_j[i, j] = lab_j[j, i] = jaccard(np.where(labels[i])[0],
                                                    np.where(labels[j])[0])
                tj = np.argsort(-P[j])[:top_k]
                top_j[i, j] = top_j[j, i] = jaccard(ti, tj)
        lab_off = lab_j[np.triu_indices(D, 1)]
        top_off = top_j[np.triu_indices(D, 1)]

        # ---- 3. cross-prediction AUROC --------------------------------
        X = np.full((D, D), np.nan)
        for i in range(D):
            for j in range(D):
                if labels[j].sum() < 5:
                    continue
                X[i, j] = roc_auc_score(labels[j], P[i])
        diag = np.diag(X)
        offd = X[~np.eye(D, dtype=bool)]

        print(f"=== {model.upper()} ===")
        print(f"  prediction rank correlation : median {np.median(off):.3f} "
              f"(IQR {np.percentile(off,25):.3f}-{np.percentile(off,75):.3f})")
        print(f"  label Jaccard overlap       : median {np.median(lab_off):.3f}")
        print(f"  top-{top_k} prediction overlap   : median {np.median(top_off):.3f}")
        print(f"  cross-prediction AUROC      : own {np.nanmedian(diag):.3f}  "
              f"other {np.nanmedian(offd):.3f}  "
              f"ratio {np.nanmedian(offd)/np.nanmedian(diag):.3f}")
        print(flush=True)

        out[model] = {
            "pred_rho_median": float(np.median(off)),
            "pred_rho_iqr": [float(np.percentile(off, 25)), float(np.percentile(off, 75))],
            "label_jaccard_median": float(np.median(lab_off)),
            f"top{top_k}_overlap_median": float(np.median(top_off)),
            "cross_auroc_own": float(np.nanmedian(diag)),
            "cross_auroc_other": float(np.nanmedian(offd)),
            "cross_auroc_ratio": float(np.nanmedian(offd) / np.nanmedian(diag)),
        }

        pd.DataFrame(rho, index=diseases, columns=diseases).round(4) \
            .to_csv(f"{RESULTS}/pred_correlation_{model}.csv")
        pd.DataFrame(X, index=diseases, columns=diseases).round(4) \
            .to_csv(f"{RESULTS}/cross_prediction_auroc_{model}.csv")
        pd.DataFrame(lab_j, index=diseases, columns=diseases).round(4) \
            .to_csv(f"{RESULTS}/label_jaccard.csv")

        # per-disease: how much better is its own model than the average
        # of all the others at finding its targets?
        rows = []
        for j, dis in enumerate(diseases):
            others = np.delete(X[:, j], j)
            rows.append({"disease": dis, "n_labels": int(labels[j].sum()),
                         "own_model_auroc": round(X[j, j], 4),
                         "other_models_median": round(float(np.nanmedian(others)), 4),
                         "advantage": round(float(X[j, j] - np.nanmedian(others)), 4)})
        pd.DataFrame(rows).sort_values("advantage") \
            .to_csv(f"{RESULTS}/disease_specificity_{model}.csv", index=False)

    json.dump(out, open(f"{RESULTS}/cross_disease_summary.json", "w"), indent=2)

    print("Interpretation guide:")
    print("  If prediction correlation markedly exceeds label overlap, the")
    print("  models impose a disease-agnostic ranking on top of whatever")
    print("  disease signal exists.")
    print("  If cross-prediction ratio is near 1.0, any disease's model")
    print("  finds another disease's targets about as well as its own.")
    print(f"\nwrote tables to {RESULTS}/", flush=True)


if __name__ == "__main__":
    main()
