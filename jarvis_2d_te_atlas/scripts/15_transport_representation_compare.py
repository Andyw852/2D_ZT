"""L0-D + L0-F: R1/R2/R3 与 T1/T2/T3 transport representation 比较。"""
import json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from pathlib import Path

root = Path(__file__).resolve().parents[1]
raw = json.loads((root / "data" / "processed" / "transport_eigenvalues_raw.json").read_text(encoding="utf-8"))

def desc(jid, carrier):
    n = "n" if carrier == "n" else "p"
    sf, cf, kf, pf = n + "seeb", n + "cond", n + "kappa", n + "pf"
    if not all(k in raw[jid] for k in (sf, cf, kf, pf)):
        return None
    S = np.array([float(x) for x in raw[jid][sf]])
    C = np.array([float(x) for x in raw[jid][cf]])
    K = np.array([float(x) for x in raw[jid][kf]])
    PF = np.array([float(x) for x in raw[jid][pf]])
    C = np.sort(C); K = np.sort(K)
    med = np.median(S); mad = np.median(np.abs(S - med))
    d = {}
    d["S_mean"] = S.mean(); d["S_median"] = med; d["S_MAD"] = mad
    d["S_range"] = S.max() - S.min(); d["S_std"] = S.std()
    pos = (S > 0).sum() if carrier == "p" else (S < 0).sum()
    d["S_sign_fraction"] = pos / 3
    # sigma
    s1, s2, s3 = C
    d["log_sigma_mean"] = np.log10(C.mean() + 1e-6)
    d["log_sigma_dom_geo"] = np.log10(np.sqrt(s2 * s3) + 1e-6) if s2 > 0 else np.nan
    d["D_sigma"] = np.log10(s2 / s1) if s1 > 1e-9 else np.nan
    d["A_sigma_dom"] = np.log10(s3 / s2) if s2 > 1e-9 else np.nan
    # kappa
    k1, k2, k3 = K
    d["log_kappa_dom_geo"] = np.log10(np.sqrt(k2 * k3) + 1e-6) if k2 > 0 else np.nan
    d["D_kappa"] = np.log10(k2 / k1) if k1 > 1e-9 else np.nan
    d["A_kappa_dom"] = np.log10(k3 / k2) if k2 > 1e-9 else np.nan
    d["PF_mean"] = PF.mean()
    return d

def build(carrier):
    rows = []
    for jid in raw:
        d = desc(jid, carrier)
        if d:
            d["jid"] = jid
            rows.append(d)
    return pd.DataFrame(rows)

R1 = ["S_mean", "log_sigma_mean"]
R2 = ["S_median", "S_MAD", "S_range", "log_sigma_mean", "D_sigma", "A_sigma_dom"]
R3 = ["S_median", "S_MAD", "log_sigma_dom_geo", "A_sigma_dom", "D_sigma"]
T1 = ["S_median", "S_MAD", "S_range", "log_sigma_dom_geo", "D_sigma", "A_sigma_dom"]
T2 = ["S_median", "S_MAD", "S_range", "log_kappa_dom_geo", "D_kappa", "A_kappa_dom"]
T3 = T1 + ["log_kappa_dom_geo", "D_kappa", "A_kappa_dom"]

def std_mat(df, cols):
    sub = df.dropna(subset=cols)
    X = StandardScaler().fit_transform(sub[cols].values)
    return sub.index.values, X

def dist_rank(Xa, Xb):
    Da = np.linalg.norm(Xa[:, None, :] - Xa[None, :, :], axis=-1)
    Db = np.linalg.norm(Xb[:, None, :] - Xb[None, :, :], axis=-1)
    iu = np.triu_indices(len(Xa), k=1)
    return stats.spearmanr(Da[iu], Db[iu])[0]

def knn_overlap(Xa, Xb, k):
    Da = np.linalg.norm(Xa[:, None, :] - Xa[None, :, :], axis=-1)
    Db = np.linalg.norm(Xb[:, None, :] - Xb[None, :, :], axis=-1)
    np.fill_diagonal(Da, np.inf); np.fill_diagonal(Db, np.inf)
    ka = np.argsort(Da, axis=1)[:, :k]; kb = np.argsort(Db, axis=1)[:, :k]
    return np.mean([len(set(ka[i]) & set(kb[i])) / k for i in range(len(Xa))])

def pca_rank(X):
    pca = PCA().fit(StandardScaler().fit_transform(X))
    evr = pca.explained_variance_ratio_
    return (evr.sum()**2) / (evr**2).sum()

comp_rows = []
for carrier in ["n", "p"]:
    df = build(carrier)
    # R1/R2/R3 比较（在共同非空子集上）
    common = df.dropna(subset=R1 + R2 + R3).index
    idx = np.arange(len(df))
    for name, cols in [("R1", R1), ("R2", R2), ("R3", R3)]:
        sub = df.loc[common, cols].dropna()
        # 用共同的完整子集
    # 统一在 common 子集上
    sub_common = df.loc[common]
    mats = {}
    for name, cols in [("R1", R1), ("R2", R2), ("R3", R3)]:
        X = StandardScaler().fit_transform(sub_common[cols].values)
        mats[name] = X
    for a, b in [("R1","R2"), ("R1","R3"), ("R2","R3")]:
        dr = dist_rank(mats[a], mats[b])
        ov = {k: knn_overlap(mats[a], mats[b], k) for k in [5,10,20,30]}
        comp_rows.append({"carrier": carrier, "pair": f"{a} vs {b}", "dist_spearman": round(dr,4),
                          "knn5": round(ov[5],4), "knn10": round(ov[10],4), "knn20": round(ov[20],4), "knn30": round(ov[30],4)})
    # T1/T2/T3
    sub_t = df.dropna(subset=T1 + T2 + T3)
    matsT = {}
    for name, cols in [("T1", T1), ("T2", T2), ("T3", T3)]:
        matsT[name] = StandardScaler().fit_transform(sub_t[cols].values)
    for a, b in [("T1","T2"), ("T1","T3"), ("T2","T3")]:
        dr = dist_rank(matsT[a], matsT[b])
        ov10 = knn_overlap(matsT[a], matsT[b], 10)
        comp_rows.append({"carrier": carrier, "pair": f"{a} vs {b}", "dist_spearman": round(dr,4),
                          "knn5": None, "knn10": round(ov10,4), "knn20": None, "knn30": None})
    # PCA rank
    for name, cols in [("T1", T1), ("T2", T2), ("T3", T3)]:
        comp_rows.append({"carrier": carrier, "pair": f"{name} pca_rank", "dist_spearman": round(pca_rank(matsT[name]),3),
                          "knn5": None, "knn10": None, "knn20": None, "knn30": None})

comp = pd.DataFrame(comp_rows)
comp.to_csv(root / "data" / "audit" / "transport_representation_comparison_n.csv" if False else root / "data" / "audit" / "transport_T1_T2_T3_comparison.csv", index=False)
print(comp.to_string(index=False))

# 分 carrier 保存
for carrier in ["n", "p"]:
    sub = comp[comp["carrier"] == carrier]
    sub.to_csv(root / "data" / "audit" / f"transport_representation_comparison_{carrier}.csv", index=False)
print("\nwrote representation comparison CSVs")
