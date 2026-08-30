"""Phase J: 相关性 / 冗余 / PCA / anisotropy 审计。"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

root = Path(__file__).resolve().parents[1]
raw = json.loads((root / "data" / "processed" / "transport_eigenvalues_raw.json").read_text(encoding="utf-8"))

# ---- 载入特征 ----
n_df = pd.read_parquet(root / "features" / "transport" / "n_transport_tensor_features.parquet")
p_df = pd.read_parquet(root / "features" / "transport" / "p_transport_tensor_features.parquet")

def corr_table(df, carrier):
    base = ["S_mean", "log_sigma_mean", "log_kappa_e_mean", "PF_mean", "sigma_anisotropy_log", "kappa_e_anisotropy_log"]
    sub = df[base].dropna()
    pear = sub.corr(method="pearson")
    spear = sub.corr(method="spearman")
    return pear, spear, sub

pears, spears = {}, {}
for carrier, df in [("n", n_df), ("p", p_df)]:
    pear, spear, sub = corr_table(df, carrier)
    pears[carrier] = pear
    spears[carrier] = spear

# 输出 correlation 矩阵（合并 n/p，标注 carrier）
def stack(cd, name):
    rows = []
    for carrier, m in cd.items():
        for i, row in m.iterrows():
            for j in m.columns:
                rows.append({"carrier": carrier, "metric": name, "var1": i, "var2": j, "value": row[j]})
    return pd.DataFrame(rows)

pear_df = stack(pears, "pearson")
spear_df = stack(spears, "spearman")
pear_df.to_csv(root / "data" / "audit" / "transport_correlation_pearson.csv", index=False)
spear_df.to_csv(root / "data" / "audit" / "transport_correlation_spearman.csv", index=False)

print("=== 关键相关系数（Pearson / Spearman）===")
for carrier, df in [("n", n_df), ("p", p_df)]:
    sub = df[["S_mean","log_sigma_mean","log_kappa_e_mean","PF_mean"]].dropna()
    for (a,b) in [("log_sigma_mean","log_kappa_e_mean"), ("PF_mean","S_mean"), ("PF_mean","log_sigma_mean"), ("PF_mean","log_kappa_e_mean"), ("S_mean","log_sigma_mean")]:
        rp = stats.pearsonr(sub[a], sub[b])[0]
        rs = stats.spearmanr(sub[a], sub[b])[0]
        print(f"  [{carrier}] {a} vs {b}: pearson={rp:.4f} spearman={rs:.4f}")

# ---- PCA diagnostics（第 36 节）----
def pca_diag(df, cols, label):
    sub = df[cols].dropna()
    X = StandardScaler().fit_transform(sub.values)
    pca = PCA(n_components=min(len(cols), X.shape[0]))
    pca.fit(X)
    evr = pca.explained_variance_ratio_
    cum = np.cumsum(evr)
    # effective rank = 参与率 (participation ratio)
    eff = (evr.sum()**2) / (evr**2).sum()
    return {"label": label, "n": X.shape[0], "dims": len(cols)} |            {f"evr_{i+1}": round(evr[i], 6) for i in range(len(evr))} |            {"cum_evr_1": round(cum[0], 6), "cum_evr_2": round(cum[1], 6) if len(cum) > 1 else None,
            "effective_rank": round(eff, 3)}

mean_cols = ["S_mean", "log_sigma_mean", "log_kappa_e_mean"]
spectrum_cols = ["S_mean", "S_std", "S_range", "log_sigma_mean", "sigma_anisotropy_log", "log_kappa_e_mean", "kappa_e_anisotropy_log"]

pca_rows = []
for carrier, df in [("n", n_df), ("p", p_df)]:
    pca_rows.append(pca_diag(df, mean_cols, f"{carrier}_mean_only"))
    pca_rows.append(pca_diag(df, spectrum_cols, f"{carrier}_spectrum"))
pca_df = pd.DataFrame(pca_rows)
pca_df.to_csv(root / "data" / "audit" / "transport_pca_diagnostics.csv", index=False)
print("\n=== PCA 诊断 ===")
print(pca_df.to_string(index=False))

# ---- anisotropy 信息价值（第 37 节）：mean-only vs spectrum 的距离/kNN 比较 ----
def dist_rank_corr(df, cols, idx):
    # idx: 用于对齐的公共行（已经 dropna 的索引）
    sub = df.loc[idx, cols]
    X = StandardScaler().fit_transform(sub.values)
    D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=-1)
    iu = np.triu_indices(len(X), k=1)
    return D[iu]

def knn_overlap(df, cols_a, cols_b, k=10):
    common = df.dropna(subset=cols_a + cols_b)
    Xa = StandardScaler().fit_transform(common[cols_a].values)
    Xb = StandardScaler().fit_transform(common[cols_b].values)
    Da = np.linalg.norm(Xa[:, None, :] - Xa[None, :, :], axis=-1)
    Db = np.linalg.norm(Xb[:, None, :] - Xb[None, :, :], axis=-1)
    np.fill_diagonal(Da, np.inf); np.fill_diagonal(Db, np.inf)
    ka = np.argsort(Da, axis=1)[:, :k]
    kb = np.argsort(Db, axis=1)[:, :k]
    overlap = np.mean([len(set(ka[i]) & set(kb[i])) / k for i in range(len(common))])
    return overlap, len(common)

print("\n=== anisotropy 信息价值（第 37 节）===")
for carrier, df in [("n", n_df), ("p", p_df)]:
    common_idx = df.dropna(subset=mean_cols + spectrum_cols).index
    d_mean = dist_rank_corr(df, mean_cols, common_idx)
    d_spec = dist_rank_corr(df, spectrum_cols, common_idx)
    rho = stats.spearmanr(d_mean, d_spec)[0]
    ov, n = knn_overlap(df, mean_cols, spectrum_cols, k=10)
    print(f"  [{carrier}] n={n}: spearman(dist_mean, dist_spectrum)={rho:.4f}, kNN(10) overlap={ov:.4f}")

# ---- Seebeck 符号违规 JID 输出（第 23 节）----
viol = []
for jid in raw:
    for f, lab, expect in [("nseeb","n","negative"),("pseeb","p","positive")]:
        if f in raw[jid]:
            e = [float(x) for x in raw[jid][f]]
            m = np.mean(e)
            bad = (expect == "negative" and m >= 0) or (expect == "positive" and m <= 0)
            if bad:
                viol.append({"jid": jid, "carrier": lab, "S_mean": round(m,3),
                             "eigs": ",".join(f"{v:.2f}" for v in e)})
vdf = pd.DataFrame(viol)
vdf.to_csv(root / "data" / "audit" / "seebeck_sign_violations.csv", index=False)
print(f"\nSeebeck sign violations (mean-based): {len(vdf)} records -> data/audit/seebeck_sign_violations.csv")
