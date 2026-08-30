"""L0-B: conductivity / kappa_e 主值谱二维性审计。

将 3 个本征值按数值排序 sigma_(1) <= sigma_(2) <= sigma_(3)（仅数值排序，非方向），
定义：
  sigma_dom_geo = sqrt(sigma2 * sigma3)   # 两个较强主通道的几何尺度
  A_total = log10(sigma3/sigma1)          # 整体谱各向异性
  D_sigma = log10(sigma2/sigma1)          # suppressed-channel contrast（含 quasi-2D 维度信息）
  A_dom   = log10(sigma3/sigma2)          # dominant-channel anisotropy
谱模式分类 C1/C2/C3/C4。
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

root = Path(__file__).resolve().parents[1]
raw = json.loads((root / "data" / "processed" / "transport_eigenvalues_raw.json").read_text(encoding="utf-8"))

def spectrum_audit(field):
    rows = []
    for jid in raw:
        if field not in raw[jid]:
            continue
        v = np.array([float(x) for x in raw[jid][field]])
        vs = np.sort(v)  # 升序：vs[0]<=vs[1]<=vs[2]
        s1, s2, s3 = vs
        rows.append({
            "jid": jid,
            f"{field}_1": s1, f"{field}_2": s2, f"{field}_3": s3,
            f"{field}_mean": v.mean(),
            f"{field}_median": np.median(v),
            f"{field}_dom_geo": np.sqrt(s2 * s3) if s2 > 0 and s3 > 0 else np.nan,
            f"{field}_geo_all": (s1 * s2 * s3) ** (1/3) if (s1 > 0 and s2 > 0 and s3 > 0) else np.nan,
            f"{field}_A_total": np.log10(s3 / s1) if s1 > 0 else np.nan,
            f"{field}_D": np.log10(s2 / s1) if s1 > 0 else np.nan,
            f"{field}_A_dom": np.log10(s3 / s2) if s2 > 0 else np.nan,
        })
    return pd.DataFrame(rows)

for carrier, f in [("n", "ncond"), ("p", "pcond")]:
    df = spectrum_audit(f)
    out = root / "data" / "audit" / f"conductivity_spectrum_audit_{carrier}.csv"
    df.to_csv(out, index=False)
    D = df[f"{f}_D"].dropna()
    Ad = df[f"{f}_A_dom"].dropna()
    At = df[f"{f}_A_total"].dropna()
    print(f"\n=== {carrier}-type sigma spectrum (n={len(df)}) ===")
    print(f"  D_sigma (log10 s2/s1): median={D.median():.3f} p10={D.quantile(0.1):.3f} p90={D.quantile(0.9):.3f}")
    print(f"  A_dom   (log10 s3/s2): median={Ad.median():.3f} p10={Ad.quantile(0.1):.3f} p90={Ad.quantile(0.9):.3f}")
    print(f"  A_total (log10 s3/s1): median={At.median():.3f}")
    # 谱模式分类（基于 D 与 A_dom 的中位数作为阈值）
    d_th = D.median(); a_th = Ad.median()
    def classify(row):
        d = row[f"{f}_D"]; a = row[f"{f}_A_dom"]
        if np.isnan(d) or np.isnan(a):
            return "C4"
        if d > d_th and a <= a_th:   # s1<<s2≈s3
            return "C1"
        if d > d_th and a > a_th:    # s1<s2<s3
            return "C2"
        if d <= d_th and a <= a_th:  # 近各向同性
            return "C3"
        return "C4"
    df["class"] = df.apply(classify, axis=1)
    print(f"  classification (th D={d_th:.3f}, A_dom={a_th:.3f}):", df["class"].value_counts().to_dict())
    # 图 D vs A_dom（joint non-NaN）
    m = df[[f"{f}_D", f"{f}_A_dom"]].dropna()
    fig, ax = plt.subplots(figsize=(6, 5.5))
    ax.scatter(m[f"{f}_D"], m[f"{f}_A_dom"], s=5, alpha=0.4, color="#4477aa")
    ax.set_xlabel("D_sigma = log10(sigma2/sigma1)  (suppressed-channel contrast)")
    ax.set_ylabel("A_sigma_dom = log10(sigma3/sigma2)  (dominant-channel anisotropy)")
    ax.set_title(f"{carrier}-type conductivity spectrum")
    ax.axvline(d_th, color="r", ls="--", lw=0.7); ax.axhline(a_th, color="r", ls="--", lw=0.7)
    fig.tight_layout(); fig.savefig(root / "figures" / f"sigma_D_vs_A_{carrier}.png", dpi=130); plt.close(fig)
    # 直方图 D
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(D, bins=60, color="#4477aa"); ax.set_xlabel("D_sigma"); ax.set_title(f"{carrier}-type D_sigma distribution")
    fig.tight_layout(); fig.savefig(root / "figures" / f"sigma_dimensionality_{carrier}.png", dpi=130); plt.close(fig)

# kappa_e 同处理，并与 sigma 比较 D/A
print("\n=== kappa_e spectrum + 与 sigma 的 D/A 相关性 ===")
for carrier, cf, kf in [("n", "ncond", "nkappa"), ("p", "pcond", "pkappa")]:
    cd = spectrum_audit(cf); kd = spectrum_audit(kf)
    cd.to_csv(root / "data" / "audit" / f"conductivity_spectrum_audit_{carrier}.csv", index=False)
    kd.to_csv(root / "data" / "audit" / f"kappa_spectrum_audit_{carrier}.csv", index=False)
    m = cd.merge(kd, on="jid", suffixes=("_sig", "_kap"))
    for key, (c, k) in {"D": (f"{cf}_D", f"{kf}_D"), "A_dom": (f"{cf}_A_dom", f"{kf}_A_dom")}.items():
        sub = m[[c, k]].dropna()
        rp = stats.pearsonr(sub[c], sub[k])[0]; rs = stats.spearmanr(sub[c], sub[k])[0]
        print(f"  [{carrier}] {key}: sigma vs kappa_e  Pearson={rp:.4f} Spearman={rs:.4f}")
    # kappa D vs A 图
    mk = kd[[f"{kf}_D", f"{kf}_A_dom"]].dropna()
    fig, ax = plt.subplots(figsize=(6, 5.5))
    ax.scatter(mk[f"{kf}_D"], mk[f"{kf}_A_dom"], s=5, alpha=0.4, color="#44aa77")
    ax.set_xlabel("D_kappa"); ax.set_ylabel("A_kappa_dom"); ax.set_title(f"{carrier}-type kappa_e spectrum")
    fig.tight_layout(); fig.savefig(root / "figures" / f"kappa_D_vs_A_{carrier}.png", dpi=130); plt.close(fig)
