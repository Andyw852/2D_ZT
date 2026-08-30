"""L0-C: Seebeck 主值谱 + 符号一致性审计。"""
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
opt = json.loads((root / "data" / "raw" / "jarvis" / "dft_2d_snapshot.json").read_text(encoding="utf-8"))
opt_by_jid = {r["attributes"]["_jarvis_jid"]: r["attributes"] for r in opt}

def gap(jid):
    a = opt_by_jid.get(jid, {})
    g = a.get("_jarvis_optb88vdw_bandgap")
    return g if g not in (None, -99999, -99999.0) else np.nan

def analyze(carrier):
    f = "nseeb" if carrier == "n" else "pseeb"
    expect_pos = (carrier == "p")  # p 期望正，n 期望负
    rows = []
    for jid in raw:
        if f not in raw[jid]:
            continue
        S = np.array([float(x) for x in raw[jid][f]])
        med = np.median(S); mean = S.mean()
        mad = np.median(np.abs(S - med))
        n_expected = int(((S > 0) if expect_pos else (S < 0)).sum())
        # majority sign: at least 2 of 3 have same sign
        pos_cnt = int((S > 0).sum()); neg_cnt = int((S < 0).sum())
        if pos_cnt >= 2:
            majority_sign = 1
        elif neg_cnt >= 2:
            majority_sign = -1
        else:
            majority_sign = 0
        rows.append({
            "jid": jid, "S_mean": mean, "S_median": med, "S_std": S.std(),
            "S_MAD": mad, "S_min": S.min(), "S_max": S.max(), "S_range": S.max() - S.min(),
            "S_abs_mean": np.abs(S).mean(), "S_abs_max": np.abs(S).max(),
            "S_relative_spread": S.std() / (abs(mean) + 1e-12),
            "S_1": S[0], "S_2": S[1], "S_3": S[2],
            "N_expected_sign": n_expected, "sign_fraction": n_expected / 3,
            "majority_sign": majority_sign,
            "mean_matches_majority": int(np.sign(mean) == majority_sign),
            "median_matches_majority": int(np.sign(med) == majority_sign),
            "gap": gap(jid),
        })
    df = pd.DataFrame(rows)
    # 分类
    def cls(r):
        n = r["N_expected_sign"]
        if carrier == "n":
            return {3: "N3", 2: "N2", 1: "N1", 0: "N0"}[n]
        else:
            return {3: "P3", 2: "P2", 1: "P1", 0: "P0"}[n]
    df["class"] = df.apply(cls, axis=1)
    return df

for carrier in ["n", "p"]:
    df = analyze(carrier)
    out = root / "data" / "audit" / f"seebeck_spectrum_{carrier}.csv"
    df.to_csv(out, index=False)
    print(f"\n=== {carrier}-type Seebeck (n={len(df)}) ===")
    print("  sign consistency class counts:", df["class"].value_counts().sort_index().to_dict())
    print("  sign_fraction mean:", round(df["sign_fraction"].mean(), 3))
    # mean vs median 稳健性
    mm = df["mean_matches_majority"].mean(); md = df["median_matches_majority"].mean()
    print(f"  mean matches majority sign: {mm:.3f}, median matches majority sign: {md:.3f}")
    # sign consistency vs band gap
    g = df.dropna(subset=["gap"])
    for cls, grp in g.groupby("class"):
        print(f"    {cls}: n={len(grp)} gap median={grp['gap'].median():.3f} mean={grp['gap'].mean():.3f}")
    # 严重违例 (N0/N1 or P0/P1) vs 正常 (N3/N2 or P3/P2) 的 gap 比较
    bad = g[g["class"].isin(["N0","N1","P0","P1"])]
    good = g[g["class"].isin(["N3","N2","P3","P2"])]
    if len(bad) and len(good):
        u = stats.mannwhitneyu(bad["gap"], good["gap"], alternative="less")
        print(f"  gap: bad-sign median={bad['gap'].median():.3f} vs good-sign median={good['gap'].median():.3f}, Mann-Whitney p={u.pvalue:.2e}")
    # 图
    fig, ax = plt.subplots(figsize=(6, 5))
    order = sorted(g["class"].unique())
    ax.boxplot([g[g["class"]==c]["gap"].dropna() for c in order], tick_labels=order)
    ax.set_ylabel("OptB88vdW band gap (eV)"); ax.set_title(f"{carrier}-type sign consistency vs band gap")
    fig.tight_layout(); fig.savefig(root / "figures" / f"seebeck_sign_vs_gap_{carrier}.png", dpi=130); plt.close(fig)
    # mean vs median 散点
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(df["S_median"], df["S_mean"], s=4, alpha=0.4, color="#4477aa")
    lim = max(abs(df[["S_median","S_mean"]].max().max()), abs(df[["S_median","S_mean"]].min().min())) * 1.05
    ax.plot([-lim, lim], [-lim, lim], "r--", lw=0.7)
    ax.set_xlabel("S_median"); ax.set_ylabel("S_mean"); ax.set_title(f"{carrier}-type S_mean vs S_median")
    fig.tight_layout(); fig.savefig(root / "figures" / f"seebeck_mean_vs_median_{carrier}.png", dpi=130); plt.close(fig)

# 输出 sign consistency CSV（含 gap）
for carrier in ["n", "p"]:
    df = analyze(carrier)
    df.to_csv(root / "data" / "audit" / f"seebeck_sign_consistency_{carrier}.csv", index=False)
print("\nwrote seebeck spectrum + sign consistency CSVs")
