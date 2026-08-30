"""L0-E: effective-mass 主值谱审计。"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

root = Path(__file__).resolve().parents[1]
raw = json.loads((root / "data" / "processed" / "transport_eigenvalues_raw.json").read_text(encoding="utf-8"))

rows = []
for jid in raw:
    if "electron_mass_300K" not in raw[jid] and "hole_mass_300K" not in raw[jid]:
        continue
    row = {"jid": jid}
    for f in ["electron_mass_300K", "hole_mass_300K"]:
        if f in raw[jid]:
            v = np.array([float(x) for x in raw[jid][f]])
            vs = np.sort(v)
            row[f"{f}_raw"] = v
            row[f"{f}_1"] = vs[0]; row[f"{f}_2"] = vs[1]; row[f"{f}_3"] = vs[2]
            row[f"{f}_mean"] = v.mean()
            row[f"{f}_median"] = np.median(v)
            row[f"{f}_min"] = vs[0]
            row[f"{f}_dom_geo"] = np.sqrt(vs[0] * vs[1])  # 两个较小(物理面内)主通道几何均值
            row[f"{f}_range"] = vs[2] - vs[0]
            row[f"{f}_spectral_ratio"] = np.log10(vs[2] / vs[0]) if vs[0] > 0 else np.nan
    rows.append(row)
df = pd.DataFrame(rows)
df.to_csv(root / "data" / "audit" / "effective_mass_spectrum.csv", index=False)

for f, lab in [("electron_mass_300K", "electron"), ("hole_mass_300K", "hole")]:
    v = df[f"{f}_mean"].dropna()
    med = df[f"{f}_median"].dropna()
    sr = df[f"{f}_spectral_ratio"].dropna()
    neg = int((df[f"{f}_1"].dropna() < 0).sum())
    print(f"\n=== {lab} effective mass (n={df[f'{f}_mean'].notna().sum()}) ===")
    print(f"  mean: median={v.median():.4f} p10={v.quantile(0.1):.4f} p90={v.quantile(0.9):.4f}")
    print(f"  median: median={med.median():.4f} p10={med.quantile(0.1):.4f} p90={med.quantile(0.9):.4f}")
    print(f"  spectral_ratio (log10 max/min): median={sr.median():.3f} p10={sr.quantile(0.1):.3f} p90={sr.quantile(0.9):.3f}")
    print(f"  negative eigenvalues: {neg}")
    # mean 被大值污染程度：median 与 mean 的比值
    ratio = (v / med.clip(lower=1e-9)).median()
    print(f"  median(mean/median) = {ratio:.3f} (>>1 说明 mean 被 out-of-plane-like 大值严重污染)")

# 图
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
for ax, f, lab in [(axes[0], "electron_mass_300K", "electron"), (axes[1], "hole_mass_300K", "hole")]:
    sr = df[f"{f}_spectral_ratio"].dropna()
    ax.hist(sr, bins=50, color="#4477aa")
    ax.set_xlabel("log10(m_max/m_min)")
    ax.set_title(f"{lab} effective-mass spectral ratio")
fig.tight_layout(); fig.savefig(root / "figures" / "effective_mass_spectrum.png", dpi=130); plt.close(fig)
print("\nwrote effective_mass_spectrum.csv + figure")
