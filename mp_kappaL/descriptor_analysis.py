"""MP 数据: kL 与单个物理描述符的直接相关 (解释层) + 图。"""
import json
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

root = Path(__file__).resolve().parents[1]
mp = root / "mp_kappaL"

meta = pd.read_parquet(mp / "processed" / "views_meta.parquet")
elec = pd.read_parquet(mp / "processed" / "electronic_jarvis.parquet")
meta = meta.merge(elec, on="material_id", how="left")
comp_frac = np.load(mp / "processed" / "comp_frac.npy").astype(float)
elems = np.load(mp / "processed" / "elem_basis.npy", allow_pickle=True)

from pymatgen.core import Element
mass = np.array([Element(e).atomic_mass.real if hasattr(Element(e).atomic_mass,'real') else float(Element(e).atomic_mass) for e in elems])
avg_mass = (comp_frac * mass).sum(axis=1)
meta["avg_mass"] = avg_mass

# 声速 proxy: 体模量/密度 (已由 longitudinal 给出, 但也算一个 sqrt(B/rho))
meta["v_s_proxy"] = np.sqrt(meta["bulk_vrh"] * 1e9 / (meta["density"] * 1000))  # m/s

meta = meta[meta["bulk_vrh"] < 1000].copy()  # 去异常

logk = np.log10(meta["clarke"].values.astype(float))
feats = ["debye", "bulk_vrh", "shear_vrh", "density", "avg_mass",
         "v_long", "v_trans", "v_s_proxy", "gap_opt", "gap_mbj", "m_elec", "m_hole"]

print("=== Spearman correlation of log10(clarke) with descriptors ===")
rows = []
for c in feats:
    s = meta[[c]].copy(); s["logk"] = logk
    s = s.dropna()
    if len(s) < 30: continue
    rho, pv = stats.spearmanr(s[c], s["logk"])
    rows.append({"feature": c, "N": int(len(s)), "spearman": round(float(rho),3), "p": float(pv)})
    print(f"  {c:18s} N={len(s):5d} Spearman={rho:+.3f} (p={pv:.2e})")
pd.DataFrame(rows).to_csv(mp / "processed" / "descriptor_corr.csv", index=False)

# ---- 图 ----
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig_dir = mp / "figures"
fig_dir.mkdir(exist_ok=True)

fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
d = meta.dropna(subset=["v_s_proxy"])
ax[0].scatter(d["v_s_proxy"], np.log10(d["clarke"]), c=d["density"], cmap="viridis", s=6, alpha=0.5)
ax[0].set_xlabel("sound-velocity proxy sqrt(B/rho) [m/s]"); ax[0].set_ylabel("log10(clarke kL) [W/m/K]")
rho,_ = stats.spearmanr(d["v_s_proxy"], np.log10(d["clarke"]))
ax[0].set_title(f"kL vs sound velocity (Spearman={rho:+.2f}, N={len(d)})")

d2 = meta.dropna(subset=["debye"])
ax[1].scatter(d2["debye"], np.log10(d2["clarke"]), c=d2["bulk_vrh"], cmap="plasma", s=6, alpha=0.5)
ax[1].set_xlabel("Debye temperature [K]"); ax[1].set_ylabel("log10(clarke kL)")
rho2,_ = stats.spearmanr(d2["debye"], np.log10(d2["clarke"]))
ax[1].set_title(f"kL vs Debye temperature (Spearman={rho2:+.2f})")

d3 = meta.dropna(subset=["gap_opt"])
ax[2].scatter(d3["gap_opt"], np.log10(d3["clarke"]), s=6, alpha=0.4)
ax[2].set_xlabel("band gap (OptB88) [eV]"); ax[2].set_ylabel("log10(clarke kL)")
rho3,_ = stats.spearmanr(d3["gap_opt"], np.log10(d3["clarke"]))
ax[2].set_title(f"kL vs band gap (Spearman={rho3:+.2f})")
plt.tight_layout(); plt.savefig(fig_dir / "kL_descriptor_scatter.png", dpi=130)
print("saved figures/kL_descriptor_scatter.png")

# 结构距离 vs kL 距离散点 (用已保存的距离? 这里重算轻量版: 抽样)
fig2, ax = plt.subplots(figsize=(5.2, 4.4))
sys.path.insert(0, str(root / "jarvis_2d_te_atlas" / "scripts"))
from graph_utils import hellinger_distance, soap_distance
soap_geo = np.load(mp / "processed" / "soap_geo.npy").astype(float)
d_geo = soap_distance(soap_geo); d_geo /= d_geo.max()
d_comp = hellinger_distance(comp_frac); d_comp /= d_comp.max()
d_struct = 0.5*d_geo + 0.5*d_comp
logk_all = np.log10(meta["clarke"].values.astype(float))
d_kL = np.abs(logk_all[:,None] - logk_all[None,:])
n = len(meta)
rng = np.random.RandomState(0)
idx = rng.choice(n, size=4000, replace=False)
sub = d_struct[np.ix_(idx,idx)]; subk = d_kL[np.ix_(idx,idx)]
iu = np.triu_indices(4000, k=1)
ax.scatter(sub[iu], subk[iu], s=1, alpha=0.15)
rho,_ = stats.spearmanr(sub[iu], subk[iu])
ax.set_xlabel("structure distance"); ax.set_ylabel("|log10 kL_i - log10 kL_j|")
ax.set_title(f"structure vs kL pairwise distance (Spearman={rho:+.3f})")
plt.tight_layout(); plt.savefig(fig_dir / "struct_vs_kL_dist.png", dpi=130)
print("saved figures/struct_vs_kL_dist.png")
