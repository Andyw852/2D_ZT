"""kappa_L verify: direct descriptor correlation + figures (interpretability layer)."""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

root = Path(__file__).resolve().parents[1]
df = pd.read_parquet(root / "features" / "kl_verify" / "kl_views.parquet")
n = len(df)
logk = np.log10(df["kL_300"].values.astype(float))

# build physics descriptors
from pymatgen.core import Composition
avg_mass = np.array([Composition(f).weight / (Composition(f).num_atoms if False else 1) for f in df["formula"]])
# average atomic mass (amu) = total weight / num atoms
avg_mass = np.array([Composition(f).weight for f in df["formula"]])
natoms = np.array([Composition(f).num_atoms for f in df["formula"]])
avg_mass = avg_mass / natoms  # per-atom average mass (amu)

df["avg_mass"] = avg_mass
df["v_s_proxy"] = np.sqrt(df["B_kv"].values.astype(float) * 1e9 / (df["density"].values.astype(float) * 1000))  # m/s (B in GPa)

feats = ["Eg_opt", "m_elec", "m_hole", "B_kv", "G_gv", "density", "avg_mass", "v_s_proxy"]
print("=== Spearman correlation of log10(kL@300K) with descriptors ===")
rows = []
for c in feats:
    s = df[[c]].copy()
    s["logk"] = logk
    s = s.dropna()
    rho, pv = stats.spearmanr(s[c], s["logk"])
    rp, pp = stats.pearsonr(s[c], s["logk"])
    rows.append({"feature": c, "N": int(len(s)), "spearman_rho": round(float(rho),3),
                 "spearman_p": float(pv), "pearson_r": round(float(rp),3)})
    print(f"  {c:12s} N={len(s):3d}  Spearman={rho:+.3f} (p={pv:.2e})  Pearson={rp:+.3f}")
pd.DataFrame(rows).to_csv(root / "data" / "audit" / "kl_descriptor_corr.csv", index=False)

# ---- figures ----
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig_dir = root / "figures"
fig_dir.mkdir(exist_ok=True)

# fig 1: log kL vs sound-velocity proxy + density (2-panel)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
d = df.dropna(subset=["v_s_proxy"])
ax[0].scatter(d["v_s_proxy"], np.log10(d["kL_300"]), c=d["density"], cmap="viridis", s=40)
ax[0].set_xlabel("sound-velocity proxy  sqrt(B/rho)  [m/s]")
ax[0].set_ylabel("log10(kL @300K)  [W/m/K]")
ax[0].set_title("kL vs elastic/sound-velocity proxy")
rho,_ = stats.spearmanr(d["v_s_proxy"], np.log10(d["kL_300"]))
ax[0].text(0.05, 0.9, f"Spearman={rho:+.2f}, N={len(d)}", transform=ax[0].transAxes)

d2 = df.dropna(subset=["density"])
ax[1].scatter(d2["density"], np.log10(d2["kL_300"]), c=d2["B_kv"], cmap="plasma", s=40)
ax[1].set_xlabel("density [g/cm^3]")
ax[1].set_ylabel("log10(kL @300K)")
ax[1].set_title("kL vs density")
rho2,_ = stats.spearmanr(d2["density"], np.log10(d2["kL_300"]))
ax[1].text(0.05, 0.9, f"Spearman={rho2:+.2f}, N={len(d2)}", transform=ax[1].transAxes)
plt.tight_layout()
plt.savefig(fig_dir / "kl_descriptor_scatter.png", dpi=130)
print("saved fig: kl_descriptor_scatter.png")

# fig 2: pairwise distance correlation (structure-dist vs kL-dist)
from graph_utils import hellinger_distance, soap_distance
soap_mean = np.load(root / "data" / "processed" / "kl_soap_geo.npy")
frac = np.load(root / "data" / "processed" / "kl_comp_frac.npy")
d_geo = soap_distance(soap_mean); d_geo /= d_geo.max()
d_comp = hellinger_distance(frac); d_comp /= d_comp.max()
d_struct = 0.5*d_geo + 0.5*d_comp
d_kL = np.abs(logk[:, None] - logk[None, :])
iu = np.triu_indices(n, k=1)
x, y = d_struct[iu], d_kL[iu]
fig, ax = plt.subplots(figsize=(5.2, 4.2))
ax.scatter(x, y, s=2, alpha=0.3)
rho,_ = stats.spearmanr(x, y)
ax.set_xlabel("structure distance (geo-SOAP + composition)")
ax.set_ylabel("kL distance  |log10 kL_i - log10 kL_j|")
ax.set_title(f"structure vs kL pairwise distance  (Spearman={rho:+.2f})")
plt.tight_layout()
plt.savefig(fig_dir / "kl_struct_vs_kL_dist.png", dpi=130)
print("saved fig: kl_struct_vs_kL_dist.png")
