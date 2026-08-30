#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Standalone analysis: build an "easy-feature -> thermoelectric target" panel from the
JARVIS dft_2d atlas artifacts, compute Spearman correlations, and produce
(1) a correlation CSV  (2) a 2D polyline correlation figure.

This is the scientific core that the multi-agent "Analyst" uses to form hypotheses.
Run with the te_manifold conda env.
"""
import json, os, sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = "/home/wangchao/work_wc/2D_ZT/jarvis_2d_te_atlas"             # jarvis_2d_te_atlas
OUT  = "/home/wangchao/work_wc/2D_ZT/matexplore/reports"
os.makedirs(OUT, exist_ok=True)

# ---- 1. load targets -------------------------------------------------------
zte = pd.read_csv(os.path.join(ROOT, "data", "processed", "ZT_e_all.csv"))
magpie = pd.read_parquet(os.path.join(ROOT, "features", "structure", "composition_magpie.parquet"))
elec = pd.read_parquet(os.path.join(ROOT, "features", "electronic", "electronic_features_v1.parquet"))

# snapshot for structural extras
with open(os.path.join(ROOT, "data", "raw", "jarvis", "dft_2d_snapshot.json")) as f:
    snap = json.load(f)
struct = []
for d in snap:
    a = d["attributes"]
    lv = a.get("lattice_vectors") or [[1,0,0],[0,1,0],[0,0,1]]
    lv = np.array(lv, dtype=float)
    # in-plane area = |a x b|, out-of-plane (vacuum) = c projection
    a1, a2, a3 = lv
    n_z = a3 / np.linalg.norm(a3)
    c = abs(np.dot(a3, n_z))
    area = abs(np.linalg.norm(np.cross(a1, a2)))
    struct.append({
        "jid": a.get("_jarvis_jid"),
        "formula": a.get("_jarvis_formula"),
        "density": a.get("_jarvis_density"),
        "formation_peratom": a.get("_jarvis_formation_energy_peratom"),
        "exfoliation": a.get("_jarvis_exfoliation_energy"),
        "ehull": a.get("_jarvis_ehull"),
        "nsites": a.get("_jarvis_nat"),
        "spg_number": a.get("_jarvis_spg_number"),
        "inplane_area": area,
        "vacuum_c": c,
        "aspect_c_over_sqrtA": c / np.sqrt(area) if area > 0 else np.nan,
    })
struct = pd.DataFrame(struct)
# -99999 -> NaN
for c in ["density","formation_peratom","exfoliation"]:
    struct[c] = pd.to_numeric(struct[c], errors="coerce").replace(-99999, np.nan)

# ---- 2. merge --------------------------------------------------------------
# zte has one row per (jid, carrier). pivot to wide so each jid has n & p targets.
w = zte.pivot_table(index="jid", columns="carrier",
                    values=["ZT_e","S_median","PF_mean","log_sigma_dom_geo"]).reset_index()
w.columns = ["jid"] + [f"{a}_{b}" for a,b in w.columns[1:]]
base = magpie.merge(elec, on="jid", how="left").merge(struct, on="jid", how="left").merge(w, on="jid", how="left")
base.to_csv(os.path.join(OUT, "feature_target_panel.csv"), index=False)

# ---- 3. feature panel (easy features only) ---------------------------------
FEATURES = [
    # electronic (still "cheap" relative to full transport; needed by hypothesis)
    ("Eg (band gap, eV)", "Eg_optb88vdw"),
    ("m* (median, m_e)", "m_star"),               # per-carrier: electron for n, hole for p
    # composition / magpie
    ("electronegativity mean", "electronegativity_mean"),
    ("electronegativity range", "electronegativity_range"),
    ("atomic mass mean", "atomic_mass_mean"),
    ("atomic mass max", "atomic_mass_max"),
    ("atomic radius mean", "atomic_radius_mean"),
    ("Z (atomic number) mean", "Z_mean"),
    ("ionization energy mean", "ionization_energy_mean"),
    ("electron affinity mean", "electron_affinity_mean"),
    ("row mean", "row_mean"),
    ("group range", "group_range"),
    # structural / thermodynamic
    ("density (g/cm3)", "density"),
    ("formation energy /atom (eV)", "formation_peratom"),
    ("exfoliation energy (eV)", "exfoliation"),
    ("n_sites", "nsites"),
    ("space group number", "spg_number"),
    ("in-plane area (Å^2)", "inplane_area"),
    ("vacuum / sqrt(area)", "aspect_c_over_sqrtA"),
]

# target list: (label, column, carrier)
def target_cols(carrier):
    return [
        ("ZT_e", f"ZT_e_{carrier}"),
        ("log10 PF", f"PF_mean_{carrier}"),     # will be log10-ed
        ("|S| (uV/K)", f"S_median_{carrier}"),  # will be abs
        ("log10 sigma", f"log_sigma_dom_geo_{carrier}"),
    ]

rows = []
for carrier in ["n", "p"]:
    d = base.copy()
    if carrier == "n":
        d["m_star"] = d["m_elec_median"]
    else:
        d["m_star"] = d["m_hole_median"]
    for feat_label, feat_col in FEATURES:
        col = feat_col
        for tlabel, tcol in target_cols(carrier):
            sub = d[[col, tcol]].dropna()
            if len(sub) < 50:
                continue
            x = sub[col].to_numpy(float)
            y = sub[tcol].to_numpy(float)
            if tlabel == "log10 PF":
                y = np.log10(np.clip(y, 1e-9, None))
            elif tlabel == "|S| (uV/K)":
                y = np.abs(y)
            rho = spearmanr(x, y).statistic
            if not np.isfinite(rho):
                continue
            rows.append(dict(carrier=carrier, feature=feat_label, target=tlabel,
                             spearman=float(rho), n=len(sub)))

corr = pd.DataFrame(rows)
corr.to_csv(os.path.join(OUT, "feature_target_correlation.csv"), index=False)

# ---- 4. 2D polyline figure --------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# order features by |mean correlation| for readability
feat_order = (corr.groupby("feature")["spearman"].apply(lambda s: s.abs().max())
              .sort_values(ascending=False).index.tolist())

carrier_color = {"n": "#1f77b4", "p": "#d62728"}
carrier_ls    = {"n": "-", "p": "--"}
targets = ["ZT_e", "log10 PF", "|S| (uV/K)", "log10 sigma"]

fig, axes = plt.subplots(2, 2, figsize=(15, 9))
for ax, tgt in zip(axes.ravel(), targets):
    sub = corr[corr.target == tgt]
    for c in ["n","p"]:
        cc = sub[sub.carrier == c].set_index("feature")["spearman"].reindex(feat_order)
        ax.plot(range(len(feat_order)), cc.to_numpy(), ls=carrier_ls[c],
                color=carrier_color[c], marker="o", ms=4, lw=1.4,
                label=f"{c}-type")
    ax.axhline(0, color="gray", lw=0.8)
    ax.axhline(0.5, color="gray", lw=0.6, ls=":"); ax.axhline(-0.5, color="gray", lw=0.6, ls=":")
    ax.set_xticks(range(len(feat_order)))
    ax.set_xticklabels([f.split("(")[0].strip() for f in feat_order], rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Spearman rho")
    ax.set_title(f"target: {tgt}")
    ax.set_ylim(-1.05, 1.05)
    ax.legend(fontsize=7); ax.grid(alpha=0.25)
fig.suptitle("Easy-feature  vs  thermoelectric-target  Spearman correlation  (JARVIS dft_2d, 1103 2D)",
             fontsize=12)
fig.tight_layout(rect=[0,0,1,0.97])
png = os.path.join(OUT, "feature_target_correlation_polyline.png")
fig.savefig(png, dpi=130)
print("saved:", png)
print("panel rows:", len(base), "feature-target pairs:", len(corr))
print(corr[corr.target=="ZT_e"].pivot_table(index="feature", columns="carrier", values="spearman").round(3).to_string())
