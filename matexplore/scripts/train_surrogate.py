#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Train cheap-feature -> thermoelectric-target surrogate (Ridge) on the JARVIS panel.
"Cheap" = composition (magpie) + structural descriptors, NO band-structure-derived
features (Eg, m*) — so it can score a brand-new POSCAR without any DFT.

Outputs: reports/surrogate_metrics.csv + skills/te-screen/models/ (JSON coeffs + scaler).
"""
import json, os
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
REP  = os.path.join(os.path.dirname(HERE), "reports")
SKILL_MODELS = os.path.join(os.path.dirname(HERE), "skills", "te-screen", "models")
os.makedirs(REP, exist_ok=True)
os.makedirs(SKILL_MODELS, exist_ok=True)

panel = pd.read_csv(os.path.join(REP, "feature_target_panel.csv"))

# cheap features: composition (magpie) + bare-lattice geometry.
# NO Eg/m*, NO DFT-derived (formation/exfoliation), NO spglib-dependent (spg_number).
# These 14 are reproducible from a bare POSCAR with element masses + a property table.
CHEAP = [
    "electronegativity_mean","electronegativity_range",
    "atomic_mass_mean","atomic_mass_max","atomic_radius_mean",
    "Z_mean","ionization_energy_mean","electron_affinity_mean",
    "row_mean","group_range",
    "density","nsites","inplane_area","aspect_c_over_sqrtA",
]

def run_target(carrier, target, ycol):
    d = panel.dropna(subset=CHEAP + [ycol]).copy()
    X = d[CHEAP].to_numpy(float)
    y = d[ycol].to_numpy(float)
    if target == "log10_PF":
        y = np.log10(np.clip(y, 1e-9, None))
    y = np.nan_to_num(y, nan=np.nanmedian(y))
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    preds = np.empty_like(y); preds[:] = np.nan
    for tr, te in kf.split(X):
        sc = StandardScaler().fit(X[tr])
        m = Ridge(alpha=1.0).fit(sc.transform(X[tr]), y[tr])
        preds[te] = m.predict(sc.transform(X[te]))
    rho = spearmanr(y, preds).statistic
    # refit on all data for deployment
    sc = StandardScaler().fit(X)
    m = Ridge(alpha=1.0).fit(sc.transform(X), y)
    blob = dict(carrier=carrier, target=target, ycol=ycol,
                features=CHEAP,
                feature_mean=sc.mean_.tolist(), feature_scale=sc.scale_.tolist(),
                coef=m.coef_.tolist(), intercept=float(m.intercept_),
                cv_spearman=float(rho), n=len(d))
    fn = f"{carrier}_{target}.json"
    with open(os.path.join(SKILL_MODELS, fn), "w") as f:
        json.dump(blob, f, indent=2)
    return blob

rows = []
for carrier, mass in [("n","m_elec_median"),("p","m_hole_median")]:
    for target, ycol in [("ZT_e", f"ZT_e_{carrier}"),
                         ("log10_PF", f"PF_mean_{carrier}")]:
        b = run_target(carrier, target, ycol)
        rows.append(dict(carrier=carrier, target=target, cv_spearman=b["cv_spearman"],
                         n=b["n"], n_features=len(CHEAP)))
        print(f"[{carrier}/{target}] CV Spearman = {b['cv_spearman']:+.3f}  (n={b['n']})")

# also a reference: how well do Eg alone and m* alone predict ZT_e (for the report)
print()
print("--- reference (band-structure features, single) ---")
d = panel.dropna(subset=["Eg_optb88vdw","ZT_e_n"])
print(f"Eg alone -> ZT_e(n): {spearmanr(d.Eg_optb88vdw, d.ZT_e_n).statistic:+.3f}")
d = panel.dropna(subset=["m_elec_median","ZT_e_n"])
print(f"m* alone -> ZT_e(n): {spearmanr(d.m_elec_median, d.ZT_e_n).statistic:+.3f}")

pd.DataFrame(rows).to_csv(os.path.join(REP, "surrogate_metrics.csv"), index=False)
print("\nsaved models to", SKILL_MODELS)
print("saved metrics to", os.path.join(REP, "surrogate_metrics.csv"))
