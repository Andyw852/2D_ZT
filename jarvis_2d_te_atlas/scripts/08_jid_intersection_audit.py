"""Phase H: JID 交集审计 + 数值质量/符号检查。"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

root = Path(__file__).resolve().parents[1]
raw = json.loads((root / "data" / "processed" / "transport_eigenvalues_raw.json").read_text(encoding="utf-8"))
opt = json.loads((root / "data" / "raw" / "jarvis" / "dft_2d_snapshot.json").read_text(encoding="utf-8"))
opt_by_jid = {r["attributes"]["_jarvis_jid"]: r["attributes"] for r in opt}

def fset(field):
    return {j for j in raw if field in raw[j] and len(raw[j][field]) == 3}

sets = {
    "S_n": fset("nseeb"),
    "S_p": fset("pseeb"),
    "sigma_n": fset("ncond"),
    "sigma_p": fset("pcond"),
    "kappa_e_n": fset("nkappa"),
    "kappa_e_p": fset("pkappa"),
    "PF_n": fset("npf"),
    "PF_p": fset("ppf"),
    "mstar": fset("electron_mass_300K"),
}
# D_Eg: 所有有 optb88vdw_bandgap 的材料
sets["Eg"] = {j for j in opt_by_jid if opt_by_jid[j].get("_jarvis_optb88vdw_bandgap") not in (None, -99999, -99999.0)}

print("=== JID set sizes ===")
for k, v in sets.items():
    print(f"  {k}: {len(v)}")

# 关键交集
print("\n=== 关键交集 ===")
D_tn = sets["S_n"] & sets["sigma_n"] & sets["kappa_e_n"]
D_tp = sets["S_p"] & sets["sigma_p"] & sets["kappa_e_p"]
print(f"|S_n ∩ sigma_n| = {len(sets['S_n'] & sets['sigma_n'])}")
print(f"|S_n ∩ kappa_e_n| = {len(sets['S_n'] & sets['kappa_e_n'])}")
print(f"|sigma_n ∩ kappa_e_n| = {len(sets['sigma_n'] & sets['kappa_e_n'])}")
print(f"|S_n ∩ sigma_n ∩ kappa_e_n| = {len(D_tn)}")
print(f"|S_p ∩ sigma_p ∩ kappa_e_p| = {len(D_tp)}")
print(f"|mstar ∩ D_transport_n| = {len(sets['mstar'] & D_tn)}")
print(f"|mstar ∩ D_transport_p| = {len(sets['mstar'] & D_tp)}")
print(f"|Eg ∩ D_transport_n| = {len(sets['Eg'] & D_tn)}")
print(f"|Eg ∩ D_transport_p| = {len(sets['Eg'] & D_tp)}")

# overlap matrix
names = list(sets.keys())
rows = []
for a in names:
    for b in names:
        inter = len(sets[a] & sets[b])
        rows.append({
            "View_A": a, "View_B": b,
            "N_A": len(sets[a]), "N_B": len(sets[b]),
            "N_intersection": inter,
            "Fraction_A": round(inter / len(sets[a]), 4) if sets[a] else 0,
            "Fraction_B": round(inter / len(sets[b]), 4) if sets[b] else 0,
        })
odf = pd.DataFrame(rows)
odf.to_csv(root / "data" / "audit" / "view_jid_overlap.csv", index=False)
print(f"\nwrote data/audit/view_jid_overlap.csv ({len(odf)} rows)")

# 数值质量检查（对每个 property，n/p，基于 3 本征值的 mean）
def numstats(vals):
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    return {
        "count": len(vals),
        "nan_inf": 0,
        "zero": int((vals == 0).sum()),
        "positive": int((vals > 0).sum()),
        "negative": int((vals < 0).sum()),
        "min": vals.min(),
        "p1": np.percentile(vals, 1),
        "median": np.median(vals),
        "p99": np.percentile(vals, 99),
        "max": vals.max(),
    }

FIELD_MEAN = {"nseeb": "S_n", "pseeb": "S_p", "ncond": "sigma_n", "pcond": "sigma_p",
              "nkappa": "kappa_e_n", "pkappa": "kappa_e_p", "npf": "PF_n", "ppf": "PF_p"}
naudit = []
for f, label in FIELD_MEAN.items():
    means = []
    for j in raw:
        if f in raw[j]:
            e = [float(x) for x in raw[j][f]]
            means.append(np.mean(e))
    st = numstats(means)
    st["property"] = label
    st["field"] = f
    naudit.append(st)
ndf = pd.DataFrame(naudit)
cols = ["property", "field", "count", "zero", "positive", "negative", "min", "p1", "median", "p99", "max"]
ndf = ndf[cols]
ndf.to_csv(root / "data" / "audit" / "transport_numeric_audit.csv", index=False)
print(f"\n=== 数值质量审计（基于 eigenvalue mean）===")
print(ndf.to_string(index=False))

# 符号检查（第 23-26 节）
print("\n=== Seebeck 符号检查 ===")
for f, label, expect in [("nseeb", "n-Seebeck", "negative"), ("pseeb", "p-Seebeck", "positive")]:
    viol_mean = 0
    viol_any = 0
    n = 0
    for j in raw:
        if f in raw[j]:
            n += 1
            e = [float(x) for x in raw[j][f]]
            m = np.mean(e)
            if expect == "negative" and m >= 0:
                viol_mean += 1
            if expect == "positive" and m <= 0:
                viol_mean += 1
            # any eigenvalue violating
            if expect == "negative" and any(v >= 0 for v in e):
                viol_any += 1
            if expect == "positive" and any(v <= 0 for v in e):
                viol_any += 1
    print(f"  {label}: n={n}, mean-sign violation={viol_mean}, any-eigenvalue violation={viol_any}")

print("\n=== 正值检查（sigma/kappa_e/PF，eigenvalue 级）===")
for f, label in [("ncond","sigma_n"),("pcond","sigma_p"),("nkappa","kappa_e_n"),("pkappa","kappa_e_p"),("npf","PF_n"),("ppf","PF_p")]:
    neg = 0; zero = 0; n = 0
    for j in raw:
        if f in raw[j]:
            e = [float(x) for x in raw[j][f]]
            n += 3
            neg += sum(1 for v in e if v < 0)
            zero += sum(1 for v in e if v == 0)
    print(f"  {label}: total eigenvalues={n}, negative={neg}, zero={zero}")
