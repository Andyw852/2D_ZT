"""Phase F - Step 7: Property Coverage Audit。

根据真实下载的 1103 条 dft_2d 数据统计每个物性的 N_available / coverage / data form，
并按第 17 节阈值给出 A/B/C/SKIPPED 决策：
- A 类（主要层）：coverage >= 50% 且单位/条件明确
- B 类（partial 层）：N >= 100 但 coverage < 50%
- C 类（exploratory）：N < 100
- SKIPPED_NOT_AVAILABLE：N = 0

同时做热电相关字段关键词搜索（第 15 节）。
"""
import json
import math
from pathlib import Path

import pandas as pd

root = Path(__file__).resolve().parents[1]
audit = root / "data" / "audit"
reports = root / "reports"
audit.mkdir(exist_ok=True)
reports.mkdir(exist_ok=True)

records = json.loads((root / "data" / "raw" / "jarvis" / "dft_2d_snapshot.json").read_text(encoding="utf-8"))
attrs = [r["attributes"] for r in records]
N = len(records)

MISSING_STR = {"", "na", "n/a", "none", "null", "not available", "nan"}
def is_missing(v):
    if v is None:
        return True
    if isinstance(v, float):
        if math.isnan(v) or v == -99999.0:
            return True
    if isinstance(v, int) and not isinstance(v, bool):
        if v == -99999:
            return True
    if isinstance(v, str):
        if v.strip().lower() in MISSING_STR:
            return True
    if isinstance(v, (list, tuple, dict)):
        if len(v) == 0:
            return True
    return False

def n_available(fields, mode="all"):
    """mode='all': 所有 field 都非缺失；mode='any': 任一 field 非缺失。"""
    cnt = 0
    for a in attrs:
        vals = [a.get(f) for f in fields]
        ok = all(not is_missing(v) for v in vals) if mode == "all" else any(not is_missing(v) for v in vals)
        if ok:
            cnt += 1
    return cnt

def decide(n, coverage):
    if n == 0:
        return "SKIPPED_NOT_AVAILABLE"
    if coverage >= 0.5:
        return "A (main layer)"
    if n >= 100:
        return "B (partial layer)"
    return "C (exploratory)"

rows = [
    # property, fields, mode, data_form, note
    ("Structure", ["lattice_vectors", "cartesian_site_positions", "species_at_sites"], "all", "structure", "lattice + cartesian coords + species (100%)"),
    ("Band gap (OptB88vdW)", ["_jarvis_optb88vdw_bandgap"], "all", "scalar (eV)", "100% incl. 0.0 = metal"),
    ("Band gap (MBJ/TBmBJ)", ["_jarvis_mbj_bandgap"], "all", "scalar (eV)", ""),
    ("Band gap (HSE06)", ["_jarvis_hse_gap"], "all", "scalar (eV)", ""),
    ("Effective mass (electron)", ["_jarvis_avg_elec_mass"], "all", "scalar", "3 eigenvalues mean, units UNVERIFIED"),
    ("Effective mass (hole)", ["_jarvis_avg_hole_mass"], "all", "scalar", "3 eigenvalues mean, units UNVERIFIED"),
    ("n-Seebeck", ["_jarvis_n-Seebeck"], "all", "scalar (muV/K)", "3 eigenvalues mean, T=600K n=1e20"),
    ("p-Seebeck", ["_jarvis_p-Seebeck"], "all", "scalar (muV/K)", "3 eigenvalues mean, T=600K n=1e20"),
    ("n-Power factor", ["_jarvis_n-powerfact"], "all", "scalar", "3 eigenvalues mean, T=600K n=1e20, const tau"),
    ("p-Power factor", ["_jarvis_p-powerfact"], "all", "scalar", "3 eigenvalues mean, T=600K n=1e20, const tau"),
    ("n-Conductivity", ["_jarvis_ncond"], "all", "scalar", "3 eigenvalues mean, T=600K n=1e20, const tau"),
    ("p-Conductivity", ["_jarvis_pcond"], "all", "scalar", "3 eigenvalues mean, T=600K n=1e20, const tau"),
    ("kappa_e (n, electronic)", ["_jarvis_nkappa"], "all", "scalar", "3 eigenvalues mean, T=600K n=1e20, raw kappa/tau"),
    ("kappa_e (p, electronic)", ["_jarvis_pkappa"], "all", "scalar", "3 eigenvalues mean, T=600K n=1e20, raw kappa/tau"),
    ("Dielectric (static epsx/y/z)", ["_jarvis_epsx", "_jarvis_epsy", "_jarvis_epsz"], "any", "scalar x3", ""),
    ("Dielectric (HF mepsx/y/z)", ["_jarvis_mepsx", "_jarvis_mepsy", "_jarvis_mepsz"], "any", "scalar x3", ""),
    ("Exfoliation energy", ["_jarvis_exfoliation_energy"], "all", "scalar", ""),
    ("E_hull", ["_jarvis_ehull"], "all", "scalar (eV/atom)", "all 0.0 (on-hull)"),
    ("Formation energy/atom", ["_jarvis_formation_energy_peratom"], "all", "scalar (eV/atom)", ""),
    ("Elastic tensor", ["_jarvis_elastic_tensor"], "all", "tensor", "NOT exported by OPTIMADE"),
    ("kappa_L (lattice)", [], "all", "scalar", "NOT in JARVIS dft_2d"),
    ("ZT", [], "all", "scalar", "NOT in JARVIS dft_2d (would be target-leak)"),
]

out_rows = []
for prop, fields, mode, form, note in rows:
    if not fields:
        n = 0
    else:
        n = n_available(fields, mode)
    coverage = n / N if N else 0.0
    d = decide(n, coverage)
    out_rows.append({
        "Property": prop,
        "N_available": n,
        "N_total": N,
        "Coverage": round(coverage, 4),
        "Data_form": form,
        "Decision": d,
        "Note": note,
    })

df = pd.DataFrame(out_rows)
csv = audit / "property_coverage.csv"
df.to_csv(csv, index=False)
print(f"Wrote {csv}")
print()
print(df.to_string(index=False))

# ---- 热电相关字段关键词搜索（第 15 节）----
KEYWORDS = ["seebeck", "power", "powerfact", "conduct", "conductivity", "sigma", "thermal",
            "kappa", "zt", "mass", "effective", "elastic", "gap", "bandgap", "boltz",
            "carrier", "mobility", "exfol", "ehull", "formation"]
all_fields = sorted({k for a in attrs for k in a.keys()})
te_lines = []
for kw in KEYWORDS:
    matches = [f for f in all_fields if kw in f.lower()]
    te_lines.append(f"[{kw}] -> {matches}")
te_txt = "\n".join(te_lines)
(reports / "te_candidate_fields.txt").write_text(
    "Thermoelectric-related field keyword search (all attribute field names):\n\n" + te_txt + "\n",
    encoding="utf-8",
)
print()
print("Thermoelectric candidate fields written to reports/te_candidate_fields.txt")
print(te_txt)
