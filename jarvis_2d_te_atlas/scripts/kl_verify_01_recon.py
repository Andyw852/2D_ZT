"""κL 验证：侦察 starrydata2 实验 κL ↔ JARVIS dft_3d 的重叠情况。"""
import json, csv, sys, ast
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np

root = Path(__file__).resolve().parents[1]
starry = root / "data" / "raw" / "external" / "starrydata2"
dft3d = root / "data" / "raw" / "external" / "jarvis_kl" / "jdft_3d-8-18-2021.json"

def parse_xy(v):
    """x/y 字段是 JSON 风格数组字符串或单个数字。"""
    v = (v or "").strip()
    if not v:
        return []
    try:
        d = json.loads(v)
    except Exception:
        d = ast.literal_eval(v) if v[0] in "[(" else None
        if d is None:
            try:
                d = float(v)
            except Exception:
                return []
    if isinstance(d, (int, float)):
        return [d]
    return [float(t) for t in d]

# ---- 1. starrydata2 curves ----
curves = starry / "ThermoelectricMaterials_curves.csv"
kappa_rows = defaultdict(list)  # composition -> [(T, kL)]
total_rows = defaultdict(list)
with open(curves, encoding="utf-8", errors="ignore") as f:
    r = csv.DictReader(f)
    for row in r:
        py = (row.get("prop_y", "") or "").strip()
        comp = (row.get("composition", "") or "").strip()
        if not comp:
            continue
        xs = parse_xy(row.get("x")); ys = parse_xy(row.get("y"))
        if len(xs) != len(ys) or not xs:
            continue
        if py == "Lattice thermal conductivity":
            kappa_rows[comp].append((np.asarray(xs), np.asarray(ys)))
        elif py == "Thermal conductivity":
            total_rows[comp].append((np.asarray(xs), np.asarray(ys)))

def at_T(curves_y, T=300, tol=50):
    """每个 sample 在 T±tol 内取最近点，返回该 sample 的 kL；无可取点为 None。"""
    vals = []
    for xs, ys in curves_y:
        m = (xs >= T - tol) & (xs <= T + tol)
        if m.any():
            idx = int(np.argmin(np.abs(xs - T)))
            vals.append(float(ys[idx]))
    return vals

kL_300 = {}
for comp, cy in kappa_rows.items():
    v = at_T(cy, 300, 50)
    if v:
        kL_300[comp] = float(np.median(v))

print("compositions with kL curves:", len(kappa_rows), " total samples:", sum(len(v) for v in kappa_rows.values()))
print("compositions with kL@300K:", len(kL_300))

# ---- 2. dft_3d ----
d3 = json.load(open(dft3d))
jid_formula = {rec["jid"]: rec.get("formula", "") for rec in d3}
print("dft_3d records:", len(d3))

# ---- 3. normalization ----
from pymatgen.core import Composition
def canon(comp):
    try:
        return Composition(comp).reduced_formula
    except Exception:
        return None

sk = {canon(c): c for c in kL_300 if canon(c)}
dft_canon = defaultdict(list)
for jid, f in jid_formula.items():
    cf = canon(f)
    if cf:
        dft_canon[cf].append(jid)

print("starry kL canon formulas:", len(sk), " (of", len(kL_300), ")")
print("dft unique canon formulas:", len(dft_canon))
overlap = sorted(set(sk) & set(dft_canon))
print("=== OVERLAP ===", len(overlap))
print(overlap[:40])

out = root / "data" / "processed" / "kl_starry_compositions.json"
out.parent.mkdir(parents=True, exist_ok=True)
json.dump({"kL_300": kL_300, "overlap": overlap, "dft_canon_size": len(dft_canon)}, open(out, "w"), ensure_ascii=False)
print("saved:", out)
