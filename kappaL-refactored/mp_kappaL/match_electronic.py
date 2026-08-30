"""从本地 JARVIS dft_3d 匹配 MP 材料的电子性质 (带隙 + 有效质量)。"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from pymatgen.core import Composition

root = Path(__file__).resolve().parents[1]
mp = root / "mp_kappaL"

# ---- JARVIS dft_3d -> reduced_formula -> 最稳定多晶型 ----
jar = json.load(open(root / "jarvis_2d_te_atlas" / "data" / "raw" / "external" / "jarvis_kl" / "jdft_3d-8-18-2021.json"))
print("jarvis records:", len(jar))

def num(v, sentinel=(-99999, "na")):
    if v is None or v == "na" or v == "nan":
        return np.nan
    try:
        f = float(v)
        return np.nan if f <= -99999 else f
    except (ValueError, TypeError):
        return np.nan

best = {}  # canon -> dict of best values
for x in jar:
    f = str(x.get("formula", "")).split(" JVASP-")[0].strip()
    if not f:
        continue
    try:
        canon = Composition(f).reduced_formula
    except Exception:
        continue
    e = x.get("formation_energy_peratom")
    fe = num(e)
    if np.isnan(fe):
        continue
    cur = best.get(canon)
    if cur is None or fe < cur["fe"]:
        best[canon] = {
            "fe": fe,
            "gap_opt": num(x.get("optb88vdw_bandgap")),
            "gap_mbj": num(x.get("mbj_bandgap")),
            "m_elec": num(x.get("avg_elec_mass")),
            "m_hole": num(x.get("avg_hole_mass")),
        }
print("jarvis unique canon:", len(best))

# ---- MP 材料 ----
meta = pd.read_parquet(mp / "processed" / "views_meta.parquet")
def canon_of(f):
    try:
        return Composition(f).reduced_formula
    except Exception:
        return None
meta["canon"] = meta["formula"].map(canon_of)
rows = []
for _, r in meta.iterrows():
    b = best.get(r["canon"])
    if b:
        rows.append({"material_id": r["material_id"], "canon": r["canon"],
                     "gap_opt": b["gap_opt"], "gap_mbj": b["gap_mbj"],
                     "m_elec": b["m_elec"], "m_hole": b["m_hole"]})
elec = pd.DataFrame(rows)
print("MP materials matched to JARVIS electronic:", len(elec))
print("gap_opt non-null:", elec["gap_opt"].notna().sum(), " gap_mbj non-null:", elec["gap_mbj"].notna().sum(),
      " m_elec non-null:", elec["m_elec"].notna().sum())
elec.to_parquet(mp / "processed" / "electronic_jarvis.parquet")
print("saved electronic_jarvis.parquet")
