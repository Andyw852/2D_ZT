"""MP 数据：构建 structure / elastic / kappa_L 视图 (不含带隙)。"""
import json
import numpy as np
from pathlib import Path

root = Path(__file__).resolve().parents[1]
mp = root / "mp_kappaL"
out = mp / "processed"
out.mkdir(parents=True, exist_ok=True)

el = json.load(open(mp / "raw" / "elasticity_all.json"))
print("docs:", len(el))

# ---- 过滤有效样本 ----
records = []
for x in el:
    sv = x.get("sound_velocity") or {}
    tc = x.get("thermal_conductivity") or {}
    bm = x.get("bulk_modulus") or {}
    gm = x.get("shear_modulus") or {}
    st = x.get("structure")
    if (st and tc.get("clarke") and tc.get("clarke") > 0
            and bm.get("vrh") and gm.get("vrh")
            and x.get("debye_temperature") and x.get("density")):
        records.append({
            "material_id": x["material_id"],
            "formula": x.get("formula_pretty"),
            "structure": st,
            "clarke": float(tc["clarke"]),
            "cahill": float(tc.get("cahill", np.nan)),
            "snyder_total": float(sv.get("snyder_total", np.nan)),
            "snyder_acoustic": float(sv.get("snyder_acoustic", np.nan)),
            "bulk_vrh": float(bm["vrh"]),
            "shear_vrh": float(gm["vrh"]),
            "debye": float(x["debye_temperature"]),
            "density": float(x["density"]),
            "nsites": int(x.get("nsites") or 0),
            "v_long": float(sv.get("longitudinal", np.nan)),
            "v_trans": float(sv.get("transverse", np.nan)),
        })
print("valid records:", len(records))

# ---- 结构 → ASE (dummy X) + 元素分数 ----
import ase
from ase import Atoms
from collections import Counter

elements_all = set()
atoms_list = []
comp_rows = []  # per-material element fraction dict
for r in records:
    st = r["structure"]
    lat = np.array(st["lattice"]["matrix"], dtype=float)
    pos = np.array([s["xyz"] for s in st["sites"]], dtype=float)
    elems = []
    frac = Counter()
    for s in st["sites"]:
        for sp in s["species"]:
            frac[sp["element"]] += float(sp.get("occu", 1.0))
    total = sum(frac.values()) or 1.0
    frac = {k: v / total for k, v in frac.items()}
    comp_rows.append(frac)
    elements_all.update(frac.keys())
    atoms_list.append(Atoms(symbols="X" * len(pos), positions=pos, cell=lat, pbc=True))

elems = sorted(elements_all)
print("n elements:", len(elems))

F = np.zeros((len(records), len(elems)))
for i, fr in enumerate(comp_rows):
    for e, v in fr.items():
        F[i, elems.index(e)] = v
np.save(out / "comp_frac.npy", F.astype(np.float32))

# ---- SOAP (dummy X) ----
from dscribe.descriptors import SOAP
soap = SOAP(species=["X"], periodic=True, r_cut=6.0, n_max=6, l_max=6, sigma=1.0, average="inner")
geo = soap.create(atoms_list, n_jobs=-1)
geo = np.asarray(geo, dtype=np.float32)
print("SOAP geo shape:", geo.shape)
np.save(out / "soap_geo.npy", geo)

# ---- 元数据 ----
import pandas as pd
meta = pd.DataFrame(records).drop(columns=["structure"])
meta.to_parquet(out / "views_meta.parquet")
meta.to_json(out / "views_meta.json")
np.save(out / "elem_basis.npy", np.array(elems))

# 汇总
print("saved. clarke n:", meta["clarke"].notna().sum(),
      "snyder_total n:", meta["snyder_total"].notna().sum())
print("clarke range:", meta["clarke"].min(), meta["clarke"].max())
print("bulk_vrh range:", meta["bulk_vrh"].min(), meta["bulk_vrh"].max())
