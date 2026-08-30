"""kappa_L verify: build structure/electronic/kL views (consistent with existing 2D workflow).

Structure view = geometry-only SOAP (dummy species X, mean-pooled) + composition fraction (Hellinger),
exactly matching scripts 20/22. Electronic = Eg + effective masses. kL = starrydata2 experimental.
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict

root = Path(__file__).resolve().parents[1]
dft3d = root / "data" / "raw" / "external" / "jarvis_kl" / "jdft_3d-8-18-2021.json"
recon = json.load(open(root / "data" / "processed" / "kl_starry_compositions.json"))
kL_300 = recon["kL_300"]
overlap = recon["overlap"]

from pymatgen.core import Composition
def canon(f):
    try: return Composition(f).reduced_formula
    except Exception: return None

_kL_bucket = defaultdict(list)
for comp, v in kL_300.items():
    cf = canon(comp)
    if cf:
        _kL_bucket[cf].append(v)
canon2kL = {cf: float(np.median(vs)) for cf, vs in _kL_bucket.items()}

d3 = json.load(open(dft3d))
by_formula = defaultdict(list)
for rec in d3:
    cf = canon(rec.get("formula", ""))
    if not cf: continue
    by_formula[cf].append(rec)

MISSING = -99999.0
def finite(v):
    try: v = float(v)
    except (TypeError, ValueError): return np.nan
    return np.nan if v <= MISSING + 1 else v

rows = []
for cf in overlap:
    recs = by_formula.get(cf, [])
    if not recs: continue
    fe = [ (finite(r.get("formation_energy_peratom")) if np.isfinite(finite(r.get("formation_energy_peratom"))) else 1e9) for r in recs ]
    rec = recs[int(np.argmin(fe))]
    rows.append({
        "jid": rec["jid"],
        "formula": rec.get("formula", cf),
        "canon": cf,
        "kL_300": canon2kL.get(cf, np.nan),
        "Eg_opt": finite(rec.get("optb88vdw_bandgap")),
        "Eg_mbj": finite(rec.get("mbj_bandgap")),
        "m_elec": finite(rec.get("avg_elec_mass")),
        "m_hole": finite(rec.get("avg_hole_mass")),
        "B_kv": finite(rec.get("bulk_modulus_kv")),
        "G_gv": finite(rec.get("shear_modulus_gv")),
        "density": finite(rec.get("density")),
        "atoms": rec["atoms"],
    })

df = pd.DataFrame(rows)
print("rows:", len(df))

# ---------- structure: geometry-only SOAP (dummy X) ----------
from ase import Atoms as AseAtoms
from ase.data import atomic_numbers
from dscribe.descriptors import SOAP

def atoms_to_ase(a, dummy=True):
    elements = a["elements"]
    if dummy:
        elements = ["X"] * len(elements)
    return AseAtoms(symbols=elements, positions=np.asarray(a["coords"], float),
                     cell=np.asarray(a["lattice_mat"], float), pbc=True)

atoms_list = [atoms_to_ase(r["atoms"], dummy=True) for r in rows]
soap = SOAP(species=["X"], r_cut=6.0, n_max=6, l_max=6, sigma=1.0, periodic=True, average="off", dtype="float64")
X = soap.create(atoms_list, n_jobs=1)
soap_mean = np.array([x.mean(axis=0) for x in X])
print("geometry-only SOAP mean shape:", soap_mean.shape)

# ---------- composition fraction ----------
uniq_elems = sorted({sym for r in rows for sym in r["atoms"]["elements"]}, key=atomic_numbers.get)
print("elements (n=%d):" % len(uniq_elems), uniq_elems)
frac = np.zeros((len(rows), len(uniq_elems)))
for i, r in enumerate(rows):
    for sym in r["atoms"]["elements"]:
        frac[i, uniq_elems.index(sym)] += 1.0
frac = frac / frac.sum(axis=1, keepdims=True)

# ---------- save ----------
out = root / "features" / "kl_verify"
out.mkdir(parents=True, exist_ok=True)
df.drop(columns=["atoms"]).to_parquet(out / "kl_views.parquet")
np.save(root / "data" / "processed" / "kl_soap_geo.npy", soap_mean)
np.save(root / "data" / "processed" / "kl_comp_frac.npy", frac)
json.dump({"elements": uniq_elems, "n": len(df), "soap_shape": list(soap_mean.shape)},
          open(out / "meta.json", "w"), indent=2)
print("saved views:", out)
print("kL_300 non-null:", int(df["kL_300"].notna().sum()),
      "| m_elec non-null:", int(df["m_elec"].notna().sum()),
      "| B_kv non-null:", int(df["B_kv"].notna().sum()))
print(df[["jid","formula","kL_300","Eg_opt","m_elec","m_hole"]].head(8).to_string())
