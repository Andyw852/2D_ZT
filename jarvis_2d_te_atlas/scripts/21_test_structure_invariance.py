"""Phase L (L15-L16): SOAP invariance tests (A/B/C/D)。"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from ase import Atoms
from dscribe.descriptors import SOAP

root = Path(__file__).resolve().parents[1]
sdf = pd.read_parquet(root / "data" / "processed" / "standardized_2d_structures.parquet")
rng = np.random.RandomState(42)

def make_atoms(row, dummy=True):
    cell = np.array(json.loads(row["lattice"])); pos = np.array(json.loads(row["positions"])); sp = json.loads(row["species"])
    if dummy: sp = ["X"]*len(sp)
    return Atoms(symbols=sp, positions=pos, cell=cell, pbc=[True, True, False])

def soap_global(atoms_list, rcut=6.0):
    soap = SOAP(r_cut=rcut, n_max=6, l_max=6, sigma=1.0, periodic=True, average="off", species=["X"], dtype="float64")
    X = soap.create(atoms_list, n_jobs=-1)
    v = np.array([x.mean(axis=0) for x in X])
    v = v / np.linalg.norm(v, axis=1, keepdims=True)  # L2 normalize
    return v

def kernel(v1, v2):
    return np.dot(v1, v2)

rows = []
# Test A: atom permutation (50 材料)
sel = rng.choice(len(sdf), 50, replace=False)
for i in sel:
    a = make_atoms(sdf.iloc[i])
    perm = rng.permutation(len(a))
    b = a[perm]
    v = soap_global([a, b])
    rows.append({"jid": sdf.iloc[i]["jid"], "test": "A_permutation", "similarity": round(kernel(v[0], v[1]), 10)})

# Test B: translation (50 材料)
for i in sel:
    a = make_atoms(sdf.iloc[i])
    b = a.copy()
    b.positions = b.positions + np.array([1.0, 1.0, 3.0])  # 面内+非周期平移
    v = soap_global([a, b])
    rows.append({"jid": sdf.iloc[i]["jid"], "test": "B_translation", "similarity": round(kernel(v[0], v[1]), 10)})

# Test C: vacuum invariance (50 材料, 15/20/25/30 A)
for i in sel:
    a = make_atoms(sdf.iloc[i])
    cell = a.cell.copy(); pos = a.positions.copy()
    versions = []
    for L in [15.0, 20.0, 25.0, 30.0]:
        c = cell.copy(); c[2] = c[2] / np.linalg.norm(c[2]) * L  # 缩放 vacuum 轴到 L
        p = pos.copy()
        # 沿 vacuum 轴 recenter
        frac = p @ np.linalg.inv(c); frac[:,2] = (frac[:,2] - frac[:,2].mean() + 0.5) % 1.0
        p = frac @ c
        versions.append(Atoms(symbols=a.get_chemical_symbols(), positions=p, cell=c, pbc=[True,True,False]))
    v = soap_global(versions)
    # 所有 pairwise similarity 最小值
    sims = [kernel(v[i], v[j]) for i in range(4) for j in range(i+1, 4)]
    rows.append({"jid": sdf.iloc[i]["jid"], "test": "C_vacuum", "similarity": round(min(sims), 10)})

# Test D: supercell invariance (20 材料, 2x2)
selD = rng.choice(len(sdf), 20, replace=False)
for i in selD:
    a = make_atoms(sdf.iloc[i])
    b = a * (2, 2, 1)  # 面内 2x2 supercell
    v = soap_global([a, b])
    rows.append({"jid": sdf.iloc[i]["jid"], "test": "D_supercell", "similarity": round(kernel(v[0], v[1]), 10)})

df = pd.DataFrame(rows)
df.to_csv(root / "data" / "audit" / "structure_descriptor_invariance.csv", index=False)
print(df.groupby("test")["similarity"].agg(["min","median","mean"]).to_string())
print(f"\nwrote structure_descriptor_invariance.csv ({len(df)} rows)")
