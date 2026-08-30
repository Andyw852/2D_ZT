"""Phase L (L10-L14): geometry-only SOAP descriptor + composition (already done)."""
import json
import numpy as np
import pandas as pd
import time
from pathlib import Path
from ase import Atoms
from dscribe.descriptors import SOAP

root = Path(__file__).resolve().parents[1]
sdf = pd.read_parquet(root / "data" / "processed" / "standardized_2d_structures.parquet")

def make_atoms(row, dummy=True):
    cell = np.array(json.loads(row["lattice"]))
    pos = np.array(json.loads(row["positions"]))
    sp = json.loads(row["species"])
    if dummy:
        sp = ["X"] * len(sp)
    return Atoms(symbols=sp, positions=pos, cell=cell, pbc=[True, True, False])

atoms_list = [make_atoms(sdf.iloc[i]) for i in range(len(sdf))]
jids = sdf["jid"].tolist()
print(f"built {len(atoms_list)} ASE atoms (geometry-only, dummy species X)")

# SOAP 参数（n_max=6, l_max=6, sigma=1.0; r_cut = 4/6/8）
params = {"n_max": 6, "l_max": 6, "sigma": 1.0, "periodic": True, "average": "off", "species": ["X"]}
rcuts = [4.0, 6.0, 8.0]

store = {"jid": jids}
for rcut in rcuts:
    soap = SOAP(r_cut=rcut, n_max=params["n_max"], l_max=params["l_max"], sigma=params["sigma"],
                periodic=params["periodic"], average=params["average"], species=params["species"], dtype="float64")
    t0 = time.time()
    X = soap.create(atoms_list, n_jobs=-1)  # list of (n_atoms, n_feat)
    means = np.array([x.mean(axis=0) for x in X])
    stds = np.array([x.std(axis=0) for x in X])
    print(f"rcut={rcut}: computed in {time.time()-t0:.1f}s, mean shape={means.shape}")
    for i in range(means.shape[1]):
        store[f"soap{int(rcut)}_mean_{i}"] = means[:, i]
    # std pooling 只保留作为 sensitivity（存到单独文件）
    # 先保存 mean
sdf_out = pd.DataFrame(store)
sdf_out.to_parquet(root / "features" / "structure" / "geometry_soap_v1.parquet", index=False)
print(f"wrote geometry_soap_v1.parquet: {sdf_out.shape}")

# metadata
meta = {
    "descriptor": "geometry-only SOAP (all atoms mapped to dummy species X)",
    "parameters": params,
    "r_cut_candidates": rcuts,
    "feature_dim_per_rcut": 147,
    "pooling": "mean over atoms",
    "normalization": "L2-normalize + SOAP kernel (applied in graph building)",
    "species_sensitive_soap": "SKIPPED_HIGH_DIMENSION (81 elements would give ~488k features)",
    "composition": "elemental fraction vector (81 elements), Hellinger distance",
}
(root / "features" / "structure" / "structure_feature_metadata.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
print("wrote structure_feature_metadata.json")
