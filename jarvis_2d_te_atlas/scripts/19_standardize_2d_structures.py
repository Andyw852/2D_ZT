"""Phase L (L5-L9): 标准化 2D 结构 + composition fraction。"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

root = Path(__file__).resolve().parents[1]
recs = json.loads((root / "data" / "raw" / "jarvis" / "dft_2d_snapshot.json").read_text(encoding="utf-8"))
audit = pd.read_csv(root / "data" / "audit" / "vacuum_axis_audit.csv").set_index("jid")
vac_axis = audit["vacuum_axis"].to_dict()

def recenter(cell, cart, vac=2):
    frac = cart @ np.linalg.inv(cell)
    frac = frac % 1.0
    # 沿 vacuum 轴把原子质心移到 0.5
    shift = 0.5 - frac[:, vac].mean()
    frac[:, vac] = (frac[:, vac] + shift) % 1.0
    return frac @ cell

rows = []
for r in recs:
    a = r["attributes"]
    jid = a["_jarvis_jid"]
    cell = np.array(a["lattice_vectors"], dtype=float)
    cart = np.array(a["cartesian_site_positions"], dtype=float)
    species = a["species_at_sites"]
    v = int(vac_axis[jid])
    if v != 2:
        # 重排: 把 v 轴换到第 2 轴
        perm = [i for i in range(3) if i != v] + [v]
        cell = cell[perm]
        # cart 坐标列也要重排
        cart = cart[:, perm]
    cart = recenter(cell, cart, vac=2)
    rows.append({
        "jid": jid, "formula": a["chemical_formula_reduced"],
        "lattice": json.dumps(cell.tolist()),
        "positions": json.dumps(cart.tolist()),
        "species": json.dumps(species),
        "nsites": len(species),
    })
sdf = pd.DataFrame(rows).sort_values("jid").reset_index(drop=True)
sdf.to_parquet(root / "data" / "processed" / "standardized_2d_structures.parquet", index=False)
print(f"standardized structures: {sdf.shape}")

# composition fraction (L8)
all_species = sorted({sp for s in sdf["species"] for sp in json.loads(s)})
print(f"distinct elements in database: {len(all_species)} -> {all_species}")
el_idx = {e: i for i, e in enumerate(all_species)}
frac_rows = []
for _, row in sdf.iterrows():
    sp = json.loads(row["species"])
    frac = np.zeros(len(all_species))
    for s in sp:
        frac[el_idx[s]] += 1.0 / len(sp)
    frac_rows.append({"jid": row["jid"], "fraction": json.dumps(frac.tolist()), "n_species": len(sp)})
fdf = pd.DataFrame(frac_rows)
fdf.to_parquet(root / "features" / "structure" / "composition_fraction.parquet", index=False)
print(f"composition fraction: {fdf.shape} (element dim = {len(all_species)})")

# 保存元素列表
meta = {"elements": all_species, "n_elements": len(all_species), "distance": "Hellinger"}
import json as J
(root / "features" / "structure" / "composition_element_list.json").write_text(J.dumps(meta), encoding="utf-8")
print("wrote composition_fraction.parquet + element list")
