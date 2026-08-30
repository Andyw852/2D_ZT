"""Phase L (L2-L4): 二维结构 vacuum axis 识别与审计。"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

root = Path(__file__).resolve().parents[1]
recs = json.loads((root / "data" / "raw" / "jarvis" / "dft_2d_snapshot.json").read_text(encoding="utf-8"))

def audit(attrs):
    cell = np.array(attrs["lattice_vectors"], dtype=float)      # 3x3, rows = lattice vectors
    cart = np.array(attrs["cartesian_site_positions"], dtype=float)
    species = attrs["species_at_sites"]
    # fractional coords: cart = frac @ cell  => frac = cart @ inv(cell)
    frac = cart @ np.linalg.inv(cell)
    frac = frac % 1.0
    lens = np.linalg.norm(cell, axis=1)
    gaps = []
    for a in range(3):
        f = np.sort(frac[:, a])
        d = np.diff(np.concatenate([[f[-1] - 1.0], f]))  # 含周期边界: f[0]-f[-1]+1, f[1]-f[0], ...
        lg = d.max() * lens[a]
        gaps.append(lg)
    gaps = np.array(gaps)
    order = np.argsort(-gaps)  # 降序
    gap1, gap2, gap3 = gaps[order[0]], gaps[order[1]], gaps[order[2]]
    conf = gap1 / (gap2 + 1e-9)
    return {
        "jid": attrs["_jarvis_jid"],
        "formula": attrs["chemical_formula_reduced"],
        "vacuum_axis": int(order[0]),
        "vacuum_gap": round(gap1, 4),
        "second_largest_gap": round(gap2, 4),
        "third_gap": round(gap3, 4),
        "vacuum_confidence": round(conf, 4),
        "status": "OK" if conf >= 1.5 else "AMBIGUOUS",
    }

rows = [audit(r["attributes"]) for r in recs]
df = pd.DataFrame(rows).sort_values("jid").reset_index(drop=True)
df.to_csv(root / "data" / "audit" / "vacuum_axis_audit.csv", index=False)

print(f"total structures: {len(df)}")
print("vacuum_axis distribution:", df["vacuum_axis"].value_counts().to_dict())
print("vacuum_confidence: median=%.3f p10=%.3f p90=%.3f min=%.3f" % (
    df["vacuum_confidence"].median(), df["vacuum_confidence"].quantile(0.1),
    df["vacuum_confidence"].quantile(0.9), df["vacuum_confidence"].min()))
print("status counts:", df["status"].value_counts().to_dict())
amb = df[df["status"] == "AMBIGUOUS"]
print(f"ambiguous (conf<1.5): {len(amb)}")
print(amb[["jid","formula","vacuum_axis","vacuum_gap","second_largest_gap","vacuum_confidence"]].head(20).to_string(index=False))
