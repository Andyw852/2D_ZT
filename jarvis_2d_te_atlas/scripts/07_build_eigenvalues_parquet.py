"""Phase G5/G6: 构建 eigenvalues parquet + 验证 OPTIMADE 均值定义 + 复数检查。"""
import json
import re
import numpy as np
import pandas as pd
from pathlib import Path

root = Path(__file__).resolve().parents[1]
raw = json.loads((root / "data" / "processed" / "transport_eigenvalues_raw.json").read_text(encoding="utf-8"))
opt = json.loads((root / "data" / "raw" / "jarvis" / "dft_2d_snapshot.json").read_text(encoding="utf-8"))
opt_by_jid = {r["attributes"]["_jarvis_jid"]: r["attributes"] for r in opt}

FIELD_MAP = {
    "nseeb": "_jarvis_n-Seebeck",
    "pseeb": "_jarvis_p-Seebeck",
    "ncond": "_jarvis_ncond",
    "pcond": "_jarvis_pcond",
    "npf": "_jarvis_n-powerfact",
    "ppf": "_jarvis_p-powerfact",
    "nkappa": "_jarvis_nkappa",
    "pkappa": "_jarvis_pkappa",
}
MASS_MAP = {
    "electron_mass_300K": "_jarvis_avg_elec_mass",
    "hole_mass_300K": "_jarvis_avg_hole_mass",
}

def to_float(s):
    s = s.strip()
    # 检查复数 / 非数值
    if re.search(r"[ijIJ]", s) and not re.search(r"[eE][+-]?\d", s):
        return None, "complex"
    try:
        return float(s), "ok"
    except ValueError:
        return None, "nonnumeric:" + s[:20]

# 1) 复数检查（第 12 节）
complex_count = 0
nonnumeric = 0
for jid, eig in raw.items():
    for f, vals in eig.items():
        for v in vals:
            fv, status = to_float(v)
            if status != "ok":
                if status == "complex":
                    complex_count += 1
                else:
                    nonnumeric += 1
print(f"complex/imaginary values: {complex_count}, non-numeric: {nonnumeric}")

# 2) 构建 eigenvalues parquet
rows = []
for jid, eig in raw.items():
    row = {"jid": jid}
    for f, vals in eig.items():
        nums = []
        for v in vals:
            fv, _ = to_float(v)
            nums.append(fv if fv is not None else np.nan)
        for i in range(3):
            row[f"{f}_{i}"] = nums[i] if i < len(nums) else np.nan
    rows.append(row)

df = pd.DataFrame(rows)
df = df.sort_values("jid").reset_index(drop=True)
out = root / "data" / "processed" / "transport_eigenvalues.parquet"
out.parent.mkdir(parents=True, exist_ok=True)
df.to_parquet(out, index=False)
print(f"wrote {out}: {df.shape}")

# 3) 验证 OPTIMADE scalar == mean(3 eigenvalues)（第 18 节）
rng = np.random.RandomState(0)
records = []
for f, optfield in FIELD_MAP.items():
    jids_with = [j for j in raw if f in raw[j] and optfield in opt_by_jid.get(j, {}) and opt_by_jid[j].get(optfield) not in (None, -99999, -99999.0)]
    if not jids_with:
        continue
    sample = rng.choice(jids_with, size=min(50, len(jids_with)), replace=False)
    for jid in sample:
        eigs = np.array([to_float(v)[0] for v in raw[jid][f]])
        mean_eig = eigs.mean()
        opt_val = float(opt_by_jid[jid][optfield])
        abs_err = abs(opt_val - mean_eig)
        rel_err = abs_err / max(abs(opt_val), 1e-12)
        records.append({"field": f, "jid": jid, "opt_scalar": opt_val, "mean_3eigs": mean_eig,
                        "abs_err": abs_err, "rel_err": rel_err, "exact": abs_err < 1e-6})

valdf = pd.DataFrame(records)
valcsv = root / "reports" / "optimade_mean_validation.csv"
valdf.to_csv(valcsv, index=False)
print(f"\n=== mean validation ({len(valdf)} samples) ===")
print(f"MAE = {valdf.abs_err.mean():.6g}, max abs err = {valdf.abs_err.max():.6g}")
print(f"exact matches (abs<1e-6): {valdf.exact.sum()} / {len(valdf)}")
print(f"mismatches: {(~valdf.exact).sum()}")
print(f"wrote {valcsv}")
