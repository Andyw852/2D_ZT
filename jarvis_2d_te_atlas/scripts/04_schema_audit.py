"""Phase E - Step 5: Schema Audit。

读取全部 1103 条 dft_2d 记录的 attributes 字段，统计每个字段：
dtype、N_total、N_nonmissing、N_missing、coverage、N_unique、sample_value。

缺失值定义（严格）：
- None / NaN
- -99999（JARVIS 官方 "not available" 哨兵值，int 或 float）
- "" / "na" / "NA" / "N/A" / "None" / "null" / "not available"
- 空 list []
注意：0 不是缺失值（band_gap=0 可能代表金属）。
"""
import json
import math
from collections import OrderedDict
from pathlib import Path

import pandas as pd

root = Path(__file__).resolve().parents[1]
raw = root / "data" / "raw" / "jarvis" / "dft_2d_snapshot.json"
audit_dir = root / "data" / "audit"
audit_dir.mkdir(exist_ok=True)

records = json.loads(raw.read_text(encoding="utf-8"))
print("records:", len(records))

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
    if isinstance(v, (list, tuple)):
        if len(v) == 0:
            return True
    if isinstance(v, dict):
        if len(v) == 0:
            return True
    return False

def dtype_name(v):
    if v is None:
        return "null"
    if isinstance(v, bool):
        return "bool"
    if isinstance(v, int):
        return "int"
    if isinstance(v, float):
        return "float"
    if isinstance(v, str):
        return "str"
    if isinstance(v, list):
        return "list"
    if isinstance(v, dict):
        return "dict"
    return type(v).__name__

# 收集所有字段
all_attrs = [r.get("attributes", {}) for r in records]
all_keys = OrderedDict()
for a in all_attrs:
    for k in a.keys():
        all_keys[k] = all_keys.get(k, 0) + 1

rows = []
for key in all_keys.keys():
    vals = [a.get(key) for a in all_attrs]
    n_total = len(vals)
    n_missing = sum(1 for v in vals if is_missing(v))
    n_nonmissing = n_total - n_missing
    coverage = n_nonmissing / n_total if n_total else 0.0
    # 类型：取第一个非缺失值
    sample = None
    dtype = "unknown"
    for v in vals:
        if not is_missing(v):
            sample = v
            dtype = dtype_name(v)
            break
    # 唯一值数量（非缺失）
    uniq = set()
    for v in vals:
        if not is_missing(v):
            try:
                uniq.add(repr(v))
            except Exception:
                uniq.add(str(v))
    n_unique = len(uniq)
    rows.append({
        "field": key,
        "dtype": dtype,
        "N_total": n_total,
        "N_nonmissing": n_nonmissing,
        "N_missing": n_missing,
        "coverage": round(coverage, 6),
        "N_unique": n_unique,
        "sample": repr(sample)[:120],
    })

df = pd.DataFrame(rows)
df = df.sort_values(["N_nonmissing", "field"], ascending=[False, True])
csv_path = audit_dir / "schema.csv"
df.to_csv(csv_path, index=False)
print(f"Wrote {csv_path} with {len(df)} fields")
print()
print(df.to_string(index=False))
