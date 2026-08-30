"""Step 5：按 material_id 重建电子性质视图（废弃化学式跨库匹配）。

v0 的 match_electronic.py 用 pymatgen reduced_formula 跨库匹配 MP↔JARVIS，导致
14.8% 的 MP 材料共享同一 canon（26 个 SiO2、21 个 Al2O3、18 个 C...）被赋予完全
相同的电子性质，人为重复抬高近邻重叠。

本脚本改用 MP 自己的 summary band_gap（raw/summary_bandgap.json，按 material_id
严格对应，零匹配噪声）——download_mp.py 早已下载却从未使用。
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import config


def main():
    sb = json.load(open(config.RAW_DIR / "summary_bandgap.json"))
    rows = []
    for mid, x in sb.items():
        rows.append({
            "material_id": x["material_id"],
            "band_gap": float(x.get("band_gap", np.nan)) if x.get("band_gap") is not None else np.nan,
            "is_metal": bool(x.get("is_metal", False)),
            "efermi": float(x.get("efermi", np.nan)) if x.get("efermi") is not None else np.nan,
            "is_gap_direct": x.get("is_gap_direct"),
            "formula_pretty": x.get("formula_pretty"),
        })
    elec = pd.DataFrame(rows).drop_duplicates(subset=["material_id"])
    assert elec["material_id"].is_unique, "material_id 不唯一"
    print("MP-native electronic by material_id:", len(elec))

    n_metal = int((elec["is_metal"] | (elec["band_gap"] == 0)).sum())
    nonmetal = elec[~elec["is_metal"] & (elec["band_gap"] > 0)]
    feat_cols = ["band_gap", "is_metal", "efermi"]
    dup_all = elec[feat_cols].duplicated().mean()
    dup_nonmetal = nonmetal[feat_cols].duplicated().mean()
    print(f"金属(带隙=0) 数量: {n_metal} / {len(elec)}")
    print(f"重复特征向量占比(全体): {dup_all:.4f}")
    print(f"重复特征向量占比(非金属): {dup_nonmetal:.4f}")

    # 对比旧化学式匹配（供 Step 5 报告新旧差异）
    old = pd.read_parquet(config.PROC_DIR / "electronic_jarvis.parquet")
    old_dup = old[["gap_opt", "gap_mbj", "m_elec", "m_hole"]].duplicated().mean()
    print("旧化学式匹配 N:", len(old), " gap_opt==0:", int((old['gap_opt'] == 0).sum()),
          " 重复特征向量占比:", round(old_dup, 4))

    elec.to_parquet(config.PROC_DIR / "electronic_by_mpid.parquet", index=False)
    print("saved processed/electronic_by_mpid.parquet")


if __name__ == "__main__":
    main()
