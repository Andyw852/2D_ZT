"""Step 6：整理可用的晶格热导率目标，并显式记录映射质量。

数据源（本地可用）：
- experimental：starrydata2 实验 κ_L@300K（权威但小，作验证集）
- Snyder-300K-model：MP ``sound_velocity.snyder_acoustic``（解析模型代理；不是 AFLOW AGL）
- Clarke-min / Cahill-min：MP 非晶极限最小 κ_L（保留作参照下界，不再作主目标）

不同方法的 κ_L 绝不混用：每个源标注 method 与 temperature。
输出 kappa_L_targets.parquet（长表）+ 数据源交叉表 + 一致性散点图。
"""
import ast
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import config
from pymatgen.core import Composition

STARRY = config.EXTERNAL_DATA_DIR / "starrydata2"


def _parse_xy(v):
    v = (v or "").strip()
    if not v:
        return []
    try:
        d = json.loads(v)
    except Exception:
        try:
            d = ast.literal_eval(v) if v[:1] in "[(" else float(v)
        except Exception:
            return []
    if isinstance(d, (int, float)):
        return [float(d)]
    return [float(t) for t in d]


def _canon(comp):
    try:
        return Composition(comp).reduced_formula
    except Exception:
        return None


def load_starry_kappa_L():
    """starrydata2 实验 κ_L@300K，先归一化组成再聚合所有曲线。

    旧实现先按原始 composition 取中位数，再用 dict comprehension 归一化为
    reduced_formula；多个原始 composition 落到同一 canon 时会被静默覆盖。
    """
    curves = STARRY / "ThermoelectricMaterials_curves.csv"
    kappa_rows = defaultdict(list)
    with open(curves, encoding="utf-8", errors="ignore") as f:
        for row in csv.DictReader(f):
            py = (row.get("prop_y", "") or "").strip()
            comp = (row.get("composition", "") or "").strip()
            if not comp or py != "Lattice thermal conductivity":
                continue
            xs = _parse_xy(row.get("x")); ys = _parse_xy(row.get("y"))
            if len(xs) == len(ys) and xs:
                kappa_rows[comp].append((np.asarray(xs), np.asarray(ys)))
    canon_values = defaultdict(list)
    canon_compositions = defaultdict(set)
    for comp, cy in kappa_rows.items():
        vals = []
        for xs, ys in cy:
            m = (xs >= 250) & (xs <= 350)
            if m.any():
                val = float(ys[int(np.argmin(np.abs(xs - 300)))])
                if np.isfinite(val) and 0 < val <= 5000:
                    vals.append(val)
        if vals:
            canon = _canon(comp)
            if canon:
                canon_values[canon].extend(vals)
                canon_compositions[canon].add(comp)
    rows = []
    for canon, vals in canon_values.items():
        rows.append({
            "canon": canon,
            "kappa_L": float(np.median(vals)),
            "n_curves": len(vals),
            "n_raw_compositions": len(canon_compositions[canon]),
            "kappa_q25": float(np.quantile(vals, 0.25)),
            "kappa_q75": float(np.quantile(vals, 0.75)),
        })
    return pd.DataFrame(rows).sort_values("canon").reset_index(drop=True)


def main():
    meta = pd.read_parquet(config.PROC_DIR / "views_meta.parquet")
    meta["canon"] = meta["formula"].map(_canon)

    starry = load_starry_kappa_L()
    print("starrydata2 归一化 canon 数(有 κ_L@300K):", len(starry))

    canon_to_mids = meta.groupby("canon", dropna=True)["material_id"].agg(list).to_dict()

    # 长表：每个源一份
    rows = []
    for _, r in meta.iterrows():
        mid, canon = r["material_id"], r["canon"]
        rows.append({"material_id": mid, "formula": r["formula"],
                     "method": "Snyder-300K-model", "temperature": 300.0,
                     "kappa_L": float(r["snyder_acoustic"]),
                     "match_quality": "exact_material_id", "n_mp_matches": 1})
        rows.append({"material_id": mid, "formula": r["formula"], "method": "Clarke-min",
                     "temperature": 300.0, "kappa_L": float(r["clarke"]),
                     "match_quality": "exact_material_id", "n_mp_matches": 1})

    # 实验数据只有组成，没有结构/MP id。仅当该组成在清洗后的 MP 表中恰好对应
    # 一个 material_id 时才赋 id；多晶型不再复制同一个实验值制造伪样本。
    for _, r in starry.iterrows():
        mids = canon_to_mids.get(r["canon"], [])
        if len(mids) == 1:
            mid, quality = mids[0], "unique_formula_to_mpid"
        elif len(mids) > 1:
            mid, quality = pd.NA, "ambiguous_formula_polymorphs"
        else:
            mid, quality = pd.NA, "no_mp_match"
        rows.append({
            "material_id": mid,
            "formula": r["canon"],
            "method": "experimental",
            "temperature": 300.0,
            "kappa_L": float(r["kappa_L"]),
            "match_quality": quality,
            "n_mp_matches": len(mids),
            "n_curves": int(r["n_curves"]),
            "kappa_q25": float(r["kappa_q25"]),
            "kappa_q75": float(r["kappa_q75"]),
        })

    targets = pd.DataFrame(rows)
    targets.to_parquet(config.PROC_DIR / "kappa_L_targets.parquet", index=False)
    starry.to_csv(config.PROC_DIR / "experimental_formula_targets.csv", index=False)

    # 交叉表
    print("=== 数据源交叉表 ===")
    for mth in ["experimental", "Snyder-300K-model", "Clarke-min"]:
        sub = targets[targets["method"] == mth]
        q = sub["kappa_L"].quantile([0.01, 0.5, 0.99]) if len(sub) else pd.Series()
        print(f"  {mth:14s} N={len(sub):5d}  "
              f"p1={q.get(0.01, np.nan):.4g} p50={q.get(0.5, np.nan):.4g} p99={q.get(0.99, np.nan):.4g}")

    # 重叠 + 一致性散点
    exp_rows = targets[
        (targets["method"] == "experimental")
        & (targets["match_quality"] == "unique_formula_to_mpid")
    ].dropna(subset=["material_id"])
    exp = exp_rows.set_index("material_id")["kappa_L"]
    slack = targets[targets["method"] == "Snyder-300K-model"].set_index("material_id")["kappa_L"]
    common = exp.index.intersection(slack.index)
    n_amb = int((targets["match_quality"] == "ambiguous_formula_polymorphs").sum())
    print(f"experimental 唯一公式→MP 映射数: {len(common)}；多晶型歧义公式: {n_amb}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if len(common) > 3:
        from scipy import stats
        x = np.log10(slack.loc[common].values)
        y = np.log10(exp.loc[common].values)
        rho, _ = stats.spearmanr(x, y)
        fig, ax = plt.subplots(figsize=(5.2, 4.6))
        ax.scatter(x, y, s=12, alpha=0.6)
        lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
        ax.plot([lo, hi], [lo, hi], "k--", lw=1, label="y=x")
        ax.set_xlabel("log10 κ_L (Snyder 300 K model, MP)"); ax.set_ylabel("log10 κ_L (experimental, starrydata2)")
        ax.set_title(f"κ_L consistency (N={len(common)}, Spearman={rho:+.2f})")
        ax.legend()
        plt.tight_layout(); plt.savefig(config.FIG_DIR / "kappa_L_consistency.png", dpi=130)
        print("saved figures/kappa_L_consistency.png")
    print("saved processed/kappa_L_targets.parquet + experimental_formula_targets.csv")


if __name__ == "__main__":
    main()
