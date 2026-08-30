"""数据工具：物理合理性清洗 + 按 material_id 的严格对齐。

- clean_records(): 一次性清洗（Step 2），输出被剔除样本清单，可审计。
- load_aligned(): 读入特征矩阵 + row_index，与 meta 做 inner join 重排（Step 3），
  禁止裸的按行号对齐。
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import config


def clean_records(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """按 config.PHYSICAL_RANGES 做一次物理合理性清洗。

    返回 (df_clean, rejected)。rejected 每行含 material_id、触发的 rule（字段）、
    原始 field_value 与 bounds，写入 processed/rejected_samples.csv 以便审计。
    清洗后对每个字段做越界断言（Step 2 要求 6）。
    """
    df = df.copy()
    # 独立计算每个字段的越界掩码（不顺序删除，保证每条规则的命中数是独立可审计的）
    bad_mask = pd.Series(False, index=df.index)
    rejected_rows = []
    for col, (lo, hi) in config.PHYSICAL_RANGES.items():
        if col not in df.columns:
            continue
        s = df[col]
        bad = s.notna() & ((s < lo) | (s > hi))
        if bad.any():
            sub = df.loc[bad, ["material_id", col]]
            for mid, val in zip(sub["material_id"], sub[col]):
                rejected_rows.append({
                    "material_id": mid,
                    "rule": col,
                    "field_value": float(val) if pd.notna(val) else np.nan,
                    "bounds": f"[{lo}, {hi}]",
                })
        bad_mask = bad_mask | bad
    df = df.loc[~bad_mask].reset_index(drop=True)

    # 断言：清洗后所有字段都落在区间内，否则 raise（Step 2 要求 6）
    for col, (lo, hi) in config.PHYSICAL_RANGES.items():
        if col in df.columns:
            s = df[col]
            assert ((s >= lo) & (s <= hi)).all(), f"清洗后字段 {col} 仍存在越界值"

    rejected = pd.DataFrame(rejected_rows)
    return df, rejected


def load_aligned(feature_paths: dict[str, str], meta: pd.DataFrame,
                 id_col: str = "material_id",
                 row_index_path=None) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    """加载若干特征矩阵，按 material_id 与 meta 做 inner join 并对齐。

    feature_paths: {name: npy_path}。所有矩阵共享同一 processed/row_index.npy
    （保存顺序与特征行一一对应）。函数内部对每一步做 assert（Step 3）。
    返回 (meta_aligned, {name: matrix_aligned})，矩阵行序与 meta_aligned 完全一致。
    """
    row_index_path = row_index_path or (config.PROC_DIR / "row_index.npy")
    assert Path(row_index_path).exists(), f"缺少 {row_index_path}（特征矩阵必须伴随显式 material_id 索引）"
    row_ids = np.load(row_index_path, allow_pickle=True).astype(str)

    mats: dict[str, np.ndarray] = {}
    for name, path in feature_paths.items():
        X = np.load(path).astype(np.float64)
        assert X.shape[0] == len(row_ids), (
            f"特征矩阵 {name} 行数 {X.shape[0]} != row_index 长度 {len(row_ids)}，"
            "禁止裸按行号对齐（Step 3）"
        )
        mats[name] = X

    meta = meta.copy()
    assert meta[id_col].is_unique, f"{id_col} 存在重复键，merge 会静默错位（bug C）"
    assert len(set(row_ids)) == len(row_ids), "row_index.npy 中存在重复 material_id"

    # 只保留 meta 与特征索引都有的 material_id，并按 meta 的 material_id 顺序重排
    idx_df = pd.DataFrame({"__idx": np.arange(len(row_ids)), id_col: row_ids})
    meta_join = meta.merge(idx_df, on=id_col, how="inner")
    assert meta_join["__idx"].notna().all(), "对齐后存在未命中的 material_id"
    assert len(meta_join) == len(meta_join[id_col].drop_duplicates()), (
        "inner join 后出现重复行，特征索引或 meta 键不唯一"
    )
    order = meta_join["__idx"].values.astype(int)
    meta_aligned = meta_join.drop(columns=["__idx"]).reset_index(drop=True)
    for name in mats:
        mats[name] = mats[name][order]
    # 最终断言：每行一一对应
    assert len(meta_aligned) == len(mats[list(mats)[0]]), "对齐后长度不一致"
    return meta_aligned, mats
