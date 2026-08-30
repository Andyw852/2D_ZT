"""快速行序不变性审计（固定 material-id 子样本）。"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import RobustScaler

ROOT = Path(__file__).resolve().parents[1]
for p in (str(ROOT), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

import config
from graph_utils import hellinger_distance, soap_distance
from mp_kappaL import metrics
from mp_kappaL.data_utils import load_aligned


def _structure_distance(soap, comp):
    dg = soap_distance(soap).astype(np.float32)
    dc = hellinger_distance(comp).astype(np.float32)
    dg = np.clip(dg / np.quantile(dg[dg > 0], 0.95), 0, 1)
    dc = np.clip(dc / np.quantile(dc[dc > 0], 0.95), 0, 1)
    return (0.5 * dg + 0.5 * dc).astype(np.float32)


def _overlap(DA, DB, ids):
    a = metrics.knn_neighbor_matrix(DA, config.K, ids=ids)
    b = metrics.knn_neighbor_matrix(DB, config.K, ids=ids)
    return float(np.mean([
        np.intersect1d(a[i], b[i]).size / config.K for i in range(len(ids))
    ]))


def main():
    meta = pd.read_parquet(config.PROC_DIR / "views_meta.parquet")
    elec = pd.read_parquet(config.PROC_DIR / "electronic_by_mpid.parquet")
    meta = meta.merge(elec[["material_id", "band_gap"]], on="material_id", how="left")
    meta, feats = load_aligned(
        {"soap": config.PROC_DIR / "soap_geo.npy",
         "comp": config.PROC_DIR / "comp_frac.npy"}, meta)
    rng = np.random.RandomState(config.SEED)
    rows = []
    for A in ["Structure", "Elastic", "Eg"]:
        valid = np.ones(len(meta), dtype=bool)
        if A == "Eg":
            valid &= meta["band_gap"].notna().to_numpy()
        pool = np.flatnonzero(valid)
        take = rng.choice(pool, min(config.ROW_ORDER_AUDIT_N, len(pool)), replace=False)
        take.sort()
        ids = meta.loc[take, "material_id"].to_numpy()
        DK = np.abs(
            np.log10(meta.loc[take, "clarke"].to_numpy(float))[:, None]
            - np.log10(meta.loc[take, "clarke"].to_numpy(float))[None, :]
        ).astype(np.float32)
        if A == "Structure":
            DA = _structure_distance(feats["soap"][take], feats["comp"][take])
        elif A == "Elastic":
            X = meta.loc[take, ["bulk_vrh", "shear_vrh", "debye"]].to_numpy(float)
            DA = squareform(pdist(RobustScaler().fit_transform(X))).astype(np.float32)
        else:
            x = meta.loc[take, "band_gap"].to_numpy(float)
            DA = np.abs(x[:, None] - x[None, :]).astype(np.float32)
        base = _overlap(DA, DK, ids)
        for shuffle in range(config.N_ORDER_SHUFFLES):
            perm = rng.permutation(len(ids))
            got = _overlap(DA[np.ix_(perm, perm)], DK[np.ix_(perm, perm)], ids[perm])
            rows.append({"shuffle": shuffle, "pair": f"{A} vs kL_clarke",
                         "N_audit": len(ids), "overlap": got, "base": base,
                         "relative_deviation": abs(got - base) / (base + 1e-12)})
    detail = pd.DataFrame(rows)
    summary = detail.groupby(["pair", "N_audit"], as_index=False).agg(
        mean=("overlap", "mean"), std=("overlap", "std"),
        base=("base", "first"), max_dev=("relative_deviation", "max"))
    summary["cv"] = summary["std"] / (summary["mean"] + 1e-12)
    summary.to_csv(config.PROC_DIR / "row_order_robustness.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
