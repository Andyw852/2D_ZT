"""MP 数据：跨视图近邻重叠 + 成对距离相关（重构版，Step 3/4/9）。

Step 3：load_aligned 按 material_id 严格对齐，禁止裸按行号。
Step 4：kNN 统一带 material 身份破平局（tiebreak）；并列诊断 + 行序稳健性。
Step 9：报置换 z + bootstrap 95% CI（不报顶下限的 p）；k 敏感性 + BH 校正。
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

ROOT = Path(__file__).resolve().parents[1]
for p in (str(ROOT), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)
import config
from mp_kappaL.data_utils import load_aligned
from mp_kappaL import metrics
from graph_utils import hellinger_distance, soap_distance


def build_view(D, ids):
    return {"D": np.asarray(D, dtype=np.float32), "ids": np.asarray(ids)}


def col_distance(meta, cols, log=False, scale=True, filt=None):
    """从 meta 列构造 (n,n) 距离矩阵 + ids。filt 为布尔掩码（如 Snyder 数值爆炸过滤）。"""
    sub = meta.copy()
    if filt is not None:
        sub = sub[filt]
    sub = sub.dropna(subset=cols).reset_index(drop=True)
    X = sub[cols].values.astype(float)
    if log:
        X = np.log10(np.clip(X, 1e-12, None))
    if scale and X.shape[1] > 1:
        from sklearn.preprocessing import RobustScaler
        X = RobustScaler().fit_transform(X)
    D = squareform(pdist(X)).astype(np.float32)
    return build_view(D, sub["material_id"].to_numpy())


def common_pair_distances(views, A, B, common=None):
    """在两个视图的共同材料上先取子矩阵，再进行任何 kNN/距离比较。

    旧实现先在各自全集建 kNN，之后丢弃不在交集内的邻居。覆盖率不同（尤其
    Eg 只有约 20% 材料）时，每行会不足 k 个邻居，却仍除以 k，系统性压低重叠。
    """
    if common is None:
        common = sorted(set(views[A]["ids"]) & set(views[B]["ids"]))
    common = np.asarray(common)
    pos_a = {mid: i for i, mid in enumerate(views[A]["ids"])}
    pos_b = {mid: i for i, mid in enumerate(views[B]["ids"])}
    ia = np.fromiter((pos_a[m] for m in common), dtype=int, count=len(common))
    ib = np.fromiter((pos_b[m] for m in common), dtype=int, count=len(common))
    DA = views[A]["D"][np.ix_(ia, ia)]
    DB = views[B]["D"][np.ix_(ib, ib)]
    return common, DA, DB


def common_pair_neighbors(views, A, B, k, common=None):
    common, DA, DB = common_pair_distances(views, A, B, common=common)
    NNA = metrics.knn_neighbor_matrix(DA, k, ids=common)
    NNB = metrics.knn_neighbor_matrix(DB, k, ids=common)
    return common, NNA, NNB, DA, DB


def main():
    meta = pd.read_parquet(config.PROC_DIR / "views_meta.parquet")
    elec = pd.read_parquet(config.PROC_DIR / "electronic_by_mpid.parquet")
    meta = meta.merge(elec, on="material_id", how="left")
    assert meta["material_id"].is_unique, "merge 后 material_id 不唯一（bug C）"

    meta, feats = load_aligned(
        {"soap_geo": config.PROC_DIR / "soap_geo.npy",
         "comp_frac": config.PROC_DIR / "comp_frac.npy"}, meta)
    soap_geo = feats["soap_geo"].astype(np.float32)
    comp_frac = feats["comp_frac"].astype(np.float32)
    print("N materials (cleaned):", len(meta))

    # ---- 距离矩阵 ----
    d_geo = soap_distance(soap_geo).astype(np.float32)
    d_comp = hellinger_distance(comp_frac).astype(np.float32)
    # robust scale + clip：max 会被单个离群点控制。
    d_geo_n = np.clip(d_geo / np.quantile(d_geo[d_geo > 0], 0.95), 0, 1)
    d_comp_n = np.clip(d_comp / np.quantile(d_comp[d_comp > 0], 0.95), 0, 1)
    d_struct = (0.5 * d_geo_n + 0.5 * d_comp_n).astype(np.float32)
    ids_all = meta["material_id"].to_numpy()

    snyder_filt = (meta["snyder_acoustic"] > 0) & (meta["snyder_acoustic"] < 1e4)

    views = {
        "Structure":      build_view(d_struct, ids_all),
        "Structure_geo":  build_view(d_geo_n, ids_all),
        "Structure_comp": build_view(d_comp_n, ids_all),
        "kL_clarke":      col_distance(meta, ["clarke"], log=True, scale=False),
        "kL_cahill":      col_distance(meta, ["cahill"], log=True, scale=False),
        "kL_snyder":      col_distance(meta, ["snyder_acoustic"], log=True, scale=False, filt=snyder_filt),
        "Elastic":        col_distance(meta, ["bulk_vrh", "shear_vrh", "debye"]),
        "Eg":             col_distance(meta, ["band_gap"], scale=False),
        "Electronic":     col_distance(meta, ["band_gap", "efermi"]),
    }

    pairs = [
        ("Structure", "kL_clarke"),
        ("Structure_geo", "kL_clarke"),
        ("Structure_comp", "kL_clarke"),
        ("Elastic", "kL_clarke"),
        ("Eg", "kL_clarke"),
        ("Electronic", "kL_clarke"),
        ("Structure", "kL_snyder"),
        ("Elastic", "kL_snyder"),
        ("Eg", "kL_snyder"),
        ("Structure", "Eg"),
        ("Electronic", "Elastic"),
        ("kL_clarke", "kL_cahill"),
    ]

    K = config.K
    rng = np.random.RandomState(config.SEED)

    # ---- 每视图并列诊断 ----
    degen_rows = []
    for name, v in views.items():
        dg = metrics.degeneracy_diagnostics(v["D"], K)
        dg["view"] = name
        dg["N"] = len(v["ids"])
        degen_rows.append(dg)
        print(f"  [{name:16s}] N={len(v['ids']):5d} zero_offdiag={dg['zero_offdiag_frac']:.4f} "
              f"k/k+1_tie={dg['k_k1_tie_frac']:.4f}")
    pd.DataFrame(degen_rows).to_csv(config.PROC_DIR / "knn_degeneracy.csv", index=False)

    # ---- 跨视图近邻重叠 + z + bootstrap CI ----
    print(f"=== cross-view neighbor overlap (k={K}) + permutation z + bootstrap CI ===")
    ov_rows = []
    for A, B in pairs:
        common = sorted(set(views[A]["ids"]) & set(views[B]["ids"]))
        if len(common) < 20:
            print(f"  {A} vs {B}: N={len(common)} <20 skip"); continue
        common, NNA, NNB, _, _ = common_pair_neighbors(views, A, B, K, common)
        r = metrics.crossview_overlap(NNA, NNB, K, rng, n_perm=config.N_PERM, n_boot=config.N_BOOTSTRAP)
        ov_rows.append({"pair": f"{A} vs {B}", **{k: r[k] for k in
                          ["n", "overlap", "null_mean", "z", "enrichment", "ci_lo", "ci_hi"]}})
        print(f"  {A:16s} vs {B:14s}: N={r['n']:5d} overlap={r['overlap']:.4f} "
              f"null={r['null_mean']:.4f} z={r['z']:7.1f} enrich={r['enrichment']:.2f}x "
              f"CI=[{r['ci_lo']:.4f},{r['ci_hi']:.4f}]")
    ov_df = pd.DataFrame(ov_rows)
    if len(ov_df):
        ov_df["enrich_q"] = metrics.benjamini_hochberg(ov_df["z"].abs().map(_z_to_p))
        ov_df.to_csv(config.PROC_DIR / "view_overlap.csv", index=False)

    # ---- 成对距离相关 + Mantel z + bootstrap CI ----
    print(f"=== cross-view distance Spearman + Mantel z + bootstrap CI ===")
    dr_rows = []
    for A, B in pairs:
        common = sorted(set(views[A]["ids"]) & set(views[B]["ids"]))
        if len(common) < 20: continue
        _, DA, DB = common_pair_distances(views, A, B, common)
        r = metrics.distance_spearman(DA, DB, rng)
        dr_rows.append({"pair": f"{A} vs {B}", **{k: r[k] for k in
                         ["n", "spearman", "z", "ci_lo", "ci_hi"]}})
        print(f"  {A:16s} vs {B:14s}: N={r['n']:5d} Spearman={r['spearman']:+.4f} "
              f"z={r['z']:6.1f} CI=[{r['ci_lo']:+.4f},{r['ci_hi']:+.4f}]")
    dr_df = pd.DataFrame(dr_rows)
    if len(dr_df):
        dr_df.to_csv(config.PROC_DIR / "view_distance_corr.csv", index=False)

    # ---- k 敏感性（Step 9）----
    print(f"=== k sensitivity (k in {config.K_SENSITIVITY}) ===")
    ksen_rows = []
    for kk in config.K_SENSITIVITY:
        for A, B in [("Structure", "kL_clarke"), ("Elastic", "kL_clarke"),
                     ("Eg", "kL_clarke"), ("Structure", "Eg")]:
            common = sorted(set(views[A]["ids"]) & set(views[B]["ids"]))
            _, NNA, NNB, _, _ = common_pair_neighbors(views, A, B, kk, common)
            r = metrics.crossview_overlap(NNA, NNB, kk, rng, n_perm=100, n_boot=200)
            ksen_rows.append({"k": kk, "pair": f"{A} vs {B}", "overlap": r["overlap"],
                              "z": r["z"], "enrichment": r["enrichment"]})
    pd.DataFrame(ksen_rows).to_csv(config.PROC_DIR / "k_sensitivity.csv", index=False)

    # ---- 行序稳健性（Step 4）----
    print(f"=== row-order robustness ({config.N_ORDER_SHUFFLES} shuffles) ===")
    headline = [("Structure", "kL_clarke"), ("Elastic", "kL_clarke"), ("Eg", "kL_clarke")]
    common_by_pair = {}
    for A, B in headline:
        ids = np.array(sorted(set(views[A]["ids"]) & set(views[B]["ids"])))
        if len(ids) > config.ROW_ORDER_AUDIT_N:
            audit_rng = np.random.RandomState(config.SEED)
            ids = np.sort(audit_rng.choice(ids, config.ROW_ORDER_AUDIT_N, replace=False))
        common_by_pair[(A, B)] = ids

    def headline_overlap(common_orders):
        res = {}
        for A, B in headline:
            _, NNA, NNB, _, _ = common_pair_neighbors(
                views, A, B, K, common_orders[(A, B)])
            res[(A, B)] = metrics.crossview_overlap(
                NNA, NNB, K, rng, n_perm=100, n_boot=100)["overlap"]
        return res

    base = headline_overlap(common_by_pair)
    rob_rows = []
    for s in range(config.N_ORDER_SHUFFLES):
        shuffled = {pair: ids[rng.permutation(len(ids))] for pair, ids in common_by_pair.items()}
        ovs = headline_overlap(shuffled)
        for (A, B), ov in ovs.items():
            rob_rows.append({"shuffle": s, "pair": f"{A} vs {B}",
                             "N_audit": len(common_by_pair[(A, B)]), "overlap": ov})
    rob_df = pd.DataFrame(rob_rows)
    if len(rob_df):
        summ = rob_df.groupby("pair")["overlap"].agg(["mean", "std"]).reset_index()
        summ["cv"] = summ["std"] / (summ["mean"] + 1e-12)
        summ["base"] = summ["pair"].map(lambda p: base[tuple(p.split(" vs "))])
        summ["max_dev"] = abs(summ["mean"] - summ["base"]) / (summ["base"] + 1e-12)
        summ.to_csv(config.PROC_DIR / "row_order_robustness.csv", index=False)
        print(summ.to_string(index=False))

    print("saved: knn_degeneracy.csv, view_overlap.csv, view_distance_corr.csv, "
          "k_sensitivity.csv, row_order_robustness.csv")


def _z_to_p(z):
    from scipy.stats import norm
    return 2.0 * norm.sf(abs(z))


if __name__ == "__main__":
    main()
