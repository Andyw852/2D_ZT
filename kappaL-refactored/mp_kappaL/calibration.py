"""Step 7：跨视图度量的生成模型敏感性曲线。

问题（实测）：Elastic 与 kL_clarke 数学上是 R²=0.98 的确定性关系，但跨视图 Spearman
只读到 +0.65。既然确定性关系的读数上限是 0.65，Structure 的 0.36 / Eg 的 0.09 就
无法解释 —— 我们不知道这把尺子的刻度。高维↔高维与高维↔1维不是同一把尺子（测度集中）。

做法：用真实视图（保留真实距离分布），构造合成目标
    y = c + noise,  c = 真实描述符块的随机线性组合
调节噪声使 y 与 c 的真实 R² 落在 {0.0,0.1,...,1.0}，对每个水平跑同一套跨视图流程，
画出「模拟 R² → 度量读数」曲线。反演结果只在相同视图和相同生成机制下成立，
不能当作真实物理关系的通用 R² 估计。

view_case:
  "structure_1d" / "structure_hd" / "elastic_1d" / "eg_1d"。

特别地，Eg↔κL 是 1D↔1D，旧代码却套用了 Structure↔1D 的曲线；由此得到的
“等效 R²≈0.48”没有统计学含义。
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform

ROOT = Path(__file__).resolve().parents[1]
for p in (str(ROOT), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)
import config
from mp_kappaL.data_utils import load_aligned
from mp_kappaL import metrics
from graph_utils import hellinger_distance, soap_distance


def _zscore(M):
    return (M - M.mean(axis=0)) / (M.std(axis=0) + 1e-8)


def _sampled_spearman(Da, Db, rng, n_samp=80_000):
    n = Da.shape[0]
    iu = np.triu_indices(n, k=1)
    nsamp = min(n_samp, len(iu[0]))
    sel = rng.choice(len(iu[0]), size=nsamp, replace=False)
    ip, jp = iu[0][sel], iu[1][sel]
    return stats.spearmanr(Da[ip, jp], Db[ip, jp]).statistic


def _make_noisy(c, q, rng):
    if q <= 0.0:
        return rng.standard_normal(len(c))
    if q >= 1.0:
        return c.copy()
    return c + rng.standard_normal(len(c)) * np.sqrt((1.0 - q) / q)


def _overlap_between(Da, Db, k, rng):
    nnA = metrics.knn_neighbor_matrix(Da, k, tiebreak_seed=config.SEED)
    nnB = metrics.knn_neighbor_matrix(Db, k, tiebreak_seed=config.SEED)
    n = Da.shape[0]
    ov = np.array([np.intersect1d(nnA[i], nnB[i]).size for i in range(n)], dtype=float) / k
    return float(ov.mean())


def build_curves():
    meta = pd.read_parquet(config.PROC_DIR / "views_meta.parquet")
    elec = pd.read_parquet(config.PROC_DIR / "electronic_by_mpid.parquet")
    meta = meta.merge(elec[["material_id", "band_gap"]], on="material_id", how="left")
    meta, feats = load_aligned(
        {"soap_geo": config.PROC_DIR / "soap_geo.npy",
         "comp_frac": config.PROC_DIR / "comp_frac.npy"}, meta)
    soap_geo = feats["soap_geo"].astype(np.float32)
    comp_frac = feats["comp_frac"].astype(np.float32)

    rng = np.random.RandomState(config.SEED)
    n_cal = min(config.N_CALIB_SAMPLE, len(meta))
    idx = rng.choice(len(meta), size=n_cal, replace=False)

    # 真实视图距离（保留真实分布）
    d_geo = soap_distance(soap_geo[idx]).astype(np.float32); d_geo /= d_geo.max()
    d_comp = hellinger_distance(comp_frac[idx]).astype(np.float32); d_comp /= d_comp.max()
    d_struct = (0.5 * d_geo + 0.5 * d_comp).astype(np.float32)
    E = meta[["bulk_vrh", "shear_vrh", "debye"]].values[idx].astype(np.float32)
    d_elast = squareform(pdist(_zscore(E))).astype(np.float32)

    eg_all = meta["band_gap"].dropna().index.to_numpy()
    n_eg = min(config.N_CALIB_SAMPLE, len(eg_all))
    idx_eg = rng.choice(eg_all, size=n_eg, replace=False)
    eg = meta.loc[idx_eg, "band_gap"].to_numpy(dtype=float)
    eg = (eg - eg.mean()) / (eg.std() + 1e-12)
    d_eg = np.abs(eg[:, None] - eg[None, :]).astype(np.float32)

    # 信号空间：完整结构特征（几何 SOAP + 组成），逐列标准化
    S_full = np.concatenate([_zscore(soap_geo[idx]), _zscore(comp_frac[idx])], axis=1).astype(np.float32)

    q_grid = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    rows = []
    for q in q_grid:
        for rep in range(config.N_CALIB_REPEAT):
            # 1d：随机方向
            w = rng.standard_normal(S_full.shape[1]); w /= np.linalg.norm(w)
            c = (S_full @ w).astype(np.float64)
            c = (c - c.mean()) / (c.std() + 1e-12)
            y = _make_noisy(c, q, rng)
            Dy = np.abs(y[:, None] - y[None, :]).astype(np.float32)
            rows.append({"q": q, "view_case": "structure_1d", "rep": rep,
                         "spearman": _sampled_spearman(d_struct, Dy, rng),
                         "overlap": _overlap_between(d_struct, Dy, config.K, rng)})

            # hd：5 维随机投影
            P = rng.standard_normal((5, S_full.shape[1])); P /= np.linalg.norm(P, axis=1, keepdims=True)
            C = (S_full @ P.T).astype(np.float64)
            C = (C - C.mean(0)) / (C.std(0) + 1e-12)
            if q <= 0.0:
                Y = rng.standard_normal(C.shape)
            elif q < 1.0:
                Y = C + rng.standard_normal(C.shape) * np.sqrt((1 - q) / q)
            else:
                Y = C.copy()
            Dyh = squareform(pdist(Y)).astype(np.float32)
            rows.append({"q": q, "view_case": "structure_hd", "rep": rep,
                         "spearman": _sampled_spearman(d_struct, Dyh, rng),
                         "overlap": _overlap_between(d_struct, Dyh, config.K, rng)})

            # elastic 正对照：随机方向在弹性空间
            we = rng.standard_normal(E.shape[1]); we /= np.linalg.norm(we)
            ce = (E @ we).astype(np.float64)
            ce = (ce - ce.mean()) / (ce.std() + 1e-12)
            ye = _make_noisy(ce, q, rng)
            Dye = np.abs(ye[:, None] - ye[None, :]).astype(np.float32)
            rows.append({"q": q, "view_case": "elastic_1d", "rep": rep,
                         "spearman": _sampled_spearman(d_elast, Dye, rng),
                         "overlap": _overlap_between(d_elast, Dye, config.K, rng)})

            # Eg 是 1D 且含大量 ties，必须用它自己的生成曲线。
            yeg = _make_noisy(eg, q, rng)
            Dyeg = np.abs(yeg[:, None] - yeg[None, :]).astype(np.float32)
            rows.append({"q": q, "view_case": "eg_1d", "rep": rep,
                         "spearman": _sampled_spearman(d_eg, Dyeg, rng),
                         "overlap": _overlap_between(d_eg, Dyeg, config.K, rng)})

    df = pd.DataFrame(rows)
    df.to_csv(config.PROC_DIR / "calibration_curve.csv", index=False)
    return df


def calibrate(rho_observed, view_case="structure_1d", curve_df=None):
    """在指定生成模型内把读数反演为“模拟等效 R²”。

    这不是通用物理 R²；调用者必须使用与实测视图相同的 ``view_case``。
    """
    if curve_df is None:
        curve_df = pd.read_csv(config.PROC_DIR / "calibration_curve.csv")
    case_col = "view_case" if "view_case" in curve_df.columns else "dim_case"
    sub = curve_df[curve_df[case_col] == view_case]
    if sub.empty:
        raise ValueError(f"校准曲线中没有 view_case={view_case!r}")
    summ = sub.groupby("q")["spearman"].agg(["mean", "std"]).reset_index()
    summ["lo"] = summ["mean"] - 1.96 * summ["std"]
    summ["hi"] = summ["mean"] + 1.96 * summ["std"]
    qs = summ["q"].values
    # Monte Carlo 波动会造成局部非单调；反演前用 isotonic regression 强制单调。
    from sklearn.isotonic import IsotonicRegression
    iso = IsotonicRegression(increasing=True, out_of_bounds="clip")
    mean_curve = iso.fit_transform(qs, summ["mean"].values)
    lo_curve = iso.fit_transform(qs, summ["lo"].values)
    hi_curve = iso.fit_transform(qs, summ["hi"].values)

    def inv(curve):
        if rho_observed >= curve[-1]:
            return 1.0
        if rho_observed <= curve[0]:
            return float(qs[0])
        return float(np.interp(rho_observed, curve, qs))

    q_hat = inv(mean_curve)
    q_lo = inv(hi_curve)   # 读数相同，若曲线整体更高则等效 R² 更低
    q_hi = inv(lo_curve)
    return (max(0.0, q_lo), q_hat, min(1.0, q_hi))


def plot_curve(df):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cases = ["structure_1d", "structure_hd", "elastic_1d", "eg_1d"]
    fig, ax = plt.subplots(2, 2, figsize=(11, 8.2))
    for i, dc in enumerate(cases):
        sub = df[df["view_case"] == dc]
        summ = sub.groupby("q")["spearman"].agg(["mean", "std"]).reset_index()
        a = ax.flat[i]
        a.errorbar(summ["q"], summ["mean"], yerr=1.96 * summ["std"], marker="o", capsize=3)
        a.set_xlabel("simulated R² parameter")
        a.set_ylabel("cross-view distance Spearman")
        a.set_title(f"generator sensitivity ({dc})")
        a.set_ylim(-0.05, 1.05)
        a.grid(alpha=0.3)
    plt.tight_layout(); plt.savefig(config.FIG_DIR / "calibration.png", dpi=130)
    print("saved figures/calibration.png")


def main():
    df = build_curves()
    plot_curve(df)
    print("=== 各 view_case 在模拟 R²=1 时的读数（仅限该生成模型）===")
    for dc in ["structure_1d", "structure_hd", "elastic_1d", "eg_1d"]:
        sub = df[df["view_case"] == dc]
        ceiling = sub[sub["q"] == 1.0]["spearman"].mean()
        print(f"  {dc:10s} ceiling(Spearman @ R²=1) = {ceiling:+.3f}")
    print("=== 仅作生成模型内敏感性对照；不可解释为物理 R² ===")
    observed = pd.read_csv(config.PROC_DIR / "view_distance_corr.csv").set_index("pair")["spearman"]
    comparisons = [
        (float(observed["Structure vs kL_clarke"]), "Structure↔kL_clarke", "structure_1d"),
        (float(observed["Eg vs kL_clarke"]), "Eg↔kL_clarke", "eg_1d"),
        (float(observed["Elastic vs kL_clarke"]), "Elastic↔kL_clarke", "elastic_1d"),
    ]
    for rho, label, dc in comparisons:
        lo, h, hi = calibrate(rho, dc, df)
        print(f"  {label:28s} rho={rho:+.3f} -> 等效 R²_hat={h:.2f} [{lo:.2f}, {hi:.2f}]")


if __name__ == "__main__":
    main()
