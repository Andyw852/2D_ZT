"""Step 8：增量信息方法 —— 控制已知描述符后，结构/几何还携带多少关于 κ_L 的新信息。

方法（都用按化学体系分组的 CV，避免同族材料泄漏）：
1. 块消融：Block-C(组成)/G(几何)/E(弹性)/X(电子) 逐块加入/剔除，报增量 R²（5-fold GroupKFold）。
2. 偏相关：每个描述符对 log κ_L 的偏 Spearman，控制 {density, avg_mass, debye}。
3. 残差互信息诊断：交叉拟合移除已知物理块后，几何 PC 与目标残差的 MI。

目标：
- log10(snyder_acoustic) —— Snyder 300 K 解析模型（大样本公式诊断，不是真值）
- log10(clarke) —— 参照（已知是弹性的确定性函数，用于对照）
- log10(κ_L experimental) —— 仅唯一公式→唯一 MP material_id 的实验验证（N=59）
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold, cross_val_predict, cross_val_score
from sklearn.feature_selection import mutual_info_regression

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import config
from mp_kappaL.data_utils import load_aligned
from pymatgen.core import Composition


def chemical_system(formula):
    try:
        return Composition(str(formula)).chemical_system
    except Exception:
        return None


def build_blocks(meta, feats):
    soap_geo = feats["soap_geo"]
    comp_frac = feats["comp_frac"]
    elems = np.load(config.PROC_DIR / "elem_basis.npy", allow_pickle=True)
    from pymatgen.core import Element
    mass = np.array([float(Element(e).atomic_mass.real) for e in elems])
    avg_mass = (comp_frac * mass).sum(axis=1)
    num_elem = (comp_frac > 1e-6).sum(axis=1)

    B_C = np.column_stack([comp_frac, avg_mass, num_elem])
    # 三个块严格互斥。density/nsites 放进已知物理控制块，避免把它们与 SOAP
    # 打包后将公式变量误称为“几何增量信息”。
    B_G = soap_geo
    B_E = meta[["bulk_vrh", "shear_vrh", "debye", "v_long", "v_trans",
              "density", "nsites"]].values.astype(float)

    blocks = {"C": B_C, "G": B_G, "E": B_E}
    # X 电子块只在有 band_gap 的材料上可用（按 id，N≈2700）
    if "band_gap" in meta.columns and meta["band_gap"].notna().sum() > 30:
        B_X = meta[["band_gap", "efermi"]].copy()
        B_X["is_metal"] = meta["is_metal"].astype(float)
        B_X = B_X.fillna(0.0).values.astype(float)
        blocks["X"] = B_X
    return blocks, avg_mass


def rf_cv(X, y, groups, n_splits=5):
    gkf = GroupKFold(n_splits=n_splits)
    # max_features 固定为全部特征：旧值 sqrt 会使“加入高维 SOAP”同时降低每次分裂
    # 看到弹性变量的概率，把模型超参数变化误当成负增量信息。
    rf = RandomForestRegressor(n_estimators=100, max_features=1.0,
                               min_samples_leaf=2, random_state=config.SEED, n_jobs=-1)
    scores = cross_val_score(rf, X, y, cv=gkf, groups=groups,
                             scoring="r2", error_score=np.nan)
    return float(np.nanmean(scores)), float(np.nanstd(scores))


def paired_geometry_increment(meta, blocks, target_col, groups):
    """同一 folds 上比较 C+E 与 C+E+G，返回逐 fold 的配对增量。"""
    raw_y = meta[target_col].to_numpy(float)
    valid = np.isfinite(raw_y) & (raw_y > 0)
    y = np.log10(raw_y[valid])
    groups = np.asarray(groups)[valid]
    base = np.column_stack([blocks["C"][valid], blocks["E"][valid]])
    full = np.column_stack([base, blocks["G"][valid]])
    cv = GroupKFold(n_splits=min(5, len(np.unique(groups))))
    rows = []
    for fold, (train, test) in enumerate(cv.split(base, y, groups)):
        model_args = dict(n_estimators=100, max_features=1.0, min_samples_leaf=2,
                          random_state=config.SEED, n_jobs=-1)
        m0 = RandomForestRegressor(**model_args).fit(base[train], y[train])
        m1 = RandomForestRegressor(**model_args).fit(full[train], y[train])
        from sklearn.metrics import r2_score
        r0 = float(r2_score(y[test], m0.predict(base[test])))
        r1 = float(r2_score(y[test], m1.predict(full[test])))
        rows.append({"target": target_col, "N": len(y), "fold": fold,
                     "r2_CE": r0, "r2_CEG": r1, "delta_r2_geometry": r1 - r0})
    return pd.DataFrame(rows)


def run_increment_only():
    """只重跑决定 Q1 的配对增量，供快速审计/复现。"""
    meta = pd.read_parquet(config.PROC_DIR / "views_meta.parquet")
    kappa = pd.read_parquet(config.PROC_DIR / "kappa_L_targets.parquet")
    exp = kappa[(kappa["method"] == "experimental") & kappa["material_id"].notna()]
    meta["kappa_exp"] = meta["material_id"].map(exp.set_index("material_id")["kappa_L"])
    meta, feats = load_aligned(
        {"soap_geo": config.PROC_DIR / "soap_geo.npy",
         "comp_frac": config.PROC_DIR / "comp_frac.npy"}, meta)
    meta["sys"] = meta["formula"].map(chemical_system)
    keep = meta["sys"].notna().to_numpy()
    meta = meta[keep].reset_index(drop=True)
    for name in feats:
        feats[name] = feats[name][keep]
    blocks, avg_mass = build_blocks(meta, feats)
    meta["avg_mass"] = avg_mass
    targets = ["snyder_acoustic", "clarke", "kappa_exp"]
    out = pd.concat([
        paired_geometry_increment(meta, blocks, t, meta["sys"].to_numpy()) for t in targets
    ], ignore_index=True)
    out.to_csv(config.PROC_DIR / "geometry_increment_folds.csv", index=False)
    summary = out.groupby(["target", "N"], as_index=False).agg(
        delta_mean=("delta_r2_geometry", "mean"),
        delta_std=("delta_r2_geometry", "std"),
        delta_min=("delta_r2_geometry", "min"),
        delta_max=("delta_r2_geometry", "max"))
    summary.to_csv(config.PROC_DIR / "geometry_increment_summary.csv", index=False)
    print(summary.to_string(index=False))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.4, 4.3))
    labels = summary["target"].tolist()
    ax.errorbar(labels, summary["delta_mean"], yerr=summary["delta_std"],
                fmt="o", capsize=5, color="#2a6fbb")
    ax.axhline(0, color="black", lw=1)
    ax.set_ylabel("paired ΔR² = R²(C+E+G) − R²(C+E)")
    ax.set_title("Incremental geometry signal (GroupKFold; mean ± fold SD)")
    for i, n in enumerate(summary["N"]):
        ax.text(i, summary.loc[i, "delta_mean"], f"  N={n}", va="bottom")
    plt.tight_layout(); plt.savefig(config.FIG_DIR / "geometry_increment.png", dpi=160)
    print("saved geometry_increment_folds.csv, geometry_increment_summary.csv, geometry_increment.png")


def block_ablation(meta, blocks, target_col, groups, use_X=False):
    """逐块消融：每个子集组合的 R²，返回组合 -> (mean, std)。"""
    raw_y = meta[target_col].values.astype(float)
    valid = np.isfinite(raw_y) & (raw_y > 0)
    y = np.log10(raw_y[valid])
    groups = np.asarray(groups)[valid]
    names = ["C", "G", "E"] + (["X"] if (use_X and "X" in blocks) else [])
    results = {}
    combos = [("C",), ("G",), ("E",), ("C", "E"), ("E", "G"), ("C", "E", "G")]
    if use_X and "X" in names:
        combos.extend([("X",), ("C", "E", "X"), ("C", "E", "G", "X")])
    for combo_names in combos:
        X = np.column_stack([blocks[n][valid] for n in combo_names])
        m, s = rf_cv(X, y, groups)
        results["+".join(combo_names)] = (m, s)
    return results


def partial_spearman(meta, target_col, features, confounders, avg_mass):
    """秩变量残差化的偏 Spearman，含截距。

    被控制变量本身的偏相关未定义，因而不再输出 density/avg_mass/debye 对自身的结果。
    """
    y = np.log10(meta[target_col].values.astype(float))
    conf_cols = [c for c in confounders if c in meta.columns]
    rows = []
    for c in features:
        if c in conf_cols:
            continue
        d = pd.DataFrame({"x": meta[c] if c in meta.columns else avg_mass,
                          "y": y})
        for cf in conf_cols:
            d[cf] = meta[cf].values
        d = d.dropna()
        if len(d) < 30:
            continue
        ranked = d[["x", "y", *conf_cols]].rank(method="average")
        Xc = np.column_stack([np.ones(len(ranked)), ranked[conf_cols].values])
        rx = ranked["x"].values - Xc @ np.linalg.lstsq(Xc, ranked["x"].values, rcond=None)[0]
        ry = ranked["y"].values - Xc @ np.linalg.lstsq(Xc, ranked["y"].values, rcond=None)[0]
        rho, pv = stats.pearsonr(rx, ry)
        rows.append({"feature": c, "N": len(d), "partial_spearman": round(float(rho), 3), "p": float(pv)})
    return pd.DataFrame(rows)


def residual_mi_proxy(meta, blocks, target_col, groups, n_perm=30):
    """交叉拟合后的非负残差 MI 诊断（不是严格的条件互信息估计量）。

    ``I(S;κ)-I(S;E)`` 不等于 ``I(S;κ|E)``；旧实现因此甚至给出 -0.50。
    这里分别用已知物理块 E 预测 y 与几何 PC，再估计两边残差的 MI，并用标签
    置换估计有限样本基线。
    """
    y = np.log10(meta[target_col].values.astype(float))
    # Structure 低维摘要：SOAP 前 5 个主成分
    from sklearn.decomposition import PCA
    S = blocks["G"]
    pca = PCA(n_components=5, random_state=config.SEED)
    S_pc = pca.fit_transform(S)
    E = blocks["E"]
    mask = np.isfinite(y) & (meta[target_col].values.astype(float) > 0)
    S_pc, E, y = S_pc[mask], E[mask], y[mask]
    groups = np.asarray(groups)[mask]
    cv = GroupKFold(n_splits=min(5, len(np.unique(groups))))
    rf_y = RandomForestRegressor(n_estimators=120, max_features=1.0,
                                 min_samples_leaf=3, random_state=config.SEED, n_jobs=-1)
    rf_s = RandomForestRegressor(n_estimators=120, max_features=1.0,
                                 min_samples_leaf=3, random_state=config.SEED + 1, n_jobs=-1)
    y_hat = cross_val_predict(rf_y, E, y, cv=cv, groups=groups, n_jobs=1)
    s_hat = cross_val_predict(rf_s, E, S_pc, cv=cv, groups=groups, n_jobs=1)
    y_res = y - y_hat
    s_res = S_pc - s_hat
    raw = float(mutual_info_regression(
        s_res, y_res, random_state=config.SEED, n_neighbors=5).sum())
    rng = np.random.RandomState(config.SEED)
    null = np.array([
        mutual_info_regression(s_res, rng.permutation(y_res),
                               random_state=config.SEED + i + 1, n_neighbors=5).sum()
        for i in range(n_perm)
    ])
    corrected = max(0.0, raw - float(null.mean()))
    return {
        "N": len(y),
        "residual_mi_raw": raw,
        "permutation_mean": float(null.mean()),
        "permutation_p95": float(np.quantile(null, 0.95)),
        "residual_mi_bias_corrected": corrected,
    }


def main():
    meta = pd.read_parquet(config.PROC_DIR / "views_meta.parquet")
    kappa = pd.read_parquet(config.PROC_DIR / "kappa_L_targets.parquet")
    exp = kappa[(kappa["method"] == "experimental") & kappa["material_id"].notna()]
    exp = exp.set_index("material_id")["kappa_L"]
    meta["kappa_exp"] = meta["material_id"].map(exp)

    meta, feats = load_aligned(
        {"soap_geo": config.PROC_DIR / "soap_geo.npy",
         "comp_frac": config.PROC_DIR / "comp_frac.npy"}, meta)
    meta["sys"] = meta["formula"].map(chemical_system)
    keep = meta["sys"].notna().to_numpy()
    meta = meta[keep].reset_index(drop=True)
    for name in feats:
        feats[name] = feats[name][keep]
    blocks, avg_mass = build_blocks(meta, feats)
    meta["avg_mass"] = avg_mass
    groups = meta["sys"].values

    print("=== 块消融（GroupKFold by chemical system, RandomForest R²）===")
    all_rows = []
    targets = ["snyder_acoustic", "clarke"]
    if meta["kappa_exp"].notna().sum() >= 30:
        targets.append("kappa_exp")
    for target in targets:
        use_X = "X" in blocks and meta["band_gap"].notna().sum() > 100
        res = block_ablation(meta, blocks, target, groups, use_X=use_X)
        print(f"  target=log10({target}):")
        for combo in sorted(res, key=lambda c: (c.count("+"), c)):
            m, s = res[combo]
            print(f"    {combo:12s} R²={m:+.3f} ± {s:.3f}")
            all_rows.append({"target": target, "combo": combo,
                             "N": int(meta[target].notna().sum()),
                             "r2_mean": round(m, 4), "r2_std": round(s, 4)})
    pd.DataFrame(all_rows).to_csv(config.PROC_DIR / "block_ablation.csv", index=False)

    print("=== 偏 Spearman（控制 density/avg_mass/debye）===")
    feats_list = ["bulk_vrh", "shear_vrh", "v_long", "v_trans", "nsites"]
    pc = partial_spearman(meta, "snyder_acoustic", feats_list, config.CONFOUNDERS, avg_mass)
    pc.to_csv(config.PROC_DIR / "partial_corr.csv", index=False)
    print(pc.to_string(index=False))

    print("=== 交叉拟合残差互信息（Structure 5-PC；非严格 CMI）===")
    cmi_rows = []
    for target in targets:
        c = residual_mi_proxy(meta, blocks, target, groups)
        c["target"] = target
        cmi_rows.append(c)
        print(f"  {target}: {c}")
    pd.DataFrame(cmi_rows).to_csv(config.PROC_DIR / "conditional_mi.csv", index=False)

    # 增量贡献条形图
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    df = pd.DataFrame(all_rows)
    fig, ax = plt.subplots(figsize=(7, 4.2))
    sub = df[df["target"] == "snyder_acoustic"]
    order = sorted(sub["combo"], key=lambda c: (c.count("+"), c))
    means = [sub[sub["combo"] == c]["r2_mean"].values[0] for c in order]
    stds = [sub[sub["combo"] == c]["r2_std"].values[0] for c in order]
    ax.bar(range(len(order)), means, yerr=stds, capsize=3)
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=45, ha="right")
    ax.set_ylabel("R² (GroupKFold CV)")
    ax.set_title("block ablation: log10(Snyder κ_L)")
    plt.tight_layout(); plt.savefig(config.FIG_DIR / "block_ablation.png", dpi=130)
    print("saved figures/block_ablation.png")


if __name__ == "__main__":
    main()
