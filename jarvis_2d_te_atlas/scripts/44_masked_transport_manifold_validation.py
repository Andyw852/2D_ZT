"""Masked-view validation for the three-view electronic-transport atlas.

Purpose
-------
Test whether composition/structure + electronic information can recover a
material's *withheld* electrical-transport performance, and whether a transport
manifold built only from the training materials adds useful information.

This is a retrospective discovery test.  For every fold, all transport
features of the test materials are removed before the transport graph is built.
The target is the JARVIS-defined log10 power factor.  No phonon, lattice thermal
conductivity, stability, or synthesizability data are used.

Outputs
-------
data/audit/masked_transport_manifold_metrics.csv
data/audit/masked_transport_manifold_summary.csv
data/processed/masked_transport_oof_predictions.csv
figures/masked_transport_manifold_validation.{png,pdf}
reports/masked_transport_manifold_validation.md
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.spatial.distance import cdist, pdist, squareform
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import average_precision_score, mean_absolute_error
from sklearn.model_selection import GroupKFold, KFold
from sklearn.preprocessing import RobustScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
from graph_utils import kNN_affinity  # noqa: E402
from multiview_utils import build_supra, joint_embedding  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
K_GRAPH = 15
K_REG = 15
EMBED_DIM = 16
RANDOM_SEED = 42
N_FOLDS = 5
N_GROUPS = 10
ANCHOR_LAMBDA = 0.3

MODEL_ORDER = [
    "structure",
    "electronic",
    "structure+electronic",
    "three-view (masked transport)",
    "electronic RF",
    "transport oracle [invalid]",
]

MODEL_LABELS = {
    "structure": "组成/结构流形",
    "electronic": "电子结构流形",
    "structure+electronic": "结构+电子流形",
    "three-view (masked transport)": "三视图流形（测试输运隐藏）",
    "electronic RF": "电子结构随机森林",
    "transport oracle [invalid]": "输运流形上限（泄漏对照）",
}

MODEL_PLOT_LABELS = {
    "structure": "Composition/structure",
    "electronic": "Electronic",
    "structure+electronic": "Structure + electronic",
    "three-view (masked transport)": "Three-view (masked T)",
    "electronic RF": "Electronic RF",
    "transport oracle [invalid]": "Transport oracle",
}


@dataclass
class CarrierData:
    carrier: str
    jids: list[str]
    formulas: list[str]
    d_structure: np.ndarray
    x_electronic: np.ndarray
    x_transport: np.ndarray
    y: np.ndarray
    pf: np.ndarray


def scale_layer(W: sparse.spmatrix) -> sparse.csr_matrix:
    """Normalize one graph layer to mean node strength one."""
    W = W.tocsr()
    mean_strength = float(np.asarray(W.sum(axis=1)).ravel().mean())
    return (W / mean_strength).tocsr() if mean_strength > 1e-12 else W


def affinity(D: np.ndarray, k: int = K_GRAPH) -> sparse.csr_matrix:
    """Deterministic local-scale kNN affinity."""
    return scale_layer(kNN_affinity(D, min(k, len(D) - 1), tiebreak_seed=RANDOM_SEED))


def diffusion_coordinates(W: sparse.spmatrix, n_components: int = EMBED_DIM) -> np.ndarray:
    """Diffusion-map coordinates from a symmetric affinity matrix."""
    W = W.tocsr()
    n = W.shape[0]
    degree = np.asarray(W.sum(axis=1)).ravel()
    dinv = sparse.diags(1.0 / np.sqrt(np.maximum(degree, 1e-12)))
    normalized = dinv @ W @ dinv
    k = min(n_components + 1, n - 1)
    vals, vecs = sparse.linalg.eigsh(normalized, k=k, which="LA")
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    # Drop the stationary direction. Eigenvalue weighting gives t=1 diffusion coordinates.
    return vecs[:, 1:k] * vals[None, 1:k]


def supra_consensus_coordinates(
    layers: list[tuple[str, sparse.spmatrix, list[str]]],
    consensus_jids: list[str],
    lam: float = ANCHOR_LAMBDA,
) -> np.ndarray:
    """Run the project's identity-anchored supra-graph alignment."""
    A, _, _ = build_supra(layers, consensus_jids, lam)
    _, vectors = joint_embedding(A, EMBED_DIM)
    if vectors is None:
        raise RuntimeError("Supra-graph eigensolver failed")
    # joint_embedding returns the trivial vector in column zero.
    return vectors[: len(consensus_jids), 1 : EMBED_DIM + 1]


def knn_regress(
    coords: np.ndarray,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    y_train: np.ndarray,
    k: int = K_REG,
) -> np.ndarray:
    D = cdist(coords[test_idx], coords[train_idx])
    kk = min(k, len(train_idx))
    nn = np.argpartition(D, kth=kk - 1, axis=1)[:, :kk]
    dn = np.take_along_axis(D, nn, axis=1)
    yn = y_train[nn]
    weights = 1.0 / np.maximum(dn, 1e-8)
    return (weights * yn).sum(axis=1) / weights.sum(axis=1)


def evaluate(y_true: np.ndarray, y_pred: np.ndarray, high_threshold: float) -> dict[str, float]:
    high = y_true >= high_threshold
    prevalence = float(high.mean())
    n_pick = max(1, int(np.ceil(0.20 * len(y_true))))
    selected = np.argsort(y_pred)[-n_pick:]
    precision = float(high[selected].mean())
    recall = float(high[selected].sum() / max(1, high.sum()))
    rho_raw = stats.spearmanr(y_true, y_pred).statistic
    # A constant extrapolation has no ranking skill; score it as zero rather
    # than silently dropping the fold from the aggregate.
    rho = 0.0 if not np.isfinite(rho_raw) else float(rho_raw)
    ap = float(average_precision_score(high.astype(int), y_pred)) if high.any() else np.nan
    return {
        "spearman": rho,
        "mae_log10_pf": float(mean_absolute_error(y_true, y_pred)),
        "average_precision": ap,
        "precision_at_20pct": precision,
        "recall_at_20pct": recall,
        "high_pf_prevalence": prevalence,
        "precision_lift": precision / prevalence if prevalence > 0 else np.nan,
    }


def load_carrier(carrier: str) -> CarrierData:
    s_df = pd.read_parquet(ROOT / "features/structure/geometry_soap_v1.parquet").sort_values("jid")
    s_jids = s_df["jid"].tolist()
    s_index = {jid: i for i, jid in enumerate(s_jids)}
    d_structure_all = np.load(ROOT / "data/processed/d_struct_baseline.npy")

    e_df = pd.read_parquet(ROOT / "features/electronic/electronic_features_v1.parquet").set_index("jid")
    t_df = pd.read_parquet(
        ROOT / f"features/transport/{carrier}_transport_features_v1.parquet"
    ).set_index("jid")
    structures = pd.read_parquet(ROOT / "data/processed/standardized_2d_structures.parquet")
    formula = structures.set_index("jid")["formula"].to_dict()

    mass_col = "m_elec_median" if carrier == "n" else "m_hole_median"
    common = sorted(
        set(t_df.index)
        & set(e_df.dropna(subset=["Eg_optb88vdw", mass_col]).index)
        & set(s_jids)
    )
    si = [s_index[jid] for jid in common]
    d_structure = d_structure_all[np.ix_(si, si)]

    x_electronic = e_df.loc[common, ["Eg_optb88vdw", mass_col]].to_numpy(float)
    x_electronic = RobustScaler().fit_transform(x_electronic)
    t_cols = [
        "S_median",
        "S_MAD",
        "S_sign_fraction",
        "log_sigma_dom_geo",
        "D_sigma",
        "A_sigma_dom",
    ]
    x_transport = t_df.loc[common, t_cols].to_numpy(float)
    pf = t_df.loc[common, "PF_mean"].to_numpy(float)
    y = np.log10(np.maximum(pf, 1e-3))
    return CarrierData(
        carrier=carrier,
        jids=common,
        formulas=[formula.get(jid, "") for jid in common],
        d_structure=d_structure,
        x_electronic=x_electronic,
        x_transport=x_transport,
        y=y,
        pf=pf,
    )


def build_fixed_manifolds(data: CarrierData) -> tuple[dict[str, sparse.csr_matrix], dict[str, np.ndarray]]:
    w_structure = affinity(data.d_structure)
    w_gap = affinity(squareform(pdist(data.x_electronic[:, [0]])))
    w_mass = affinity(squareform(pdist(data.x_electronic[:, [1]])))
    jids = data.jids
    graphs = {
        "structure": w_structure,
        "gap": w_gap,
        "mass": w_mass,
    }
    coords = {
        "structure": diffusion_coordinates(w_structure),
        "electronic": supra_consensus_coordinates(
            [("gap", w_gap, jids), ("mass", w_mass, jids)], jids
        ),
        "structure+electronic": supra_consensus_coordinates(
            [
                ("structure", w_structure, jids),
                ("gap", w_gap, jids),
                ("mass", w_mass, jids),
            ],
            jids,
        ),
    }
    return graphs, coords


def transport_graph_train(data: CarrierData, train_idx: np.ndarray) -> sparse.csr_matrix:
    """Build a transport graph with test materials completely absent."""
    scaler = RobustScaler().fit(data.x_transport[train_idx])
    X = scaler.transform(data.x_transport[train_idx])
    W_train = affinity(squareform(pdist(X)))
    W_train = W_train.tocoo()
    rows, cols = W_train.row, W_train.col
    W = sparse.coo_matrix(
        (W_train.data, (train_idx[rows], train_idx[cols])),
        shape=(len(data.jids), len(data.jids)),
    )
    return W.tocsr()


def transport_graph_oracle(data: CarrierData) -> sparse.csr_matrix:
    X = RobustScaler().fit_transform(data.x_transport)
    return affinity(squareform(pdist(X)))


def make_splits(data: CarrierData, se_coords: np.ndarray) -> dict[str, list[tuple[np.ndarray, np.ndarray]]]:
    random_cv = list(
        KFold(n_splits=N_FOLDS, shuffle=True, random_state=RANDOM_SEED).split(data.y)
    )
    groups = KMeans(n_clusters=N_GROUPS, n_init=30, random_state=RANDOM_SEED).fit_predict(se_coords)
    blocked_cv = list(GroupKFold(n_splits=N_FOLDS).split(data.y, groups=groups))
    return {"random": random_cv, "manifold-blocked": blocked_cv}


def run_carrier(data: CarrierData) -> tuple[pd.DataFrame, pd.DataFrame]:
    fixed_graphs, fixed_coords = build_fixed_manifolds(data)
    fixed_layers = [
        ("structure", fixed_graphs["structure"], data.jids),
        ("gap", fixed_graphs["gap"], data.jids),
        ("mass", fixed_graphs["mass"], data.jids),
    ]
    oracle_coords = supra_consensus_coordinates(
        fixed_layers + [("transport", transport_graph_oracle(data), data.jids)],
        data.jids,
    )
    splits = make_splits(data, fixed_coords["structure+electronic"])

    rows: list[dict] = []
    oof_rows: list[dict] = []
    for split_name, folds in splits.items():
        for fold, (train_idx, test_idx) in enumerate(folds, start=1):
            y_train, y_test = data.y[train_idx], data.y[test_idx]
            high_threshold = float(np.quantile(y_train, 0.80))

            predictions: dict[str, np.ndarray] = {}
            for model in ["structure", "electronic", "structure+electronic"]:
                predictions[model] = knn_regress(
                    fixed_coords[model], train_idx, test_idx, y_train
                )

            w_t_masked = transport_graph_train(data, train_idx)
            train_jids = [data.jids[i] for i in train_idx]
            masked_coords = supra_consensus_coordinates(
                fixed_layers
                + [("transport", w_t_masked[train_idx][:, train_idx], train_jids)],
                data.jids,
            )
            predictions["three-view (masked transport)"] = knn_regress(
                masked_coords, train_idx, test_idx, y_train
            )

            rf = RandomForestRegressor(
                n_estimators=400,
                min_samples_leaf=5,
                max_features=1.0,
                random_state=RANDOM_SEED + fold,
                n_jobs=-1,
            )
            rf.fit(data.x_electronic[train_idx], y_train)
            predictions["electronic RF"] = rf.predict(data.x_electronic[test_idx])

            # Deliberately invalid upper-bound control: test transport is present in this graph.
            predictions["transport oracle [invalid]"] = knn_regress(
                oracle_coords, train_idx, test_idx, y_train
            )

            for model, y_pred in predictions.items():
                metric = evaluate(y_test, y_pred, high_threshold)
                rows.append(
                    {
                        "carrier": data.carrier,
                        "split": split_name,
                        "fold": fold,
                        "model": model,
                        "n_train": len(train_idx),
                        "n_test": len(test_idx),
                        "high_pf_threshold_log10": high_threshold,
                        **metric,
                    }
                )
                if split_name == "random":
                    for local_i, pred in zip(test_idx, y_pred):
                        oof_rows.append(
                            {
                                "carrier": data.carrier,
                                "jid": data.jids[local_i],
                                "formula": data.formulas[local_i],
                                "fold": fold,
                                "model": model,
                                "observed_log10_pf": data.y[local_i],
                                "observed_pf": data.pf[local_i],
                                "predicted_log10_pf": float(pred),
                            }
                        )
    return pd.DataFrame(rows), pd.DataFrame(oof_rows)


def summarize(metrics: pd.DataFrame) -> pd.DataFrame:
    measures = [
        "spearman",
        "mae_log10_pf",
        "average_precision",
        "precision_at_20pct",
        "recall_at_20pct",
        "high_pf_prevalence",
        "precision_lift",
    ]
    summary = metrics.groupby(["carrier", "split", "model"])[measures].agg(["mean", "std"])
    summary.columns = [f"{a}_{b}" for a, b in summary.columns]
    return summary.reset_index()


def make_figure(summary: pd.DataFrame, oof: pd.DataFrame) -> None:
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9})
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.5), constrained_layout=True)
    valid = [m for m in MODEL_ORDER if "oracle" not in m]
    colors = plt.cm.Set2(np.linspace(0, 1, len(valid)))

    for ax, split_name, title in [
        (axes[0, 0], "random", "(a) Random masked-view validation"),
        (axes[0, 1], "manifold-blocked", "(b) Manifold-blocked extrapolation"),
    ]:
        x = np.arange(len(valid))
        width = 0.36
        for ci, carrier in enumerate(["n", "p"]):
            sub = summary[(summary.carrier == carrier) & (summary.split == split_name)].set_index("model")
            vals = [sub.loc[m, "spearman_mean"] for m in valid]
            errs = [sub.loc[m, "spearman_std"] for m in valid]
            ax.bar(x + (ci - 0.5) * width, vals, width, yerr=errs, capsize=2, label=f"{carrier}-type", alpha=0.88)
        ax.set_xticks(x, [MODEL_PLOT_LABELS[m] for m in valid], rotation=22, ha="right")
        ax.set_ylabel("Spearman ρ: predicted vs observed log PF")
        ax.set_title(title, loc="left")
        ax.axhline(0, color="0.45", lw=0.8)
        ax.legend(frameon=False)

    ax = axes[1, 0]
    x = np.arange(len(valid))
    width = 0.36
    for ci, carrier in enumerate(["n", "p"]):
        sub = summary[(summary.carrier == carrier) & (summary.split == "random")].set_index("model")
        vals = [sub.loc[m, "precision_lift_mean"] for m in valid]
        errs = [sub.loc[m, "precision_lift_std"] for m in valid]
        ax.bar(x + (ci - 0.5) * width, vals, width, yerr=errs, capsize=2, label=f"{carrier}-type", alpha=0.88)
    ax.set_xticks(x, [MODEL_PLOT_LABELS[m] for m in valid], rotation=22, ha="right")
    ax.axhline(1, color="0.45", lw=0.8, ls="--")
    ax.set_ylabel("Top-20% precision lift over prevalence")
    ax.set_title("(c) High-PF enrichment without test transport", loc="left")

    ax = axes[1, 1]
    model = "three-view (masked transport)"
    for carrier, marker in [("n", "o"), ("p", "^")]:
        sub = oof[(oof.carrier == carrier) & (oof.model == model)]
        ax.scatter(
            sub.observed_log10_pf,
            sub.predicted_log10_pf,
            s=14,
            alpha=0.45,
            marker=marker,
            label=f"{carrier}-type",
        )
    lims = [
        min(oof.observed_log10_pf.min(), oof.predicted_log10_pf.min()),
        max(oof.observed_log10_pf.max(), oof.predicted_log10_pf.max()),
    ]
    ax.plot(lims, lims, color="0.35", lw=0.9, ls="--")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Observed log10(PF), JARVIS convention")
    ax.set_ylabel("Out-of-fold predicted log10(PF)")
    ax.set_title("(d) Honest out-of-fold predictions", loc="left")
    ax.legend(frameon=False)

    fig.suptitle(
        "Three-view electronic-transport manifold: masked-view discovery test\n"
        "Test-material transport is hidden; no phonon, κL, stability, or synthesis variables",
        fontsize=12,
    )
    fig.savefig(ROOT / "figures/masked_transport_manifold_validation.png", dpi=220)
    fig.savefig(ROOT / "figures/masked_transport_manifold_validation.pdf")
    plt.close(fig)


def fmt(value: float) -> str:
    return f"{value:.3f}" if pd.notna(value) else "NA"


def write_report(summary: pd.DataFrame, oof: pd.DataFrame, data_sizes: dict[str, int]) -> None:
    lines = [
        "# 三视图流形的隐藏输运验证",
        "",
        "> 研究边界：仅使用组成/结构、电子结构和 n/p 电输运。目标是数据库定义的 PF，",
        "> 不使用声子、晶格热导率、稳定性或可合成性，也不宣称预测真实 ZT。",
        "",
        "## 验证问题",
        "",
        "对每个测试材料，先完全隐藏其 6 个输运特征和 PF，再用其组成/结构与电子结构将其放入材料流形。",
        "训练材料的输运图可以作为第三视图，但测试材料不能出现在输运图中。这样模拟“已知材料有输运，",
        "候选材料只有结构和电子信息”的发现情景。",
        "",
        f"- n 型完整三视图交集：{data_sizes['n']} 个材料。",
        f"- p 型完整三视图交集：{data_sizes['p']} 个材料。",
        "- 指标：log10(PF) 的 Spearman 相关、MAE，以及实际高 PF（训练集 top 20% 阈值）的富集倍数。",
        "- random：普通 5 折隐藏；manifold-blocked：将结构+电子流形分区后整区留出，测试跨区域外推。",
        "",
        "## 核心结果",
        "",
        "| carrier | split | model | Spearman ρ | MAE logPF | high-PF lift |",
        "|---|---|---|---:|---:|---:|",
    ]
    for carrier in ["n", "p"]:
        for split_name in ["random", "manifold-blocked"]:
            sub = summary[(summary.carrier == carrier) & (summary.split == split_name)].set_index("model")
            for model in MODEL_ORDER:
                r = sub.loc[model]
                lines.append(
                    f"| {carrier} | {split_name} | {MODEL_LABELS[model]} | "
                    f"{fmt(r['spearman_mean'])}±{fmt(r['spearman_std'])} | "
                    f"{fmt(r['mae_log10_pf_mean'])} | {fmt(r['precision_lift_mean'])}× |"
                )

    lines += ["", "## 结论判据", ""]
    for carrier in ["n", "p"]:
        rnd = summary[(summary.carrier == carrier) & (summary.split == "random")].set_index("model")
        blk = summary[(summary.carrier == carrier) & (summary.split == "manifold-blocked")].set_index("model")
        m = "three-view (masked transport)"
        se = "structure+electronic"
        delta_rnd = rnd.loc[m, "spearman_mean"] - rnd.loc[se, "spearman_mean"]
        delta_blk = blk.loc[m, "spearman_mean"] - blk.loc[se, "spearman_mean"]
        lines.append(
            f"- {carrier} 型：加入仅含训练材料的输运视图后，Spearman 相对结构+电子流形变化 "
            f"{delta_rnd:+.3f}（random）/{delta_blk:+.3f}（blocked）。"
        )

    lines += [
        "",
        "输运 oracle 把测试材料的输运特征放回图中，属于有意设置的泄漏上限，只用于显示“看过输运以后”",
        "图谱能达到的结果，不能当作发现性能。有效结论必须只看前五个模型。",
        "",
        "## 隐藏后找回的高 PF 材料",
        "",
        "以下排序来自 random 5-fold 的 out-of-fold 预测；每个材料被预测时，其输运特征均未参与建图。",
        "",
        "| carrier | jid | formula | predicted logPF | observed PF | observed percentile |",
        "|---|---|---|---:|---:|---:|",
    ]
    model = "three-view (masked transport)"
    for carrier in ["n", "p"]:
        sub = oof[(oof.carrier == carrier) & (oof.model == model)].copy()
        sub["pct"] = sub.observed_pf.rank(pct=True) * 100
        pred_cut = sub.predicted_log10_pf.quantile(0.80)
        obs_cut = sub.observed_pf.quantile(0.80)
        selected = sub[sub.predicted_log10_pf >= pred_cut]
        hits = selected[selected.observed_pf >= obs_cut]
        precision = len(hits) / len(selected)
        recall = len(hits) / max(1, int((sub.observed_pf >= obs_cut).sum()))
        lines.append(
            f"| {carrier} | **OOF 汇总** | predicted top20% 中命中 {len(hits)}/{len(selected)} | "
            f"precision={precision:.3f} | recall={recall:.3f} | lift={precision/0.2:.2f}× |"
        )
        for _, r in hits.nlargest(12, "predicted_log10_pf").iterrows():
            lines.append(
                f"| {carrier} | {r.jid} | {r.formula} | {r.predicted_log10_pf:.3f} | "
                f"{r.observed_pf:.2f} | {r.pct:.1f}% |"
            )

    lines += [
        "",
        "## 可以与不可以声称的结果",
        "",
        "可以声称：现有三视图是否能够在未知测试输运的条件下富集高 PF 材料，以及这种能力在跨流形区域",
        "外推时是否保留。",
        "",
        "不可以声称：预测真实 ZT；PF 是电输运性能指标，而且 JARVIS PF 还受本征值配对约定影响。",
        "",
    ]
    (ROOT / "reports/masked_transport_manifold_validation.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )


def main() -> None:
    all_metrics = []
    all_oof = []
    sizes = {}
    for carrier in ["n", "p"]:
        data = load_carrier(carrier)
        sizes[carrier] = len(data.jids)
        print(f"{carrier}-type: N={len(data.jids)}")
        metrics, oof = run_carrier(data)
        all_metrics.append(metrics)
        all_oof.append(oof)

    metrics = pd.concat(all_metrics, ignore_index=True)
    oof = pd.concat(all_oof, ignore_index=True)
    summary = summarize(metrics)

    metrics.to_csv(ROOT / "data/audit/masked_transport_manifold_metrics.csv", index=False)
    summary.to_csv(ROOT / "data/audit/masked_transport_manifold_summary.csv", index=False)
    oof.to_csv(ROOT / "data/processed/masked_transport_oof_predictions.csv", index=False)
    make_figure(summary, oof)
    write_report(summary, oof, sizes)

    show = summary[
        summary.model.isin(["structure+electronic", "three-view (masked transport)", "electronic RF"])
    ][["carrier", "split", "model", "spearman_mean", "mae_log10_pf_mean", "precision_lift_mean"]]
    print("\n", show.to_string(index=False))
    print("\nWrote masked-view metrics, OOF predictions, report, and figure.")


if __name__ == "__main__":
    main()
