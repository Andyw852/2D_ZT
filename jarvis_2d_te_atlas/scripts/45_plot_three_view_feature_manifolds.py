"""Plot cluster and frozen-feature maps in the full three-view manifold.

The n/p full manifolds use the cached lambda=0.3 supra-graph solutions with
composition/structure, band gap, carrier-specific effective mass, and the six
electrical-transport features.  Every panel reuses exactly the same Phi_1/Phi_2
coordinates within a carrier so that color patterns are directly comparable.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


ROOT = Path(__file__).resolve().parents[1]
SEED = 42
N_CLUSTER_DIMS = 10

FEATURE_SPECS = [
    ("cluster", "Cluster", None, None),
    ("Eg_optb88vdw", "Band gap Eg", "eV", "viridis"),
    ("carrier_mass", "Carrier effective mass", "m_e", "plasma"),
    ("S_median", "S median", "uV/K", "coolwarm"),
    ("logPF", "log10(PF)", "JARVIS convention", "magma"),
    ("S_MAD", "S MAD", "uV/K", "viridis"),
    ("S_sign_fraction", "Sign consistency", "fraction", "viridis"),
    ("log_sigma_dom_geo", "log sigma dominant", "log scale", "viridis"),
    ("D_sigma", "Suppressed-channel contrast", "log ratio", "cividis"),
    ("A_sigma_dom", "Dominant-channel anisotropy", "log ratio", "cividis"),
]


def load_full_space(carrier: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    cache = np.load(ROOT / "data/processed/joint_coords_cache.npz", allow_pickle=True)
    coords = cache[f"{carrier}_Full_0.3"]
    structure_nodes = pd.read_csv(ROOT / "graphs/G_structure_v1_nodes.csv")
    jids = sorted(structure_nodes.jid.tolist())
    consensus = coords[: len(jids), :20]

    df = pd.DataFrame(consensus, columns=[f"Phi_{i}" for i in range(1, 21)])
    df.insert(0, "jid", jids)

    structures = pd.read_parquet(ROOT / "data/processed/standardized_2d_structures.parquet")
    formulas = structures.set_index("jid")["formula"]
    electronic = pd.read_parquet(ROOT / "features/electronic/electronic_features_v1.parquet").set_index("jid")
    transport = pd.read_parquet(
        ROOT / f"features/transport/{carrier}_transport_features_v1.parquet"
    ).set_index("jid")
    mass_col = "m_elec_median" if carrier == "n" else "m_hole_median"
    tcols = [
        "S_median",
        "S_MAD",
        "S_sign_fraction",
        "log_sigma_dom_geo",
        "D_sigma",
        "A_sigma_dom",
        "PF_mean",
    ]
    df = df.join(formulas.rename("formula"), on="jid")
    df = df.join(electronic[["Eg_optb88vdw", mass_col]], on="jid")
    df = df.rename(columns={mass_col: "carrier_mass"})
    df = df.join(transport[tcols], on="jid")
    df["logPF"] = np.log10(np.maximum(df["PF_mean"], 1e-3))

    X = df[[f"Phi_{i}" for i in range(1, N_CLUSTER_DIMS + 1)]].to_numpy()
    scores = []
    fits = {}
    for k in range(3, 7):
        km = KMeans(n_clusters=k, n_init=50, random_state=SEED).fit(X)
        score = silhouette_score(X, km.labels_)
        scores.append({"carrier": carrier, "k": k, "silhouette": score})
        fits[k] = km
    best_k = max(scores, key=lambda row: row["silhouette"])["k"]
    labels = fits[best_k].labels_

    # Stable numbering from left-to-right cluster centroid in the display plane.
    means = pd.DataFrame({"label": labels, "x": df.Phi_1, "y": df.Phi_2}).groupby("label").mean()
    order = means.sort_values(["x", "y"]).index.tolist()
    remap = {old: new + 1 for new, old in enumerate(order)}
    df["cluster"] = [remap[x] for x in labels]

    profile_cols = [x[0] for x in FEATURE_SPECS[1:]]
    profiles = df.groupby("cluster")[profile_cols].median(numeric_only=True)
    profiles.insert(0, "n_materials", df.groupby("cluster").size())
    profiles = profiles.reset_index()
    profiles.insert(0, "carrier", carrier)
    profiles["selected_k"] = best_k
    profiles["silhouette"] = max(x["silhouette"] for x in scores)
    return df, pd.DataFrame(scores), profiles


def robust_norm(values: pd.Series) -> Normalize:
    x = values.dropna().to_numpy(float)
    if len(x) == 0:
        return Normalize(0, 1)
    lo, hi = np.quantile(x, [0.02, 0.98])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(np.nanmin(x)), float(np.nanmax(x) + 1e-12)
    return Normalize(lo, hi, clip=True)


def draw_panel(ax, df: pd.DataFrame, spec, carrier: str, panel_index: int) -> None:
    col, label, unit, cmap = spec
    ax.scatter(df.Phi_1, df.Phi_2, s=5, c="0.82", alpha=0.34, linewidths=0, rasterized=True)
    if col == "cluster":
        k = int(df.cluster.max())
        palette = plt.get_cmap("tab10")
        for cluster in range(1, k + 1):
            sub = df[df.cluster == cluster]
            ax.scatter(
                sub.Phi_1,
                sub.Phi_2,
                s=8,
                color=palette(cluster - 1),
                alpha=0.72,
                linewidths=0,
                rasterized=True,
            )
            ax.text(
                sub.Phi_1.median(),
                sub.Phi_2.median(),
                f"C{cluster}",
                ha="center",
                va="center",
                fontsize=8,
                fontweight="bold",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.2},
            )
        title = f"{carrier}-type: clusters (k={k})"
    else:
        sub = df.dropna(subset=[col])
        norm = robust_norm(sub[col])
        marks = ax.scatter(
            sub.Phi_1,
            sub.Phi_2,
            s=8,
            c=sub[col],
            cmap=cmap,
            norm=norm,
            alpha=0.78,
            linewidths=0,
            rasterized=True,
        )
        cb = ax.figure.colorbar(marks, ax=ax, fraction=0.045, pad=0.02)
        cb.ax.tick_params(labelsize=7, length=2)
        title = f"{carrier}-type: {label} ({unit}), N={len(sub)}"
    ax.set_title(title, loc="left", fontsize=9)
    ax.set_xlabel("Phi 1")
    ax.set_ylabel("Phi 2")
    ax.tick_params(labelsize=7, length=2)
    ax.text(
        0.01,
        0.99,
        chr(ord("a") + panel_index),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )


def make_combined_figure(frames: dict[str, pd.DataFrame]) -> None:
    fig, axes = plt.subplots(4, 5, figsize=(18, 14), constrained_layout=True)
    for row_block, carrier in enumerate(["n", "p"]):
        df = frames[carrier]
        for i, spec in enumerate(FEATURE_SPECS):
            row = row_block * 2 + i // 5
            col = i % 5
            draw_panel(axes[row, col], df, spec, carrier, row * 5 + col)
    fig.suptitle(
        "Full three-view manifold: clusters and frozen physical features\n"
        "Same point coordinates within each carrier; lambda=0.3; gray = feature missing",
        fontsize=14,
    )
    fig.savefig(ROOT / "figures/three_view_feature_manifolds.png", dpi=240)
    fig.savefig(ROOT / "figures/three_view_feature_manifolds.pdf")
    plt.close(fig)


def write_report(scores: pd.DataFrame, profiles: pd.DataFrame) -> None:
    lines = [
        "# 三视图空间聚类与逐特征流形图",
        "",
        "- 坐标：包含组成/结构、电子结构和载流子对应输运视图的 Full supra-graph consensus，λ=0.3。",
        "- n/p 型分别建空间；每个载流子内部所有面板使用完全相同的 Phi_1/Phi_2 坐标。",
        "- 聚类：在前 10 个流形坐标上对 k=3–6 扫描 KMeans，以 silhouette 最大值选 k。",
        "- 灰点表示材料存在于公共空间，但该面板对应的性质缺失。连续色标截断在 2%–98% 分位，避免异常值压扁颜色。",
        "- PF 仅用于着色，不参与独立发现能力证明；这张图是描述性图谱。",
        "",
        "## 聚类选择",
        "",
        "| carrier | k | silhouette | selected |",
        "|---|---:|---:|---|",
    ]
    for carrier in ["n", "p"]:
        sub = scores[scores.carrier == carrier]
        best = int(sub.loc[sub.silhouette.idxmax(), "k"])
        for _, row in sub.iterrows():
            lines.append(
                f"| {carrier} | {int(row.k)} | {row.silhouette:.4f} | {'yes' if int(row.k)==best else ''} |"
            )
    lines += ["", "## 聚类中位特征", ""]
    for carrier in ["n", "p"]:
        lines += [f"### {carrier} 型", "", profiles[profiles.carrier == carrier].to_markdown(index=False), ""]
    (ROOT / "reports/three_view_feature_manifolds.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    frames = {}
    score_parts = []
    profile_parts = []
    json_rows = []
    for carrier in ["n", "p"]:
        df, scores, profiles = load_full_space(carrier)
        frames[carrier] = df
        score_parts.append(scores)
        profile_parts.append(profiles)
        keep = ["jid", "formula", "Phi_1", "Phi_2"] + [x[0] for x in FEATURE_SPECS]
        out = df[keep].copy()
        out.insert(0, "carrier", carrier)
        json_rows.extend(out.where(pd.notna(out), None).to_dict(orient="records"))

    scores = pd.concat(score_parts, ignore_index=True)
    profiles = pd.concat(profile_parts, ignore_index=True)
    points = pd.concat(
        [frames[c].assign(carrier=c) for c in ["n", "p"]], ignore_index=True
    )
    points.to_parquet(ROOT / "manifolds/three_view_feature_manifold_points.parquet", index=False)
    scores.to_csv(ROOT / "data/audit/three_view_cluster_selection.csv", index=False)
    profiles.to_csv(ROOT / "data/processed/three_view_cluster_profiles.csv", index=False)
    (ROOT / "data/processed/three_view_feature_manifold_points.json").write_text(
        json.dumps(json_rows, ensure_ascii=False, separators=(",", ":")), encoding="utf-8"
    )
    make_combined_figure(frames)
    write_report(scores, profiles)
    print(scores.to_string(index=False))
    print("\nCluster profiles:\n", profiles.to_string(index=False))
    print("\nWrote combined feature-manifold figure, points, profiles, and report.")


if __name__ == "__main__":
    main()
